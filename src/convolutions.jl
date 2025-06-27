# Direct convolution is faster than FFT-based convolution when the
# kernel is small. This is the case for the convolution to simulate
# instrumental resolution since the instrumental sampling is typically
# ~2 times netter than the instrument resolution, resulting in kernels
# being 21 points long (for a 5sigma kernel)

# See https://www.nv5geospatialsoftware.com/docs/CONVOL.html for "edge" mode meaning
_edge_left( x::Vector{<: Real}, i::Int, value::Real) = value
_edge_right(x::Vector{<: Real}, i::Int, value::Real) = value

_edge_left( x::Vector{T}, i::Int, mode::Val{:edge_zero}) where T <: Real = zero(T)
_edge_right(x::Vector{T}, i::Int, mode::Val{:edge_zero}) where T <: Real = zero(T)

_edge_left( x::Vector{<: Real}, i::Int, mode::Val{:edge_extend}) = x[1]    # this corresponds to EDGE_TRUNCATE in IDL
_edge_right(x::Vector{<: Real}, i::Int, mode::Val{:edge_extend}) = x[end]  # this corresponds to EDGE_TRUNCATE in IDL

_edge_left( x::Vector{<: Real}, i::Int, mode::Val{:edge_mirror}) = x[1 - i]
_edge_right(x::Vector{<: Real}, i::Int, mode::Val{:edge_mirror}) = x[end-(i - length(x) - 1)]

_edge_left( x::Vector{<: Real}, i::Int, mode::Val{:edge_reflect}) = x[2 - i]
_edge_right(x::Vector{<: Real}, i::Int, mode::Val{:edge_reflect}) = x[end-(i - length(x))]


function direct_conv1d!(z::Vector{T}, x::Vector{T}, y::Vector{T}, edge_mode=Val(:edge_zero)) where T <: Real
    if length(x) < length(y)
        return direct_conv1d!(z, y, x, edge_mode)
    end
    @assert length(x) == length(z)
    @assert isodd(length(y))

    # z(t) = \int x(t - τ) y(τ) dτ
    fill!(z, zero(T))
    N = div(length(y) - 1, 2)
	for τ in -N:N
		@inbounds vy = y[τ + N + 1]
        @batch_when_threaded for t in intersect(eachindex(x), eachindex(x) .+ τ)
            # @printf "%4d  %4d\n" t t-τ
			@inbounds z[t] += x[t - τ] * vy
        end
    end

    for τ in 1:N
		@inbounds vy = y[τ + N + 1]
        for t in 1:τ # setdiff(eachindex(x), eachindex(x) .+ τ)
            # @printf "left  %4d  %4d\n" t t-τ
			@inbounds z[t] += _edge_left(x, t - τ, edge_mode) * vy
        end
    end

	for τ in -N:-1
		@inbounds vy = y[τ + N + 1]
        for t in (length(z) + τ + 1):length(z) # setdiff(eachindex(x), eachindex(x) .+ τ)
            # @printf "right %4d  %4d\n" t t-τ
			@inbounds z[t] += _edge_right(x, t - τ, edge_mode) * vy
        end
    end

    return z
end

direct_conv1d(x::Vector{T}, y::Vector{T}, edge_mode=Val(:edge_zero)) where T = direct_conv1d!(similar(x, length(x)), x, y, edge_mode)



function direct_conv1d!(z::Vector{T}, x::Vector{T}, y::Matrix{T}) where T
    @assert length(x) >  size(y)[1]
    @assert length(x) == size(y)[2]
    @assert length(x) == length(z)
    @assert isodd(size(y)[1])
    fill!(z, zero(T))
    N = div(size(y)[1] - 1, 2)
	for τ in -N:N
        @batch_when_threaded for t in intersect(eachindex(x), eachindex(x) .+ τ)
            # @printf "%4d %4d %4d -> %4d\n" τ t-τ t
		    @inbounds vy = y[τ + N + 1, t - τ]
			@inbounds z[t] += x[t - τ] * vy
        end
    end
    return z
end



function gauss_broadening(x, y, σ_kms)
    Rsampling = 2 * maximum(sampling_resolutions(x))
    grid_x, grid_y = interpolate_on_logregular_grid(x, y, Rsampling)
    kernel = gauss_kernel(Rsampling, 3e5/σ_kms)
    @assert isodd(length(kernel))
    convolved = direct_conv1d(grid_y, kernel)
    # Interpolate back to original domain
    return Dierckx.Spline1D(grid_x, convolved, k=1, bc="extrapolate")(x)
end
