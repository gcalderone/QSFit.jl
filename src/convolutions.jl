# Compute direct convolution in 1D
# z(t) = \int x(t - τ) y(τ) dτ
function direct_conv1d!(z::Vector{T}, x::Vector{T}, y::Vector{T}) where T
    @assert length(x) >  length(y)
    @assert length(x) == length(z)
    @assert isodd(length(y))
    fill!(z, zero(T))
    N = div(length(y) - 1, 2)
	for τ in -N:N
		@inbounds vy = y[τ + N + 1]
        @batch_when_threaded for t in intersect(eachindex(x), eachindex(x) .+ τ)
            # @printf "%4d %4d %4d -> %4d\n" τ t-τ t
			@inbounds z[t] += x[t - τ] * vy
        end
    end
    return z
end


function direct_conv1d(x::Vector{T}, y::Vector{T}) where T
    if length(x) < length(y)
        return conv1d(y, x)
    end
    @assert isodd(length(y))
    z = similar(x, length(x))
    direct_conv1d!(z, x, y)
    return z
end


function gauss_broadening(x, y, σ_kms)
    Rsampling = 2 * maximum(sampling_resolutions(x))
    grid_x, grid_y = interpolate_on_logregular_grid(x, y, Rsampling)
    kernel = gauss_kernel(Rsampling, σ_kms)
    @assert isodd(length(kernel))
    convolved = direct_conv1d(grid_y, kernel)
    # Interpolate back to original domain
    return Dierckx.Spline1D(grid_x, convolved, k=1, bc="extrapolate")(x)
end
