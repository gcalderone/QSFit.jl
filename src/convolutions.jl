#=
A regular log-λ grid is characterized by:
log10(λ_i+1) - log10(λ_i)  =  log10(λ_i+1 / λ_i)  =  costant step

hence:
λ_i+1 / λ_i  =  cost.
λ_i+1 / λ_i - 1  =  cost.
(λ_i+1 - λ_i) / λ_i  =  Δλ / λ  =  1/R   is constant ∀i

i.e. it has a constant spectral resolution given by:
R = c / σ_kms

The step of the grid is:
log10(σ_kms / c + 1)
=#

function conv_gauss(λ, data, σ_kms; oversampling::Int=1)
    @assert oversampling >= 1

    # Interpolate input data on the regular log-λ grid with proper resolution
    # check with: extrema([(grid_λ[i+1] - grid_λ[i]) / grid_λ[i] * 3e5 for i in 1:length(grid_λ)-1]) .* oversampling, σ_kms
    step = log10(σ_kms / 3e5 + 1) / oversampling
    grid_λ = 10. .^ range(log10(minimum(λ)), log10(maximum(λ)), step=step)
    grid_y = Spline1D(λ, data, k=1)(grid_λ)

    # Convolution with Gaussian kernel and extraction of proper subset
    kernel = gauss.(-5:1. / oversampling:5, 0., 1.)
    i1 = div(length(kernel), 2) + 1
    i2 = i1 + length(grid_y) - 1
    out = conv(grid_y, kernel)[i1:i2] ./ oversampling

    # Interpolate back to original domain
    out = Spline1D(grid_λ, out, k=1, bc="extrapolate")(λ)

    # Replace data close to the edges
    ee = (1 + 2σ_kms / 3e5) * minimum(λ)
    i = findall(λ .< ee)
    out[i] = Spline1D(grid_λ, grid_y, k=1, bc="extrapolate")(λ[i])

    ee = (1 - 2σ_kms / 3e5) * maximum(λ)
    i = findall(λ .> ee)
    out[i] = Spline1D(grid_λ, grid_y, k=1, bc="extrapolate")(λ[i])

    # @gp    λ data "w lp t 'Data'"
    # @gp :- λ out  "w lp t 'Convolution'"
    return out
end


function delta_conv_gauss(λ, data, σ_kms; oversampling::Int=1)
    step = log10(σ_kms / 3e5 + 1) / oversampling
    grid_λ = 10. .^ range(log10(minimum(λ) * (1 - 5 * σ_kms/3e5)),
                          log10(maximum(λ) * (1 + 5 * σ_kms/3e5)),
                          step=step)
    out = grid_λ .* 0.

    for ii in 1:length(data)
        out .+= data[ii] .* gauss(grid_λ, λ[ii], σ_kms / 3e5 * λ[ii])
    end

    return grid_λ, out
end



#=

x, y = QSFit.delta_conv_gauss([1e3, 3e3, 1e4, 3e4], [0., 1., 1., 0.], 3000.)
QSFit.int_tabulated(x, y)  # should be ~ 2
@gp xlog=true x y "w lp"

y2 = QSFit.conv_gauss(x, y, 10000.)
@gp :- x y2 "w lp"
QSFit.int_tabulated(x, y)  # still ~ 2

=#


function convol(v, _k)
    k = reverse(_k)
    nk = length(k)
    @assert nk < length(v)
    @assert mod(nk, 2) == 1
    r = div(nk, 2)
    out = fill(0., length(v))
    for i in r+1:length(v)-r
        out[i-r:i+r] .+= v[i] .* k
    end
    return out ./ sum(abs.(k))
end
