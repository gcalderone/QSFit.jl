using SparseArrays

#=
Notes on spectral resolution
A regular log-λ grid is characterized by a constant step:
logstep = log10(λ_i+1) - log10(λ_i)  =  log10(λ_i+1 / λ_i)

hence:
λ_i+1 / λ_i  =  const.

by subtracting a constant 1 to both sides
λ_i+1 / λ_i - 1  =  (λ_i+1 - λ_i) / λ_i  =  Δλ / λ  = const.

We interpret the constant as the reciprocal of the resolution:
Δλ / λ  =  1/R

i.e.
λ_i+1 / λ_i = 1/R + 1
logstep = log10(1/R + 1)

The resolution can also be expressed in km/s:
R = c / σ_kms
=#

# Calculate logstep corresponding to resolution R = λ / Δλ
logstep(R) = log10(1/R + 1)

function logregular_grid(x, R)
    ee = extrema(x)
    ee = (ee[1] - ee[1]/R, ee[2] + ee[2]/R)
    return 10. .^collect(log10(ee[1]):logstep(R):log10(ee[2]))
end


import GModelFit: model_domain, apply_ir!

mutable struct GaussConv <: GModelFit.AbstractInstrumentResponse
    oversampling::Int
    Rspec::Float64 # λ / Δλ
    kernel::Vector{Float64}
    buffer::Vector{Float64}
    M::SparseMatrixCSC{Float64, Int64}
    function GaussConv(Rspec; oversampling::Int=2)
        nsigma = 5
        grid = -ceil(2 * nsigma):ceil(2 * nsigma)
        kernel = gauss(grid, 0., 2.)
        return new(oversampling, Rspec, kernel, Vector{Float64}())
    end
end

# Sampling resolution is twice the spectral resolution
function model_domain(IR::GaussConv, data_domain::GModelFit.AbstractDomain)
    d = Domain(logregular_grid(coords(data_domain), IR.oversampling * IR.Rspec))
    append!(IR.buffer, fill(0., length(d)))

    x = coords(d)
    X = coords(data_domain)
    @assert issorted(x)
    M = fill(0., length(x), length(X))
    for j in 1:length(X)
        i1 = findlast( x .< X[j])
        i2 = findfirst(x .> X[j])
        f = (X[j] - x[i1]) / (x[i2] - x[i1])
        M[i1, j] = 1. - f
        M[i2, j] = f
    end
    IR.M = sparse(collect(M'))
    return d
end

function apply_ir!(IR::GaussConv,
                   data_domain::GModelFit.AbstractDomain, folded::Vector,
                   model_domain::GModelFit.AbstractDomain, unfolded::Vector)
    direct_conv1d!(IR.buffer, unfolded, IR.kernel, Val(:edge_mirror))
    folded .= IR.M * IR.buffer
end

#=
function mydelta(x)
    out = x .* 0.
    mm = div(length(out), 6)
    out[argmin(abs.(x .- 1250))] = 1.;
    out[argmin(abs.(x .- 1750))] = 1.;
    out[argmin(abs.(x .- 2600))] = 1.;
    return out
end

x = 10 .^(3:0.0001:3.5);

model = Model(:orig => @fd(x -> mydelta(x)), :convol => QSFit.GaussConv(:orig, 2000))
meval = GModelFit.ModelEval(model, Domain(x))
GModelFit.scan_model!(meval)
c = GModelFit.update_eval!(meval)

ox = coords(meval.cevals[:orig].domain)
oy = meval.cevals[:orig].tpar.buffer
QSFit.int_tabulated(ox, oy), QSFit.int_tabulated(x, c)  # <-- these should be equal
@gp xlog=true xr=extrema(x) ox oy "w lp" x c "w lp" ox QSFit.gauss_broadening(ox, oy, 150.) "w lp"

i = argmin(abs.(x .- 1250))
QSFit.estimate_fwhm(x[(i-20):(i+20)], c[(i-20):(i+20)]) / x[i] * 3e5 / 2.355  # should be ~150 km/s, i.e. 3e5 / 2000
=#
