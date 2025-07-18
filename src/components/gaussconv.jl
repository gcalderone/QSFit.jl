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


import GModelFit: unfolded_domain, prepare!, apply_ir!

mutable struct GaussConv <: GModelFit.AbstractInstrumentResponse
    oversampling::Int
    Rspec::Float64 # λ / Δλ
    kernel::Vector{Float64}
    buffer::Vector{Float64}
    M::SparseMatrixCSC{Float64, Int64}
    function GaussConv(Rspec; oversampling::Int=2)
        @assert oversampling >= 2 "Oversampling must be at least 2"
        grid = -5:(1. / oversampling):5
        kernel = gauss(grid, 0., 1) ./ oversampling # divided by "oversampling" to ensure sum(kernel) = 1
        return new(oversampling, Rspec, kernel, Vector{Float64}(), sparse([0.]))
    end
end

function prepare!(IR::GaussConv, folded_domain::GModelFit.AbstractDomain)
    d = unfolded_domain(IR, folded_domain)
    empty!(IR.buffer)
    append!(IR.buffer, fill(0., length(d)))

    # Pre-compute interpolation matrix
    x = coords(d)
    X = coords(folded_domain)
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
end


# Sampling resolution is twice the spectral resolution
unfolded_domain(IR::GaussConv, folded_domain::GModelFit.AbstractDomain) =
    Domain(logregular_grid(coords(folded_domain), IR.oversampling * IR.Rspec))

function apply_ir!(IR::GaussConv,
                   folded_domain::GModelFit.AbstractDomain, folded::Vector,
                   unfolded_domain::GModelFit.AbstractDomain, unfolded::Vector)
    direct_conv1d!(IR.buffer, unfolded, IR.kernel, Val(:edge_mirror))
    folded .= IR.M * IR.buffer
end
