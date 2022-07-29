# ____________________________________________________________________
# SpecLineAsymmGauss
#
mutable struct SpecLineAsymmGauss <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    asymm::Parameter
    index::Vector{Int}  # optimization

    function SpecLineAsymmGauss(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  Parameter(0),
                  Vector{Float64}(), true)

        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

function prepare!(comp::SpecLineAsymmGauss, domain::Domain{1})
    comp.index = collect(1:length(domain))
    return fill(NaN, length(domain))
end

function evaluate!(buffer::Vector{Float64}, comp::SpecLineAsymmGauss, x::Domain{1},
                   norm, center, fwhm, voff, asymm)
    buffer[comp.index] .= 0.
    empty!(comp.index)

    x0 = center - (voff / 3.e5) * center
    hwhm = fwhm / 3.e5 * center / 2  # Note: this is in `center` units

    sigma0 = hwhm / (2.355 / 2)
    sigma = 2. * sigma0 ./ (1 .+ exp.(asymm .* (x .- x0) ./ 2 ./ sigma0))
    X = (x[:] .- x0) ./ sigma
    i = findall(abs.(X) .< 4) # optimization
    append!(comp.index, i)
    buffer[i] .= norm * exp.(-X[i].^2 ./ 2) ./ sqrt(2pi) * sigma0
end
