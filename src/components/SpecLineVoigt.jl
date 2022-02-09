# ____________________________________________________________________
# SpecLineVoigt
#
mutable struct SpecLineVoigt <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    log_a::Parameter
    voff::Parameter
    spec_res_kms::Float64

    function SpecLineVoigt(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  Parameter(0),
                  0.)
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.log_a.low  = -2
        out.log_a.high =  2
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

function prepare!(comp::SpecLineVoigt, domain::Domain{1})
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::SpecLineVoigt, x::Domain{1},
                   norm, center, fwhm, log_a, voff)
    x0 = center - (voff / 3.e5) * center
    σ_res = comp.spec_res_kms / 3.e5 * center
    σ, γ = voigt_σγ(fwhm      / 3.e5 * center, log_a)
    σ = sqrt(σ^2 + σ_res^2)
    X = x .- x0

    function profile(x)
        return norm * voigt.(x, σ, γ)
    end
    map!(profile, buffer, X)
end
