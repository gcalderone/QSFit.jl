# ____________________________________________________________________
# SpecLineGauss
#
mutable struct SpecLineGauss <: AbstractSpecLineComp
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter

    function SpecLineGauss(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0))
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.center.fixed = true
        return out
    end
end

function evaluate!(ceval::CompEval{SpecLineGauss, Domain{1}},
                   norm, center, fwhm, voff)
    x0 = center - (voff / 3.e5) * center
    σ     = fwhm            / 2.355 / 3.e5 * center
    X = (coords(ceval.domain) .- x0) ./ σ

    function profile(x)
        (abs(x) > 4)  &&  (return 0.)
        return norm * exp(-x^2 / 2) / sqrt(2pi) / σ
    end
    map!(profile, ceval.buffer, X)
end
