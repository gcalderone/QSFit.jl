# ____________________________________________________________________
# SpecLineGauss
#
mutable struct SpecLineGauss <: AbstractSpecLineComp
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    span::Float64

    function SpecLineGauss(center::Number; span=5)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  span)
        @assert center > 0
        @assert span > 0
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
    σ  = fwhm   / 2.355 / 3.e5  * center
    x = coords(ceval.domain)

    for i in 1:length(x)
        X = x[i] - x0
        if abs(X) < ceval.comp.span * σ
            ceval.buffer[i] = norm * exp(-(X / σ)^2 / 2) / sqrt(2pi) / σ
        else
            ceval.buffer[i] = 0.
        end
    end
end
