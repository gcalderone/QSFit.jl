# ____________________________________________________________________
# SpecLineLorentz
#
mutable struct SpecLineLorentz <: AbstractSpecLineComp
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    span::Float64

    function SpecLineLorentz(center::Number; span=10)
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

function evaluate!(comp::SpecLineLorentz, domain::Domain{1}, output::Vector,
                   norm, center, fwhm, voff)
    x0 = center - (voff / 3.e5) * center
    γ  = fwhm / 2       / 3.e5  * center
    x = coords(domain)
    for i in 1:length(x)
        X = x[i] - x0
        if abs(X) < comp.span * γ
            output[i] = norm * γ / (pi * (X^2. + γ^2.))
        else
            output[i] = 0.
        end
    end
end
