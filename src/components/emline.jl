# ____________________________________________________________________
# emline
#
mutable struct emline <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    profile::Symbol
    function emline(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  :Gaussian)
        out.norm.low = 0
        out.center.low = 0        
        out.fwhm.low = 0        
        out.voff.low = 0        
        out.center.free = false
        return out
    end
end

ceval_data(domain::Domain_1D, comp::emline) = (nothing, length(domain))

function maxvalue(comp::emline)
    d = Domain([comp.center.val])
    ceval = CompEval(d, comp)
    evaluate(ceval, comp.norm.val, comp.center.val, comp.fwhm.val, comp.voff.val)
    return ceval.eval[1]
end

function evaluate(c::CompEval{Domain_1D, emline},
                  norm, center, fwhm, voff)
    x = c.domain[1]

    if c.comp.profile == :Lorentzian
        x0 = center - (voff / 3.e5) * center
        xx = (x .- x0) ./ (fwhm / 3.e5 * center / 2.)
        c.eval .= norm ./ (1 .+ xx.^2.)
        return
    end

    @assert c.comp.profile == :Gaussian
    x0 = center - (voff / 3.e5) * center
    sigma = (fwhm  / 3.e5) * center / 2.35
    #improve performance (SQRT(10*2.) = 4.4721360)
    i = findall(abs.(x .- x0) .< 4.4721360 * sigma)
    ee = ((x[i] .- x0) ./ sigma).^2. ./ 2.
    line = norm * exp.(-ee) ./ 2.50663 ./ sigma #SQRT(2*!PI) = 2.50663
    c.eval .= 0.
    c.eval[i] .= line
end

