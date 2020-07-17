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
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

ceval_data(domain::Domain_1D, comp::emline) = (collect(1:length(domain)), length(domain))

function maxvalue(comp::emline)
    d = Domain([comp.center.val])
    ceval = CompEval(d, comp)
    evaluate(ceval, comp.norm.val, comp.center.val, comp.fwhm.val, comp.voff.val)
    return ceval.eval[1]
end

function evaluate(c::CompEval{Domain_1D, emline},
                  norm, center, fwhm, voff)
    c.eval[c.cdata] .= 0.
    empty!(c.cdata)

    x = c.domain[1]
    x0 = center - (voff / 3.e5) * center
    hwhm = fwhm / 3.e5 * center / 2  # Note: this is in `center` units

    if c.comp.profile == :Lorentzian
        X = (x .- x0) ./ hwhm
        i = findall(abs.(X) .< 20)
        append!(c.cdata, i)
        c.eval[i] .= (norm / pi / hwhm) ./ (1 .+ X[i].^2.)
    else
        @assert c.comp.profile == :Gaussian
        sigma = hwhm / (2.355 / 2)
        X = (x .- x0) ./ sigma
        i = findall(abs.(X) .< 4)
        append!(c.cdata, i)
        c.eval[i] .= (norm / sqrt(2pi) / sigma) * exp.(-X[i].^2 ./ 2)
    end
end


#=
    x = Domain(500:1:1500.)
    comp = QSFit.emline(1000.)
    comp.fwhm.val = 3e4
    ceval = GFit.CompEval(x, comp)
    evaluate(ceval)
    @gp x[1] ceval.eval ./ maximum(ceval.eval) "w l"
    comp.profile = :Lorentzian
    evaluate(ceval)
    @gp :- x[1] ceval.eval ./ maximum(ceval.eval) "w l"
=#
