# ____________________________________________________________________
# SpecLineGauss
#
mutable struct SpecLineGauss <: AbstractComponent
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
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

compeval_cdata(comp::SpecLineGauss, domain::Domain_1D) = collect(1:length(domain))
compeval_array(comp::SpecLineGauss, domain::Domain_1D) = fill(NaN, length(domain))

function maxvalue(comp::SpecLineGauss)
    ceval = CompEval(comp, Domain([comp.center.val]))
    evaluate_cached(ceval)
    return ceval.buffer[1]
end

function evaluate(c::CompEval{SpecLineGauss, Domain_1D},
                  norm, center, fwhm, voff)
    c.buffer[c.cdata] .= 0.
    empty!(c.cdata)

    x = c.domain[1]
    x0 = center - (voff / 3.e5) * center
    hwhm = fwhm / 3.e5 * center / 2  # Note: this is in `center` units

    sigma = hwhm / (2.355 / 2)
    X = (x .- x0) ./ sigma
    i = findall(abs.(X) .< 4)  # optimization
    append!(c.cdata, i)
    c.buffer[i] .= (norm / sqrt(2pi) / sigma) * exp.(-X[i].^2 ./ 2)
end


#=
    x = Domain(500:1:1500.)
    comp = QSFit.SpecLineGauss(1000.)
    comp.fwhm.val = 3e4
    ceval = GFit.CompEval(comp, x)
    GFit.evaluate_cached(ceval)
    @gp x[1] ceval.buffer ./ maximum(ceval.buffer) "w l"
=#
