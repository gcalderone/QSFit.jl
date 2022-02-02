# ____________________________________________________________________
# SpecLineLorentz
#
mutable struct SpecLineLorentz <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    norm_integrated::Bool

    function SpecLineLorentz(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  true)
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

function prepare!(comp::SpecLineLorentz, domain::Domain{1})
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::SpecLineLorentz, x::Domain{1},
                   norm, center, fwhm, voff)
    x0 = center - (voff / 3.e5) * center
    hwhm = fwhm / 3.e5 * center / 2  # Note: this is in the same units as `center`
    map!(x -> begin
         (abs.(x) .> 20)  &&  (return 0.)
         ret = norm ./ (1 .+ x.^2.)
         comp.norm_integrated  &&  (ret /= (pi * hwhm))
         return ret
         end,
         buffer, (x .- x0) ./ hwhm)
end


#=
    x = Domain(500:1:1500.)
    comp = QSFit.SpecLineLorentz(1000.)
    comp.fwhm.val = 3e4
    ceval = GFit.CompEval(comp, x)
    evaluate!(ceval)
    @gp x[:] ceval.buffer ./ maximum(ceval.buffer) "w l"
=#
