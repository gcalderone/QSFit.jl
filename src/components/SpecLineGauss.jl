# ____________________________________________________________________
# SpecLineGauss
#
mutable struct SpecLineGauss <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    spec_res_kms::Float64
    norm_integrated::Bool

    function SpecLineGauss(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  0., true)

        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

function prepare!(comp::SpecLineGauss, domain::Domain{1})
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::SpecLineGauss, x::Domain{1},
                   norm, center, fwhm, voff)
    x0 = center - (voff / 3.e5) * center
    sigma_line = fwhm / 2.355      / 3.e5 * center
    sigma_spec = comp.spec_res_kms / 3.e5 * center
    sigma = sqrt(sigma_line^2 + sigma_spec^2)

    map!(x -> begin
         (abs.(x - x0) .> 4sigma)  &&  (return 0.)
         ret = norm * exp(-((x - x0) / sigma)^2 / 2)
         comp.norm_integrated  &&  (ret /= (sqrt(2pi) * sigma))
         return ret
         end, buffer, x.coords[:,1])
end


#=
    x = Domain(500:1:1500.)
    comp = QSFit.SpecLineGauss(1000.)
    comp.fwhm.val = 3e4
    ceval = GFit.CompEval(comp, x)
    GFit.evaluate_cached(ceval)
    @gp x[:] ceval.buffer ./ maximum(ceval.buffer) "w l"
=#
