function voigt(x, σ, γ)
    faddeeva(z) = erfcx(-im * z)
    z = (x + im * γ) / σ / sqrt(2)
    return real(faddeeva(z)) / (σ * sqrt(2pi))
    return ret
end
#=
@gp "unset grid" xr=[-10,10]
x = -100:0.01:100
σ = 1.53;  γ = 0. ;  vp = voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp))
σ = 1.30;  γ = 0.5;  vp = voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp))
σ = 0.01;  γ = 1.8;  vp = voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp))
σ = 1.  ;  γ = 1. ;  vp = voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp))
=#

# ____________________________________________________________________
# SpecLineVoigt
#

mutable struct SpecLineVoigt <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm_g::Parameter
    fwhm_l::Parameter
    voff::Parameter
    norm_integrated::Bool

    function SpecLineVoigt(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(3000),
                  Parameter(0),
                  true)

        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm_g.low = 0
        out.fwhm_l.low = 0
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

function prepare!(comp::SpecLineVoigt, domain::Domain{1})
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::SpecLineVoigt, x::Domain{1},
                   norm, center, fwhm_g, fwhm_l, voff)

    x0 = center - (voff / 3.e5) * center
    sigma = fwhm_g / 3.e5 * center / 2.355
    gamma = fwhm_l / 3.e5 * center / 2.
    buffer .= norm .* voigt.(x .- x0, sigma, gamma)
    comp.norm_integrated  ||  (buffer .*= (sigma * sqrt(2pi)))
end

#=
    x = Domain(500:1:1500.)
    comp = QSFit.SpecLineVoigt(1000.)
    comp.fwhm_g.val = 3e4
    comp.fwhm_l.val = 3e4
    ceval = GFit.CompEval(comp, x)
    evaluate!(ceval)
    @gp x[:] ceval.buffer ./ maximum(ceval.buffer) "w l"
=#

