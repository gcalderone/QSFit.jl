function voigt(x, σ, γ)
    z = (x + im * γ) / σ / sqrt(2)
    return real(faddeeva(z)) / (σ * sqrt(2pi))
end

# J.J.Olivero and R.L. Longbothum in Empirical fits to the Voigt line width: A brief review, JQSRT 17, P233, 1977
# https://ui.adsabs.harvard.edu/abs/1977JQSRT..17..233O/abstract
# http://snst-hu.lzu.edu.cn/zhangyi/ndata/Voigt_profile.html
function voigt_fwhm(σ, γ)
    @assert σ > 0
    @assert γ > 0
    fg = σ * 2.355   #  2 * sqrt(2 * log(2)) = 2.3548200450309493
    fl = γ * 2
    return fl * 0.5346 + sqrt(0.2166 * fl^2 + fg^2)
    #=
    Introducing dampening parameter:
    a = fl / fg
    fwhm = fl * (0.5346 + sqrt(0.2166 + (1/a)^2))
    =#
end

# The inverse function is
function voigt_σγ(fwhm, log_a)
    a = 10. ^log_a
    fl = fwhm / (0.5346 + sqrt(0.2166 + (1/a)^2))
    fg = fl / a
    return (fg / 2.355, fl / 2)
end
#=
σ, γ = 1.2, 3.4
log_a = log10((γ * 2) / (σ * 2.35482));
fwhm = QSFit.voigt_fwhm(σ, γ)
@info "" σ γ log_a fwhm
QSFit.voigt_σγ(fwhm, log_a), (σ, γ)
=#


# ____________________________________________________________________
# SpecLineVoigt
#
mutable struct SpecLineVoigt <: AbstractSpecLineComp
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    log_a::Parameter
    voff::Parameter

    function SpecLineVoigt(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  Parameter(0))
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.log_a.low  = -2
        out.log_a.high =  2
        out.center.fixed = true
        return out
    end
end

function evaluate!(ceval::CompEval{SpecLineVoigt, Domain{1}},
                   norm, center, fwhm, log_a, voff)
    x0 = center - (voff / 3.e5) * center
    σ, γ = voigt_σγ(fwhm / 3.e5 * center, log_a)
    X = coords(ceval.domain) .- x0
    # ceval.buffer .= 0.
    # i = findall(abs.(X) .< 2 * fwhm)
    ceval.buffer .= norm .* voigt.(X, σ, γ)
end
