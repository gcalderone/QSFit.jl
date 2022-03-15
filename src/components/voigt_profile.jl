function voigt(x, σ, γ)
    faddeeva(z) = erfcx(-im * z)
    z = (x + im * γ) / σ / sqrt(2)
    return real(faddeeva(z)) / (σ * sqrt(2pi))
end

# J.J.Olivero and R.L. Longbothum in Empirical fits to the Voigt line width: A brief review, JQSRT 17, P233, 1977
# http://snst-hu.lzu.edu.cn/zhangyi/ndata/Voigt_profile.html
function voigt_fwhm(σ, γ)
    @assert σ > 0
    @assert γ > 0
    fg = σ * 2.355
    fl = γ * 2
    return fl * 0.5346 + sqrt(0.2166 * fl^2 + fg^2)
    #=
    Introducing dampening parameter:
    a = fl / fg
    fwhm = fl * (0.5346 + sqrt(0.2166 + (1/a)^2))
    =#
end

function voigt_σγ(fwhm, log_a)
    a = 10. ^log_a
    fl = fwhm / (0.5346 + sqrt(0.2166 + (1/a)^2))
    fg = fl / a
    return (fg / 2.355, fl / 2)
end
# log_a = log10((γ * 2) / (σ * 2.355));  QSFit.voigt_σγ(QSFit.voigt_fwhm(σ, γ), log_a), (σ, γ)


function voigt_γ(fwhm, σ)
    A = 0.5346
    B = 0.2166
    C = (σ * 2.355)^2
    # fwhm = A * fl + sqrt(B * fl^2 + C)
    fl = A * fwhm / (A^2 - B) - sqrt((A^2 * C - B * C + B * fwhm^2) / ((A^2 - B)^2))
    γ = fl / 2
    @assert (fwhm - voigt_fwhm(σ, γ)) / fwhm < 1e-5
    return γ
end
# QSFit.voigt_γ(QSFit.voigt_fwhm(σ, γ), σ), γ
