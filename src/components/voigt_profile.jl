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


#=
using SpecialFunctions, GFit, QSFit, Gnuplot
@gp "set grid" xr=[-10,10]
x = -100:0.01:100
σ = 1.53;  γ = 1.e-3;  vp = QSFit.voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp)); print(QSFit.voigt_fwhm(σ, γ) - QSFit.estimate_fwhm(x, vp))
σ = 1.30;  γ = 0.5;    vp = QSFit.voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp)); print(QSFit.voigt_fwhm(σ, γ) - QSFit.estimate_fwhm(x, vp))
σ = 0.01;  γ = 1.8;    vp = QSFit.voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp)); print(QSFit.voigt_fwhm(σ, γ) - QSFit.estimate_fwhm(x, vp))
σ = 1.  ;  γ = 1. ;    vp = QSFit.voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp)); print(QSFit.voigt_fwhm(σ, γ) - QSFit.estimate_fwhm(x, vp))
σ = 1.  ;  γ = 5. ;    vp = QSFit.voigt.(x, σ, γ);  @gp :- x vp "w l"; println(int_tabulated(x, vp)); print(QSFit.voigt_fwhm(σ, γ) - QSFit.estimate_fwhm(x, vp))


x = Domain(950:0.01:1050.)
center = 1e3
fwhm = 2e3

function plot_line(x, type, title; style="", log_a=nothing, spec_res=0.)
    comp = type(1e3)
    comp.fwhm.val = 2e3
    comp.resolution = spec_res
    isnothing(log_a)  ||  (comp.log_a.val = log_a)
    ceval = GFit.CompEval(comp, x)
    evaluate!(ceval);
    @gp :- x[:] ceval.buffer "w l t '$title' $style"
    println(int_tabulated(x[:], ceval.buffer))
    QSFit.estimate_fwhm(x[:], ceval.buffer) / comp.center.val * 3e5
end

@gp "set grid"
plot_line(x, QSFit.SpecLineGauss  , "Gaussian")
plot_line(x, QSFit.SpecLineLorentz, "Lorentzian")
plot_line(x, QSFit.SpecLineVoigt  , "Voigt (log_a=-1)", log_a=-1, style="dt 2 lw 2")
plot_line(x, QSFit.SpecLineVoigt  , "Voigt (log_a=0)" , log_a=0 , style="lw 2")
plot_line(x, QSFit.SpecLineVoigt  , "Voigt (log_a=1)" , log_a=1 , style="dt 2 lw 2")


plot_line(x, spec_res=1e3, QSFit.SpecLineGauss  , "Gaussian")
plot_line(x, spec_res=1e3, QSFit.SpecLineLorentz, "Lorentzian")
plot_line(x, spec_res=1e3, QSFit.SpecLineVoigt  , "Voigt (log_a=1)" , log_a=1 , style="dt 2 lw 2")
plot_line(x, spec_res=1e3, QSFit.SpecLineVoigt  , "Voigt (log_a=0)" , log_a=0 , style="lw 2")
plot_line(x, spec_res=1e3, QSFit.SpecLineVoigt  , "Voigt (log_a=-1)", log_a=-1, style="dt 2 lw 2")
=#
