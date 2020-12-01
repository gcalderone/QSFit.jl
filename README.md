# QSFit

Quasar Spectral FITting package - http://qsfit.inaf.it/

** WORK IN PROGRESS **


```julia
center = 1000.
fwhm = 5e3

x = Domain(range(center - 3 * fwhm/3e5*center,
                 center + 3 * fwhm/3e5*center,
                 length=1000))

@gp "set grid"

comp = QSFit.SpecLineGauss(center);
comp.fwhm.val = fwhm;
ceval = GFit.CompEval(comp, x);
GFit.evaluate_cached(ceval);
@gp :- x[1] ceval.buffer "w l tit 'Gaussian' lw 2"
@info "FWHM=" QSFit.estimate_fwhm(x[1], ceval.buffer) / center * 3e5 / fwhm
@info "Area=" QSFit.estimate_area(x[1], ceval.buffer)

comp = QSFit.SpecLineAsymmGauss(center);
comp.fwhm.val = fwhm;
comp.asymm.val = -1;
ceval = GFit.CompEval(comp, x);
GFit.evaluate_cached(ceval);
@gp :- x[1] ceval.buffer "w l tit 'Asymm. Gaussian' lw 2"
@info "FWHM=" QSFit.estimate_fwhm(x[1], ceval.buffer) / center * 3e5 / fwhm
@info "Area=" QSFit.estimate_area(x[1], ceval.buffer)

comp = QSFit.SpecLineLorentz(center);
comp.fwhm.val = fwhm;
ceval = GFit.CompEval(comp, x);
GFit.evaluate_cached(ceval);
@gp :- x[1] ceval.buffer "w l tit 'Lorentzian' lw 2"
@info "FWHM=" QSFit.estimate_fwhm(x[1], ceval.buffer) / center * 3e5 / fwhm
@info "Area=" QSFit.estimate_area(x[1], ceval.buffer)

```

```julia
using Revise, GFit, QSFit

cd("/home/gcalderone/my/work/2020/ChangingLook/CL_1ES_1927p654")
source = QSO{DefaultRecipe}("1ES 1927+654", 0.019422, ebv=0.077)
add_spec!(source, Spectrum(Val(:ASCII), "AT2018zf/AT2018zf_optspec_20180514.dat", columns=[1,2]))
for id in 1:length(source.domain)
    source.data[id].unc .= 0.05 .* source.data[id].val;
    source.data[id].val .*= 1e17;
    source.data[id].unc .*= 1e17;
end

(model, bestfit) = fit!(source)
```
