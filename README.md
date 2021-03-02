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
@gp :- x[:] ceval.buffer "w l tit 'Gaussian' lw 2"
@info "FWHM=" QSFit.estimate_fwhm(x[:], ceval.buffer) / center * 3e5
@info "Area=" QSFit.int_tabulated(x[:], ceval.buffer)

comp = QSFit.SpecLineAsymmGauss(center);
comp.fwhm.val = fwhm;
comp.asymm.val = -1;
ceval = GFit.CompEval(comp, x);
GFit.evaluate_cached(ceval);
@gp :- x[:] ceval.buffer "w l tit 'Asymm. Gaussian' lw 2"
@info "FWHM=" QSFit.estimate_fwhm(x[:], ceval.buffer) / center * 3e5
@info "Area=" QSFit.int_tabulated(x[:], ceval.buffer)

comp = QSFit.SpecLineLorentz(center);
comp.fwhm.val = fwhm;
ceval = GFit.CompEval(comp, x);
GFit.evaluate_cached(ceval);
@gp :- x[:] ceval.buffer "w l tit 'Lorentzian' lw 2"
@info "FWHM=" QSFit.estimate_fwhm(x[:], ceval.buffer) / center * 3e5
@info "Area=" QSFit.int_tabulated(x[:], ceval.buffer)

```
