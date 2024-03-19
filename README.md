# QSFit

Quasar Spectral FITting package - http://qsfit.inaf.it/

** Warning **: the software is still under development...


## Install
```julia
using Pkg
Pkg.add("GModelFit.jl")
Pkg.add(url="https://github.com/lnicastro/GModelFitViewer.jl", rev="master")
Pkg.add(url="https://github.com/gcalderone/QSFit.jl", rev="master")
```

If you want to update from a previous version:
```julia
using Pkg
Pkg.update("GModelFit")
Pkg.update("GModelFitViewer")
Pkg.update("QSFit")
```



## Example
```julia
using QSFit

using QSFit.LineFitRecipes
spec = Spectrum(Val(:SDSS_DR10), "spec-0752-52251-0323.fits", label="My SDSS source")
recipe = Recipe(InteractiveLineFit, redshift=0.3806)
res = analyze(recipe, spec)
display(res.bestfit)

using Gnuplot
@gp res

using GModelFitViewer
viewer(res)


using QSFit.QSORecipes
spec = Spectrum(Val(:SDSS_DR10), "spec-0752-52251-0323.fits", label="My SDSS source")
recipe = Recipe(Type1, redshift=0.3806)
res = analyze(recipe, spec)
display(res.bestfit)
@gp res; viewer(res)

using GModelFit
abstract type MyRecipe <: Type1 end
import QSFit.QSORecipes.add_qso_continuum!
function add_qso_continuum!(recipe::Recipe{T}, resid::GModelFit.Residuals) where T <: MyRecipe
	@invoke add_qso_continuum!(recipe::Recipe{<: supertype(T)}, resid)
    resid.meval.model[:QSOcont].alpha.fixed = true
end
spec = Spectrum(Val(:SDSS_DR10), "spec-0752-52251-0323.fits", label="My SDSS source")
recipe = Recipe(MyRecipe, redshift=0.3806)
res = analyze(recipe, spec)
display(res.bestfit)


using GModelFit, Statistics
function add_qso_continuum!(recipe::Recipe{T}, resid::GModelFit.Residuals) where T <: MyRecipe
    λ = coords(domain(resid.data))
    comp = QSFit.sbpl(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(resid.data))
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    resid.meval.model[:QSOcont] = comp
    push!(resid.meval.model[:Continuum].list, :QSOcont)
    GModelFit.update!(resid.meval)
end
spec = Spectrum(Val(:SDSS_DR10), "spec-0752-52251-0323.fits", label="My SDSS source")
recipe = Recipe(MyRecipe, redshift=0.3806)
res = analyze(recipe, spec)
display(res.bestfit)
```
