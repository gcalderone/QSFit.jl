# QSFit

Quasar Spectral FITting package - http://qsfit.inaf.it/

** Warning **: the software is still under development...


## Install
```julia
using Pkg
Pkg.add("GModelFit.jl")
Pkg.add(url="https://github.com/lnicastro/GModelFitViewer.jl")
Pkg.add(url="https://github.com/gcalderone/QSFit.jl", rev="master")
```

## Example
```julia
using QSFit, Gnuplot, GModelFitViewer

using QSFit.LineFitRecipes
source = QSFit.Source("My SDSS source", z=0.3806,
	Spectrum(Val(:SDSS_DR10), "/home/gcalderone/my/work/software/qsfit/data/spec-0752-52251-0323.fits"))
recipe = RRef(InteractiveLineFitRecipe)
res = analyze(recipe, source)
@gp res; viewer(res)

using QSFit.QSORecipes
source = QSFit.Source("My SDSS source", z=0.3806,
	Spectrum(Val(:SDSS_DR10), "/home/gcalderone/my/work/software/qsfit/data/spec-0752-52251-0323.fits"))
recipe = QSFit.RRef(Type1Recipe)
res = analyze(recipe, source)
@gp res; viewer(res)



abstract type MyRecipe <: DefaultRecipe end
import QSFit.add_qso_continuum!
function add_qso_continuum!(recipe::QSFit.RRef{T}, state::QSFit.State) where T <: MyRecipe
	@invoke add_qso_continuum!(recipe::QSFit.RRef{<: supertype(T)}, state)
    state.model[:QSOcont].alpha.fixed = true
end


using GModelFit, Statistics
function add_qso_continuum!(recipe::QSFit.RRef{T}, state::QSFit.State) where T <: MyRecipe
    λ = coords(domain(state.model))
    comp = QSFit.sbpl(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(state.pspec.data))
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    state.model[:QSOcont] = comp
    push!(state.model[:Continuum].list, :QSOcont)
    GModelFit.update!(state.model)
end


myrecipe = QSFit.RRef(MyRecipe)
res = QSFit.analyze(myrecipe, source)
viewer(res)



```
