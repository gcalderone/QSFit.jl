# Recipes

## Introduction

The essential steps to analyze a QSO spectrum typically involve the identification of relevant features (such as emission lines, continuum emission, host galaxy emission, etc.) and their association with a corresponding *component* in a model[^1] to fit the empirical data.  The resulting best fit parameter values (and associated uncertainties) will allow further analysis.

[^1]: Note: the terms *component* and *model* used here have the same meaning as in the [**GmodelFit.jl**](https://gcalderone.github.io/GModelFit.jl/concepts/) documentation.

While the above steps are always the same, the finer details of the analysis depend on several aspects such as the spectrograph being used, its spectral resolution, the wavelength coverage, the target redshift, the peculiar features in the spectrum, the accuracy required in the outcomes, etc.

To address such complexity **QSFit.jl** provides a few built-in *analysis recipes*, namely a sequence of steps aimed to build a spectral [model](https://gcalderone.github.io/GModelFit.jl/concepts/), and to fit it against an empirical 1D spectrum.  Each recipe comes with a number of pre-defined options which can easily be customized to address specific cases.  E.g., to check default options for the `Type1` recipe:
```@example abc
using QSFit, QSFit.QSORecipes
recipe = CRecipe{Type1}()
```
You can modify any option with (e.g.):
```julia
recipe.use_host_template = false
```

In case the customization via options is not enought, it is possible to [Define new recipes](@ref) by extending the built-in ones.
 

## Built-in recipes

The built-in recipes are defined in the `QSFit.LineFitRecipes` and `QSFit.QSORecipes` submodules.

### `InteractiveLineFit`
Interactive emission line fitting. The user is asked to provide the initial guess wavelengths of emission lines by clicking on a plot (provided by the [Gnuplot.jl](https://gcalderone.github.io/Gnuplot.jl/stable/index.html) package), as well as to decide the line profile to be used. Example:
```julia
using QSFit, QSFit.LineFitRecipes

filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = CRecipe{InteractiveLineFit}(redshift=0.3806)
res = analyze(recipe, spec)
```
After interactive selection the above code will print the relevant options to replicate the same analysis using the `LineFit` recipe (see below);

### `LineFit`
Non-interactive line fitting.  The user should provide the relevant options such as the redshift, the wavelength range, the guess wavelengths of the emission lines, etc.  The latter can be obtained using the above mentioned `InteractiveLineFit` recipe.

Example:
```julia
using QSFit, QSFit.LineFitRecipes

filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = CRecipe{LineFit}(redshift=0.3806)
recipe.wavelength_range = [4530.90809628009, 5392.50547045952]
recipe.lines = QSFit.SpecLineSet()
QSFit.add_line!(recipe, recipe.lines, 4864.77, NarrowLine, BroadLine)
QSFit.add_line!(recipe, recipe.lines, 5010.88, ForbiddenLine)
res = analyze(recipe, spec)
```

### `Type1`
Automatic spectral analysis of Type1 AGN and QSO at redshifts <~ 2.1. 

Example:
```julia
using QSFit, QSFit.QSORecipes
filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = CRecipe{Type1}(redshift=0.3806, Av=0.21)
res = analyze(recipe, spec)
```

 








## Define new recipes

As anticipated, *recipes* are essentially a sequence of steps to carry out the spectral analysis.  They are implemented adopting a pattern which allows to easily define new recipes by inheriting the functionalities of the existiung ones, and overriding only the specific steps which needs to be customized or improved.

From the implementation point of view a *recipe* is an abstract type inheriting from `AbstractRecipe`, and the steps are functions accepting that type as an argument.  Consider the following example where we define a recipe `A`, a few steps and substeps, as well as a main `analyze` function triggering execution of steps in the proper order:

```@example abc
abstract type A <: AbstractRecipe end

step1(t::Type{T}) where T <: A = println("Invoked step1 from recipe A")

substep2a(::Type{T}) where T <: A = println("  Invoked substep2a from recipe A")
substep2b(::Type{T}) where T <: A = println("  Invoked substep2b from recipe A")

function step2(::Type{T}) where T <: A
    println("Invoked step2 from recipe A")
    substep2a(T)
    substep2b(T)
end

function analyze(t::Type{T}) where T <: A
    println("Invoked analyze from recipe A")
    step1(t)
    step2(t)
end

analyze(A)
```

Suppose we want to reimplement `step1`, and to reuse `substep2b` by performing some extra operation after it has been executed.  We'll need to define a new recipe `B` inheriting from `A`, and implement new methods for the `step1` and `substep2b` functions:
```@example abc
abstract type B <: A end

step1(t::Type{T}) where T <: B = println("Invoked step1 from recipe B")

function substep2b(t::Type{T}) where T <: B
    @invoke substep2b(t::Type{<: supertype(T)})
    println("  Additional operation performed in substep2b from recipe B")
end

analyze(B)
```

The above pattern allows to maximize code reuse (i.e., avoid reinventing the wheel), as well as avoiding conditional branches (`if-then-else`) since the identification of the steps to be executed is performed via multiple dispatch.

The actual recipes code follows the above pattern except for two details:
- we use the `CRecipe{}` type in place of `::Type{}` since the latter is just an abstract type and can not contain any data, while the fomer is a concrete one and can store the recipe options;
- by convention, each recipe step should invoke the `@track_recipe` macro as first statement.  This way it is possible to track down which method is actually being executed (this requires invoking `QSFit.track_recipe(true)`);



In the following example we will define a new recipe named `MyRecipe` inheriting from `Type1` and
- replace the power law continuum component with a smoothly broken power law;
- smooth the host galaxy template with a 50 points boxcar average;

```julia
using GModelFit, QSFit, QSFit.QSORecipes, Statistics
import QSFit.QSORecipes: getdomain, getmodel, getdata

abstract type MyRecipe <: Type1 end

import QSFit.QSORecipes.add_qso_continuum!
function add_qso_continuum!(recipe::CRecipe{T}, fp::GModelFit.FitProblem, ith::Int) where T <: MyRecipe
    @track_recipe
    λ = coords(getdomain(fp, ith))

    comp = QSFit.sbpl(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(getdata(fp, ith)))
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.delta.val = 0.001
    comp.delta.fixed = true
	getmodel(fp, ith)[:QSOcont] = comp
    push!(getmodel(fp, ith)[:Continuum].list, :QSOcont)
end

import QSFit.QSORecipes.add_host_galaxy!
function add_host_galaxy!(recipe::CRecipe{T}, fp::GModelFit.FitProblem, ith::Int) where T <: MyRecipe
    @track_recipe
    @invoke add_host_galaxy!(recipe::CRecipe{<: Type1}, fp, ith)
	model = getmodel(fp, ith)
    t = model[:Galaxy].template
    n = 50
	model[:Galaxy].template[1+n:length(t)-n] .= [mean(t[i-n:i+n]) for i in 1+n:length(t)-n]
end

QSFit.track_recipe(true)
filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = CRecipe{MyRecipe}(redshift=0.3806, Av=0.21)
res = analyze(recipe, spec)
display(res.bestfit)
```
