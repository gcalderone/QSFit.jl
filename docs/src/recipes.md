# Recipes

**QSFit.jl** provides the following ready-to-use recipes:

### `InteractiveLineFit`
Interactive emission line fitting. The user is asked to provide the initial guess wavelengths of emission lines by clicking on a plot (provided by the [Gnuplot.jl](https://gcalderone.github.io/Gnuplot.jl/stable/index.html) package), as well as to decide the line profile to be used. Example:
```julia
using QSFit, QSFit.LineFitRecipes

filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = Recipe(InteractiveLineFit, redshift=0.3806)
res = analyze(recipe, spec)
```
After interactive selection the above code will print the relevant options to replicate the same analysis using the `LineFit` recipe (see below);
  
### `LineFit`
Non-interactive line fitting.  The user should provide the relevant options such as the redshift, the wavelength range, the guess wavelengths of the emission lines, etc.  The latter can be obtained using the above mentioned `InteractiveLineFit` recipe.  Example:
```julia
using QSFit, QSFit.LineFitRecipes

filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = Recipe(LineFit, redshift=0.3806)
recipe.wavelength_range = [4610.1157888384, 5130.00163811462]
add_line!(recipe, QSFit.ATL.UnidentifiedTransition(4866.45), NarrowLine,BroadLine)
add_line!(recipe, QSFit.ATL.UnidentifiedTransition(5008.98), ForbiddenLine)
res = analyze(recipe, spec)
```

### `Type1`
Automatic spectral analysis of Type1 AGN and QSO at redshifts <~ 2.1.  Example:
```
using QSFit, QSFit.QSORecipes
filename = download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")
spec = Spectrum(Val(:SDSS_DR10), filename)
recipe = Recipe(Type1, redshift=0.3806, Av=0.21)
res = analyze(recipe, spec)
```


## Customizing built-in recipes
To be written.

```julia
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



## Define new recipes
To be written.

```julia
using GModelFit, QSFit
import QSFit: init_recipe!, preprocess_spec!, analyze

abstract type MyRecipe <: AbstractRecipeSpec end

function init_recipe!(recipe::Recipe{T}) where T <: MyRecipe
    @invoke init_recipe!(recipe::Recipe{<: AbstractRecipeSpec})
    recipe.wavelength_range = [1215, 7.3e3]
end

function preprocess_spec!(recipe::Recipe{T}, spec::QSFit.Spectrum) where T <: MyRecipe
    @invoke preprocess_spec!(recipe::Recipe{<: AbstractRecipeSpec}, spec)
end

function analyze(recipe::Recipe{T}, spec::QSFit.Spectrum, resid::GModelFit.Residuals) where T <: MyRecipe
    resid.mzer.config.ftol = 1.e-6
    model = resid.meval.model
    return GModelFit.minimize!(resid)
end
```
