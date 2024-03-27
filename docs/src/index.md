```@setup abc
include("setup.jl")
```

# QSFit.jl

**QSO spectral fitting made easy !!**

[![Stars](https://img.shields.io/github/stars/gcalderone/QSFit.jl?style=social)](https://github.com/gcalderone/QSFit.jl)


!!! warn
    This software is under active development and details may change at any time without notice.  Also, documentation is not yet exhaustive.

QSFit started as an attempt to perform automatic analysis of optical spectra of AGNs and QSOs in a simple, replicable and shareable way. The first implementation in IDL language (see repo [here](https://github.com/gcalderone/qsfit) has been used to analyze 71261 AGN and QSO spectra from SDSS DR10.  The resulting spectral properties (such as emission line luminosities and widths, continuum slopes, etc.) are collected in a catalog in FITS format, as well as being publicly available here: [https://qsfit.inaf.it/](https://qsfit.inaf.it/).

The details of the spectral analysis are presented in a paper: [Calderone et al. 2017](http://adsabs.harvard.edu/abs/2017MNRAS.472.4051C) (also available on [arXiv](https://arxiv.org/abs/1612.01580)).

The code has now been ported to [Julia](https://julialang.org/) and the new package is dubbed **QSFit.jl**.  It comes with the following advantages with respect to the IDL version:
- Julia is released with a [MIT license](https://en.wikipedia.org/wiki/MIT_License) and doesn't require a paid license to be executed;
- It provides better performances and allows to distribute workload on multiple CPUs or multiple host;
- It provides several reusable [Components](@ref) and a few [Recipes](@ref) to perform automatic spectral analysis as well as interactive emissione line fitting;
- It exploits the concept of *customizable recipes* (see [Recipes](@ref)) to customize the analysis for specific purposes;

**QSFit.jl** relies on the [GModelFit.jl](https://gcalderone.github.io/GModelFit.jl/) package to perform spectral fitting.  A basic knowledge of such package is required to get the most out of **QSFit.jl**.


## Installation

In the Julia REPL type:
```julia
using Pkg
Pkg.add("GModelFit.jl")
Pkg.add("GModelFitViewer.jl")
Pkg.add(url="https://github.com/gcalderone/QSFit.jl", rev="0.1.0")
```

To test the package type `Pkg.test("QSFit")`.


## Basic usage

The most important types used in **QSFit.jl** are:
- `Spectrum`: it represent an observed spectrum of a AGN or QSO.  It contains both the observed wavelengths and flux densities (with associated undertainties), as well as a map of "good" (i.e. reliable) spectral bins and an indication of spectral resolution.  The `Spectrum` object can operate on spectra taken from any instrument (not only SDSS);

- `Recipe`: it is a container for a specific *recipe* to be used to analyze a `Spectrum`, and for all recipe-specific options.  **QSFit.jl** provides a few ready-to-use recipes to analyze the spectrum, which can optionally be customized by the user.  New recipes can also be implemented.

The typical workflow for a spectral analysis is as follows:
```@example abc
using QSFit, QSFit.QSORecipes
# Download a spectrum
filename = "/home/gcalderone/my/work/software/qsfit/data/spec-0752-52251-0323.fits" 
# download("http://dr10.sdss3.org/sas/dr10/sdss/spectro/redux/26/spectra/0752/spec-0752-52251-0323.fits")

# Read spectrum into a Spectrum Object
spec = Spectrum(Val(:SDSS_DR10), filename, label="My SDSS source")

# Create a Recipe object based on the Type1 recipe to analyze the spectrum
recipe = Recipe(Type1, redshift=0.3806, Av=0.21)

# Analyze the spectrum with the above recipe
res = analyze(recipe, spec)

# Display best fit parameter values:
show(res.bestfit)
```

The best fit model can be displayed using either [Gnuplot.jl](https://gcalderone.github.io/Gnuplot.jl/stable/index.html) (to open a plot window):
```@example abc
using Gnuplot
@gp res yrange=[0, 0.25]
saveas("example1") # hide
```
![](assets/example1.png)


or [GModelFitViewer.jl](https://github.com/lnicastro/GModelFitViewer.jl) (to display the plot in a browser):
```@example abc
using GModelFitViewer
viewer(res)
println() # hide
```
