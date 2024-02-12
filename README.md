# QSFit

Quasar Spectral FITting package - http://qsfit.inaf.it/

** Warning **: the software is still under development...


## Install
```julia
using Pkg
Pkg.add(url="https://github.com/gcalderone/GModelFit.jl")
Pkg.add(url="https://github.com/lnicastro/GModelFitViewer.jl")
Pkg.add(url="https://github.com/gcalderone/QSFit.jl", rev="master")
```

## Example
```julia
using QSFit, GModelFitViewer

source = QSFit.Source("My SDSS source", 0.3806, ebv=0.)
add_spec!(source, Spectrum(Val(:SDSS_DR10), "/home/gcalderone/my/work/software/qsfit/data/spec-0752-52251-0323.fits"))
recipe = QSFit.RRef(DefaultRecipe)
res = QSFit.analyze(recipe, source)
viewer(res.bestfit, res.fitstats, res.pspec.data)
```
