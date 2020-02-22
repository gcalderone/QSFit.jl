#module QSFIT

import DataFitting: AbstractDomain, Domain_1D, Domain_2D,
    Parameter, AbstractComponent, AbstractComponentData,
    cdata, evaluate!

export QSFit, add_data!, read_sdss_dr10, run, plot

using CMPFit, DataFitting, Gnuplot, ReusePatterns
using Serialization, Statistics, DataFrames, DelimitedFiles, Interpolations, Printf
using Unitful, UnitfulAstro, Parameters
#using FFTW

DataFitting.@enable_CMPFit
DataFitting.showsettings.fixedpars = false
const showstep = true

include("utils.jl")
include("ccm_unred.jl")
include("components/emline.jl")
include("components/powerlaw.jl")
include("components/hostgalaxy.jl")
include("components/ironopt.jl")
include("components/ironuv.jl")
include("components/balmercont.jl")
#include("components/test_components.jl")
include("cosmology.jl")
include("Spectrum.jl")
include("spectral_lines.jl")

@quasiabstract mutable struct Options
    # The wavelength range used to for fitting.  Wavelengths outside
    # the range are ignored.
    λ_range::NTuple{2, Float64}
    
    Options(λ_range=(1210., 1e5)) = new(λ_range)
end


@quasiabstract struct QSFit
    name::String
    z::Float64
    ebv::Float64
    flux2lum::Float64
    lines::Vector{SpectralLine}
    log::IO
    options::concretetype(Options)
    domain::Vector{DataFitting.Domain_1D}
    data::Vector{DataFitting.Measures_1D}

    function QSFit(name, z, ebv; log="", cosmo=qsfit_cosmology())
        @assert z > 0
        @assert ebv > 0
        ld = luminosity_dist(cosmo, float(z)) # * UnitfulAstro.Gpc
        ld = uconvert(u"cm", ld)
        flux2lum = 4pi * ld^2 * unit_flux() / unit_lum()
        @assert typeof(flux2lum) == Float64
        stream = stdout
        (log != "")  &&  (stream = open(log, "w"))
        return new(string(name), float(z), float(ebv), flux2lum,
                   known_lines(), stream, Options(),
                   Vector{DataFitting.Domain_1D}(), Vector{DataFitting.Measures_1D}())
    end
end


add_data!(qsfit::QSFit, data::Spectrum) = error("No recipe for a generic `QSFit` object")

fit!(qsfit::QSFit) = error("No recipe for a generic `QSFit` object")

#end  # module
