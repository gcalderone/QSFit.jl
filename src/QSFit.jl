module QSFit

export Source, Spectrum, add_spec!, fit!, plot

import DataFitting: AbstractDomain, Domain_1D, Domain_2D,
    Parameter, AbstractComponent, AbstractComponentData,
    cdata, evaluate!, fit!

using CMPFit, DataFitting, Gnuplot, ReusePatterns, StructC14N
using Statistics, DataFrames, DelimitedFiles, Interpolations, Printf
using Unitful, UnitfulAstro, Parameters
#using FFTW

DataFitting.@enable_CMPFit
DataFitting.showsettings.fixedpars = false

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

@quasiabstract struct Source
    name::String
    z::Float64
    ebv::Float64
    flux2lum::Float64
    lines::Vector{SpectralLine}
    log::IO
    domain::Vector{DataFitting.Domain_1D}
    data::Vector{DataFitting.Measures_1D}

    function Source(name, z, ebv; log="", cosmo=qsfit_cosmology())
        @assert z > 0
        @assert ebv > 0
        ld = luminosity_dist(cosmo, float(z)) # * UnitfulAstro.Gpc
        ld = uconvert(u"cm", ld)
        flux2lum = 4pi * ld^2 * unit_flux() / unit_lum()
        @assert typeof(flux2lum) == Float64
        stream = stdout
        (log != "")  &&  (stream = open(log, "w"))
        return new(string(name), float(z), float(ebv), flux2lum,
                   known_lines(), stream,
                   Vector{DataFitting.Domain_1D}(), Vector{DataFitting.Measures_1D}())
    end
end


add_spec!(qsfit::Source, data::Spectrum) = error("No recipe for a generic `Source` object")
fit!(qsfit::Source) = error("No recipe for a generic `Source` object")

include("recipes/TypeI/general/module.jl")
include("recipes/TypeII/MCG-03-58-007/module.jl")

end  # module
