module QSFit

export QSO, parent_recipe, add_spec!, logio, close_logio, PreparedSpectrum

import GFit: Domain, CompEval,
    Parameter, AbstractComponent, prepare!, evaluate!, fit!

using CMPFit, GFit
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, Dierckx, Printf, DataStructures
using Unitful, UnitfulAstro
using Dates
using Gnuplot
using TextParse, DataFrames

include("cosmology.jl")
include("ccm_unred.jl")
include("components/powerlaw.jl")
include("components/cutoff_powerlaw.jl")
include("components/hostgalaxy.jl")
include("components/ironopt.jl")
include("components/ironuv.jl")
include("components/balmercont.jl")
include("components/SpecLineGauss.jl")
include("components/SpecLineAsymmGauss.jl")
include("components/SpecLineLorentz.jl")
include("utils.jl")
include("convolutions.jl")
include("Spectrum.jl")


abstract type AbstractRecipe end

function default_options(::Type{T}) where T <: AbstractRecipe
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict{Symbol, Float64}(:default => 0.6)
    out[:skip_lines] = Vector{Symbol}()
    return out
end


struct QSO{T <: AbstractRecipe}
    name::String
    z::Float64
    mw_ebv::Float64
    cosmo::Cosmology.AbstractCosmology
    flux2lum::Float64
    logfile::Union{Nothing, String}
    options::OrderedDict{Symbol, Any}
    specs::Vector{Spectrum}
end

function QSO{T}(name, z; ebv=0., logfile=nothing, cosmo=default_cosmology()) where T <: AbstractRecipe
    @assert z > 0
    @assert ebv >= 0
    ld = uconvert(u"cm", luminosity_dist(cosmo, float(z)))
    flux2lum = 4pi * ld^2 * (scale_flux() * unit_flux()) / (scale_lum() * unit_lum())
    return QSO{T}(string(name), float(z), float(ebv), cosmo, flux2lum, logfile,
                  default_options(T), Vector{Spectrum}())
end

parent_recipe(source::QSO{T}) where T <: AbstractRecipe =
    QSO{supertype(T)}(getfield.(Ref(source), fieldnames(typeof(source)))...)

add_spec!(source::QSO, spec::Spectrum) =
    push!(source.specs, spec)

const logio_streams = Dict{String, IOStream}()

function logio(source::QSO)
    if isnothing(source.logfile)
        GFit.showsettings.plain = false
        return stdout
    end
    if !haskey(logio_streams, source.logfile)
        if isfile(source.logfile)
            f = open(source.logfile, "a")
        else
            f = open(source.logfile, "w")
        end
        logio_streams[source.logfile] = f
        println(f, "Timestamp: ", now())
        GFit.showsettings.plain = true
    end
    return logio_streams[source.logfile]
end

function close_logio(source::QSO)
    if haskey(logio_streams, source.logfile)
        close(logio_streams[source.logfile])
        delete!(logio_streams, source.logfile)
        GFit.showsettings.plain = false
    end
end

include("SpectralLines.jl")

struct PreparedSpectrum
    id::Int
    orig::Spectrum
    domain::GFit.Domain{1}
    data::GFit.Measures{1}
    lcs::OrderedDict{Symbol, LineComponent}
end


include("DefaultRecipe.jl")
#include("DefaultRecipe_multi.jl")
include("reduce.jl")
include("viewer.jl")
include("gnuplot.jl")
include("interactive_guess.jl")

end  # module
