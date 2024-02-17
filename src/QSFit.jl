module QSFit

export add_spec!, close_log

import GModelFit: Domain, CompEval,
    Parameter, AbstractComponent, prepare!, evaluate!

using CMPFit, GModelFit, SortMerge
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, Printf, DataStructures
using Unitful, UnitfulAstro
using Dates
using InteractiveUtils
using SpecialFunctions
using Gnuplot
using TextParse
using Cosmology

import Dierckx  # import (instead of using) avoid conflicts with evaluate()

import Base.show


include("ccm_unred.jl")
include("components/powerlaw.jl")
include("components/sbpl.jl")
include("components/cutoff_powerlaw.jl")
include("components/hostgalaxy.jl")
include("components/ironopt.jl")
include("components/ironuv.jl")
include("components/balmercont.jl")
include("components/voigt_profile.jl")

abstract type AbstractSpecLineComp <: AbstractComponent end
include("components/SpecLineGauss.jl")
include("components/SpecLineLorentz.jl")
include("components/SpecLineVoigt.jl")

include("utils.jl")
include("convolutions.jl")
include("Spectrum.jl")


qsfit_data() = artifact"qsfit_data"


# ====================================================================
struct Source
    name::String
    z::Union{Nothing, Float64}
    MW_ebv::Union{Nothing, Float64}
    specs::Vector{Spectrum}
    function Source(name; z=nothing, ebv=nothing)
        return new(string(name), z, ebv, Vector{Spectrum}())
    end
end

function show(io::IO, source::Source)
    println(io, "Source: ", source.name, ", z=", source.z, ", E(B-V) from Milky Way: ", source.MW_ebv)
    for i in 1:length(source.specs)
        print(io, "    ")
        show(io, source.specs[i])
    end
end

add_spec!(source::Source, spec::Spectrum) = push!(source.specs, spec)


# ====================================================================
struct PreparedSpectrum
    origlength::Int
    resolution::Float64
    flux2lum::Float64
    domain::GModelFit.Domain{1}
    data::GModelFit.Measures{1}
end

function PreparedSpectrum(source::Source, logio=stdout;
                          id=1,
                          dered::Union{Nothing, Function}=nothing,
                          cosmology::Union{Nothing, Cosmology.AbstractCosmology}=nothing)
    data = deepcopy(source.specs[id])
    println(logio, "Source: " * data.label)
    println(logio, "  spectrum ID: ", id)
    goodfraction = count(data.good) / length(data.good)
    println(logio, "  good fraction:: ", goodfraction)
    println(logio, "  resolution: ", @sprintf("%.4g", data.resolution), " km / s (FWHM)")

    # De-reddening
    if !isnothing(source.MW_ebv)  &&  !isnothing(dered)
        tmp = dered([1450, 3000, 5100.], source.MW_ebv)
        println(logio, "De-reddening factors @ 1450, 3000, 5100 AA: ", tmp)
        tmp = dered(data.λ, source.MW_ebv)
        data.flux .*= tmp
        data.err  .*= tmp
    end

    # Compute rest frame spectrum
    flux2lum = NaN
    if !isnothing(source.z)  &&  !isnothing(cosmology)
        ld = uconvert(u"cm", luminosity_dist(cosmology, source.z))
        flux2lum = 4pi * ld^2 * (scale_flux() * unit_flux()) / (scale_lum() * unit_lum())
        println(logio, "Using flux-to-lum. conversion factor: ", flux2lum)
        data.λ    ./= (1 + source.z)
        data.flux .*= flux2lum * (1 + source.z)
        data.err  .*= flux2lum * (1 + source.z)
    end

    ii = findall(data.good)
    domain = Domain(data.λ[ii])
    return PreparedSpectrum(length(data.flux), data.resolution, flux2lum,
                            domain, Measures(domain, data.flux[ii], data.err[ii]))
end


# ====================================================================
mutable struct State
    logfile::Union{Nothing, String}
    logio::Union{IOStream, Base.TTY}
    pspec::PreparedSpectrum
    model::GModelFit.Model
    user::OrderedDict{Symbol, Any}
end


function select_good!(state::State, good::Vector{Bool})
    i = findall(good)
    domain = Domain(coords(state.pspec.domain)[i])
    meas = Measures(domain,
                    values( state.pspec.data)[i],
                    uncerts(state.pspec.data)[i])
    state.pspec = PreparedSpectrum(state.pspec.origlength,
                                   state.pspec.resolution,
                                   state.pspec.flux2lum,
                                   domain, meas)
    state.model = Model(domain)
end


struct Results
    timestamp::DateTime
    elapsed::Float64
    source::Source
    pspec::Union{Nothing, PreparedSpectrum}
    bestfit::GModelFit.ModelSnapshot
    fitstats::GModelFit.FitStats
    reduced::OrderedDict{Symbol, Any}
end

function show(io::IO, res::Results)
    show(io, res.fitstats)
    println(io)
end


# ====================================================================
abstract type AbstractRecipe end

struct RRef{T <: AbstractRecipe}
    options::OrderedDict{Symbol, Any}

    function RRef(::Type{T}; kws...) where T <: AbstractRecipe
        out = new{T}(default_options(T))
        # Set options provided as keywords
        for (k, v) in kws
            out.options[Symbol(k)] = v
        end
        return out
    end
end

default_options(::Type{<: AbstractRecipe}) = OrderedDict{Symbol, Any}()
get_cosmology(::Type{<: AbstractRecipe}) = cosmology(h=0.70, OmegaM=0.3)   #S11
get_MW_deredd_function(::Type{<: AbstractRecipe}) = ccm_unred

reduce(recipe::RRef{<: AbstractRecipe}, state::State) = OrderedDict{Symbol, Any}()

function analyze(recipe::RRef{T}, source::Source; logfile=nothing, overwrite=false) where T <: AbstractRecipe
    timestamp = now()
    starttime = time()
    if isnothing(logfile)
        GModelFit.showsettings.plain = false
        logio = stdout
    else
        if isfile(logfile)  &&  !overwrite
            error("Logfile: $logfile already exists.")
        end
        logio = open(logfile, "w")
        GModelFit.showsettings.plain = true
    end
    println(logio, "Timestamp: ", now())
    println(logio, "Using $T recipe with options:")
    show(logio, "text/plain", recipe.options)
    println(logio)
    println(logio)

    pspecs = [PreparedSpectrum(source, logio, id=id,
                               cosmology=get_cosmology(T),
                               dered=get_MW_deredd_function(T))
              for id in 1:length(source.specs)]
    if length(pspecs) == 1
        state = State(logfile, logio, pspecs[1], Model(pspecs[1].domain), OrderedDict{Symbol, Any}())
    else
        error("Multi spectrum fit is not yet supported")
    end
    bestfit, fitstats = analyze(recipe, state)
    reduced = reduce(recipe, state)

    out = Results(timestamp,
                  time() - starttime,
                  source, state.pspec,
                  bestfit, fitstats, reduced)

    println(logio, "\nTotal elapsed time: $(out.elapsed) s")
    if !isnothing(logfile)
        close(logio)
        GModelFit.showsettings.plain = false
    end

    return out
end

include("SpectralLines.jl")
include("DefaultRecipe.jl")
include("viewer.jl")
include("gnuplot.jl")

end  # module
