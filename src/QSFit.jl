module QSFit

export add_spec!, close_log

import GModelFit: Domain, CompEval,
    Parameter, AbstractComponent, prepare!, evaluate!

using CMPFit, GModelFit, SortMerge
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, Printf, DataStructures
using Unitful, UnitfulAstro
using Dates
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
include("components/SpecLineGauss.jl")
include("components/SpecLineLorentz.jl")
include("components/SpecLineVoigt.jl")
include("utils.jl")
include("convolutions.jl")
include("Spectrum.jl")
include("SpectralLines.jl")


qsfit_data() = artifact"qsfit_data"


struct Source
    name::String
    z::Float64
    mw_ebv::Float64
    specs::Vector{Spectrum}
    function Source(name, z; ebv=0.)
        @assert z >= 0
        @assert ebv >= 0
        return new(string(name), float(z), float(ebv), Vector{Spectrum}())
    end
end

function show(io::IO, source::Source)
    println(io, "Source: ", source.name, ", z=", source.z, ", E(B-V) from Milky Way: ", source.mw_ebv)
    for i in 1:length(source.specs)
        print(io, "    ")
        show(io, source.specs[i])
    end
end

add_spec!(source::Source, spec::Spectrum) =
    push!(source.specs, spec)


abstract type AbstractRecipe end

struct RRef{T <: AbstractRecipe}
    options::OrderedDict{Symbol, Any}
end

function RRef(::Type{T}; kws...) where T <: AbstractRecipe
    out = RRef{T}(OrderedDict{Symbol, Any}())
    set_default_options!(out)
    for (k, v) in kws
        out.options[Symbol(k)] = v
    end
    return out
end


struct PreparedSpectrum
    resolution::Float64
    domain::GModelFit.Domain{1}
    data::GModelFit.Measures{1}
    lcs::OrderedDict{Symbol, EmLineComponent}
    flux2lum::Float64
end


mutable struct State
    logfile::Union{Nothing, String}
    logio::Union{IOStream, Base.TTY}
    pspec::Union{Nothing, PreparedSpectrum}
    model::Union{Nothing, GModelFit.Model}
end


struct Results
    elapsed::Float64
    logfile::Union{Nothing, String}
    pspec::Union{Nothing, PreparedSpectrum}
    model::Union{Nothing, GModelFit.Model}
    bestfit::GModelFit.ModelSnapshot
    fitstats::GModelFit.FitStats
    reduced::OrderedDict{Symbol, Any}
end


function show(io::IO, res::Results)
    show(io, res.fitstats)
    println(io, "\nTotal elapsed time: $(res.elapsed) s")
end


function analyze(recipe::RRef{T}, source::Source; logfile=nothing, overwrite=false) where T <: AbstractRecipe
    starttime = time()
    if isnothing(logfile)
        GModelFit.showsettings.plain = false
        logio = stdout
    else
        if isfile(logfile)  &&  !overwrite
            error("Logfile: $logfile already exists.")
        end
        logio = open(logfile, "w")
        println(logio, "Timestamp: ", now())
        GModelFit.showsettings.plain = true
    end

    state = State(logfile, logio, nothing, nothing)
    state.pspec = PreparedSpectrum(recipe, state, source, id=1)
    state.model = Model(state.pspec.domain)
    bestfit, fitstats = analyze(recipe, state)
    reduced = reduce(recipe, state)
    
    if !isnothing(logfile)
        close(logio)
        GModelFit.showsettings.plain = false
    end

    out = Results(time() - starttime, state.logfile, state.pspec, state.model,
                  bestfit, fitstats, reduced)
    return out
end


include("DefaultRecipe.jl")

# Use DefaultRecipe when no explicit recipe is provided
analyze(source::Source; kws...) = analyze(RRef(DefaultRecipe), source; kws...)


# TODO include("viewer.jl")
include("gnuplot.jl")
# TODO include("interactive_guess.jl")

end  # module
