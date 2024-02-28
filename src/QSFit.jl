module QSFit

export RRef, AbstractRecipe, analyze


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

import Dierckx  # `import` (rather than `using`) to avoid conflicts with evaluate()

import Base.show


include("ccm_unred.jl")
include("components/powerlaw.jl")
include("components/sbpl.jl")
include("components/cutoff_powerlaw.jl")
include("components/hostgalaxy.jl")
include("components/ironopt.jl")
include("components/ironuv.jl")
include("components/balmercont.jl")

abstract type AbstractSpecLineComp <: AbstractComponent end
include("components/SpecLineGauss.jl")
include("components/SpecLineLorentz.jl")
include("components/SpecLineVoigt.jl")

include("utils.jl")
include("convolutions.jl")
include("Spectrum.jl")


qsfit_data() = artifact"qsfit_data"


# ====================================================================
mutable struct State
    logfile::Union{Nothing, String}
    logio::Union{IOStream, Base.TTY}
    spec::AbstractSpectrum
    data::Union{Nothing, GModelFit.Measures{1}}
    model::Union{Nothing, GModelFit.Model}
    user::OrderedDict{Symbol, Any}
end


function update_data!(state::State)
    ii = findall(state.spec.good)
    domain = Domain(state.spec.x[ii])
    state.data = Measures(domain,
                          state.spec.y[ii]   .* state.spec.dered[ii],
                          state.spec.err[ii] .* state.spec.dered[ii])
    state.model = Model(domain)
    return state
end


struct Results
    timestamp::DateTime
    elapsed::Float64
    spec::AbstractSpectrum
    data::Measures{1}
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
get_cosmology(::RRef{<: AbstractRecipe}) = cosmology(h=0.70, OmegaM=0.3)
get_dered_function(::RRef{<: AbstractRecipe}) = ccm_unred
convert_units!(::RRef{<: AbstractRecipe}, spec::Spectrum)          = convert_units!(spec, 1. * u"angstrom", 1e-17 * u"erg" / u"s" / u"cm"^2 / u"angstrom")
convert_units!(::RRef{<: AbstractRecipe}, spec::RestFrameSpectrum) = convert_units!(spec, 1. * u"angstrom", 1e42  * u"erg" / u"s"           / u"angstrom")

function prepare_state!(recipe::RRef{<: AbstractRecipe}, state::State)
    if !isnothing(state.spec.ebv)
        dered = get_dered_function(recipe)
        tmp = dered([1450, 3000, 5100.], state.spec.ebv)
        println(state.logio, "De-reddening factors @ 1450, 3000, 5100 AA: ", tmp)
        deredden!(state.spec, dered)
    end
    if !isnothing(state.spec.z)
        state.spec = RestFrameSpectrum(state.spec, get_cosmology(recipe))
    end
    convert_units!(recipe, state.spec)
    round_unit_scales!(state.spec)
    show(state.logio, state.spec)
    update_data!(state)
    return state
end

reduce(recipe::RRef{<: AbstractRecipe}, state::State) = OrderedDict{Symbol, Any}()

function analyze(recipe::RRef{T}, spec::AbstractSpectrum; logfile=nothing, overwrite=false) where T <: AbstractRecipe
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

    state = State(logfile, logio, deepcopy(spec), nothing, nothing, OrderedDict{Symbol, Any}())
    prepare_state!(recipe, state)
    bestfit, fitstats = analyze(recipe, state)
    reduced = reduce(recipe, state)

    out = Results(timestamp,
                  time() - starttime,
                  state.spec, state.data,
                  bestfit, fitstats, reduced)

    println(logio, "\nTotal elapsed time: $(out.elapsed) s")
    if !isnothing(logfile)
        close(logio)
        GModelFit.showsettings.plain = false
    end

    return out
end

include("SpectralLines.jl")
include("recipes/LineFitRecipes.jl")
include("recipes/QSORecipes.jl")
include("viewer.jl")
include("gnuplot.jl")

end  # module
