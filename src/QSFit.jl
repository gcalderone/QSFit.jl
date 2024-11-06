module QSFit

export CRecipe, AbstractRecipe, analyze, @track_recipe

import GModelFit: Domain, CompEval, Residuals,
    Parameter, AbstractComponent, dependencies, prepare!, evaluate!

import Base: propertynames, getproperty, setproperty!, show

using CMPFit, GModelFit, SortMerge
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, Printf, DataStructures
using Unitful, UnitfulAstro
using FITSIO
using Dates
using InteractiveUtils
using SpecialFunctions
using Polyester
using Gnuplot
using TextParse
using Cosmology
using DustExtinction

import Dierckx  # `import` (rather than `using`) to avoid conflicts with evaluate()

include("components/powerlaw.jl")
include("components/sbpl.jl")
include("components/cutoff_powerlaw.jl")
include("components/hostgalaxy.jl")
include("components/ironopt.jl")
include("components/ironuv.jl")
include("components/balmercont.jl")
include("components/gaussconv.jl")
include("components/interpolator.jl")

abstract type AbstractSpecLineComp <: AbstractComponent end
include("components/SpecLineGauss.jl")
include("components/SpecLineAsymmGauss.jl")
include("components/SpecLineLorentz.jl")
include("components/SpecLineVoigt.jl")

include("utils.jl")
include("convolutions.jl")
include("Spectrum.jl")


qsfit_data() = artifact"qsfit_data"

global _track_recipe = false
macro track_recipe()
    return :(QSFit._track_recipe  &&  printstyled("* ", stacktrace()[1], "\n", color=:light_black))
end

# ====================================================================
abstract type AbstractRecipe end

struct CRecipe{T <: AbstractRecipe}
    dict::OrderedDict{Symbol, Any}
    function CRecipe{T}(; kws...) where T <: AbstractRecipe
        @track_recipe
        out = new{T}(OrderedDict{Symbol, Any}())
        init_recipe!(out)
        for (key, value) in kws  # set options provided as keywords
            setproperty!(out, key, value)
        end
        return out
    end
end

propertynames(recipe::CRecipe) = collect(keys(getfield(recipe, :dict)))
getproperty(recipe::CRecipe, key::Symbol) = getfield(recipe, :dict)[key]
setproperty!(recipe::CRecipe, key::Symbol, value) = getfield(recipe, :dict)[key] = value
function show(io::IO, recipe::CRecipe)
    @track_recipe
    println(io, typeof(recipe))
    tmp = IOBuffer()
    show(tmp, "text/plain", getfield(recipe, :dict))
    s = String(take!(tmp))
    s = join(split(s, "\n")[2:end], "\n")
    println(io, s)
end

function init_recipe!(recipe::CRecipe{T}) where T <: AbstractRecipe
    @track_recipe
    recipe.cosmology = cosmology(h=0.70, OmegaM=0.3)
    recipe.redshift = missing
    recipe.extlaw = OD94(Rv=3.1)
    recipe.Av = missing
    recipe.unit_x    = 1. * u"angstrom"
    recipe.unit_flux = 1e-17 * u"erg" / u"s" / u"cm"^2 / u"angstrom"
    recipe.unit_lum  = 1e42  * u"erg" / u"s"           / u"angstrom"
end


# ====================================================================
struct Results
    timestamp::DateTime
    elapsed::Float64
    spec::Spectrum
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
function preprocess_spec!(recipe::CRecipe{<: AbstractRecipe}, spec::Spectrum)
    @track_recipe
    show(spec)
    if !ismissing(recipe.Av)
        @printf "      Av: %8.3f  (%s)\n" recipe.Av join(split(string(recipe.extlaw), "\n"), ",")
        deredden!(spec, recipe.extlaw, recipe.Av)
    end
    if ismissing(recipe.redshift)
        convert_units!(spec, recipe.unit_x, recipe.unit_flux)
    else
        @printf "       z: %8.3f  (%s)\n" recipe.redshift string(recipe.cosmology)
        torestframe!(spec, recipe.cosmology, recipe.redshift)
        convert_units!(spec, recipe.unit_x, recipe.unit_lum)
    end
    round_unit_scales!(spec)
end


reduce(recipe::CRecipe{<: AbstractRecipe}, resid::Residuals) = OrderedDict{Symbol, Any}()

function analyze(_recipe::CRecipe{T}, _spec::Spectrum) where T <: AbstractRecipe
    @track_recipe
    timestamp = now()
    starttime = time()

    recipe = deepcopy(_recipe)
    spec = deepcopy(_spec)
    println("Timestamp: ", timestamp)
    display(recipe)
    println()

    preprocess_spec!(recipe, spec)

    # Create GModelFit objects
    ii = findall(spec.good)
    domain = Domain(spec.x[ii])
    data = Measures(domain, spec.y[ii], spec.err[ii])
    meval = GModelFit.ModelEval(Model(), domain)
    resid = Residuals(meval, data, GModelFit.cmpfit())

    bestfit, stats = analyze(recipe, spec, resid)
    reduced = reduce(recipe, resid)

    out = Results(timestamp,
                  time() - starttime,
                  spec, data,
                  bestfit, stats, reduced)

    println("\nTotal elapsed time: $(out.elapsed) s")
    return out
end

include("SpectralLines.jl")
include("recipes/LineFitRecipes.jl")
include("recipes/QSORecipes.jl")
include("viewer.jl")
include("gnuplot.jl")

end  # module
