module QSFit

export CRecipe, AbstractRecipe, analyze, @track_recipe

import GModelFit: Domain, CompEval,
    Parameter, AbstractComponent, dependencies, prepare!, evaluate!

import Base: propertynames, getproperty, setproperty!, haskey, keys, values, iterate, getindex, setindex!, delete!, show

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

abstract type AbstractSpecLineComp <: AbstractComponent end
include("components/SpecLineGauss.jl")
include("components/SpecLineAsymmGauss.jl")
include("components/SpecLineLorentz.jl")
include("components/SpecLineVoigt.jl")


# The following is needed to allow build model incrementally. More
# specifically, it allows to evaluiate a model with a SumReducer
# component having no dependency
function evaluate!(::GModelFit.SumReducer, ::AbstractDomain, output)
    output .= 0.
end


include("utils.jl")
include("convolutions.jl")
include("Spectrum.jl")


qsfit_data() = artifact"qsfit_data"

global _track_recipe::Bool = false
function track_recipe(enable::Bool)
    global _track_recipe
    _track_recipe = enable
end
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
    s = join(sort(split(s, "\n")[2:end]), "\n")
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
    data::GModelFit.AbstractMeasures
    bestfit::GModelFit.ModelSnapshot
    fsumm::GModelFit.FitSummary
    post::OrderedDict{Symbol, Any}
end

struct MultiResults
    timestamp::DateTime
    elapsed::Float64
    spec::Vector{Spectrum}
    data::Vector{<: GModelFit.AbstractMeasures}
    bestfit::Vector{GModelFit.ModelSnapshot}
    fsumm::GModelFit.FitSummary
    post::Vector{OrderedDict{Symbol, Any}}
end

function show(io::IO, res::Union{Results, MultiResults})
    show(io, res.fsumm)
    println(io)
end

# ====================================================================
function preprocess_spec!(recipe::CRecipe{<: AbstractRecipe}, spec::Spectrum)
    @track_recipe
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


function spec2data(recipe::CRecipe{<: AbstractRecipe}, spec::Spectrum)
    ii = findall(spec.good)
    dom = Domain(spec.x[ii])
    data = Measures(dom, spec.y[ii], spec.err[ii])
    return data
end


function analyze(_recipe::CRecipe{<: AbstractRecipe}, _spec::Spectrum)
    @track_recipe
    tstart = now();
    println("Timestamp: ", tstart)
    recipe = deepcopy(_recipe)
    display(recipe); println()
    spec = deepcopy(_spec)
    show(spec)

    preprocess_spec!(recipe, spec)
    data = spec2data(recipe, spec)
    bestfit, fsumm = analyze(recipe, data)
    post = postanalysis(recipe, bestfit)

    out = Results(tstart,
                  Dates.value(convert(Millisecond, now() - tstart)) / 1000.,
                  spec, data, bestfit, fsumm, post)
    println("\nTotal elapsed time: $(out.elapsed) s")
    return out
end

function analyze(_recipe::CRecipe{<: AbstractRecipe}, _specs::Vector{Spectrum})
    @track_recipe
    tstart = now();
    println("Timestamp: ", tstart)
    recipe = deepcopy(_recipe)
    display(recipe); println()
    specs = deepcopy(_specs)
    for i in 1:length(specs)
        show(specs[i])
    end

    preprocess_spec!.(Ref(recipe), specs)
    data = spec2data.(Ref(recipe), specs)
    bestfit, fsumm = analyze(recipe, data)
    post = [postanalysis(recipe, b) for b in bestfit]

    out = MultiResults(tstart,
                       Dates.value(convert(Millisecond, now() - tstart)) / 1000.,
                       specs, data, bestfit, fsumm, post)
    println("\nTotal elapsed time: $(out.elapsed) s")
    return out
end

postanalysis(recipe::CRecipe{<: AbstractRecipe}, bestfit::GModelFit.ModelSnapshot) = OrderedDict{Symbol, Any}()

include("SpectralLines.jl")
include("recipes/LineFitRecipes.jl")
include("recipes/QSORecipes.jl")
include("viewer.jl")
include("serialize.jl")
include("gnuplot.jl")

end  # module
