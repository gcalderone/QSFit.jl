module QSFit

export add_spec!, close_log

import GFit: Domain, CompEval,
    Parameter, AbstractComponent, prepare!, evaluate!, fit!

using CMPFit, GFit, SortMerge
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, Dierckx, Printf, DataStructures
using Unitful, UnitfulAstro
using Dates
using SpecialFunctions
using Gnuplot
using TextParse
using Cosmology

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


version() = v"0.1.0"
qsfit_data() = artifact"qsfit_data"


struct Source
    name::String
    z::Float64
    mw_ebv::Float64
    specs::Vector{Spectrum}
    function Source(name, z; ebv=0.)
        @assert z > 0
        @assert ebv >= 0
        return new(string(name), float(z), float(ebv), Vector{Spectrum}())
    end
end

add_spec!(source::Source, spec::Spectrum) =
    push!(source.specs, spec)


abstract type AbstractRecipe end


abstract type Job{T <: AbstractRecipe} end
struct       cJob{T <: AbstractRecipe} <: Job{T}
    options::OrderedDict{Symbol, Any}
    logfile::Union{Nothing, String}
    logio::Union{IOStream, Base.TTY}
end

function Job{T}(; logfile=nothing) where T <: AbstractRecipe
    if isnothing(logfile)
        GFit.showsettings.plain = false
        logio = stdout
    else
        if isfile(logfile)
            @warn "Logfile: $logfile already exists, overwriting..."
        end
        logio = open(logfile, "w")
        println(logio, "Timestamp: ", now())
        GFit.showsettings.plain = true
    end
    return cJob{T}(Options(T), logfile, logio)
end


function close_log(job::Job{T}) where T <: AbstractRecipe
    if !isnothing(job.logfile)
        close(job.logio)
        GFit.showsettings.plain = false
    end
    nothing
end


struct StdSpectrum{T <: AbstractRecipe}
    resolution::Float64
    domain::GFit.Domain{1}
    data::GFit.Measures{1}
    lcs::OrderedDict{Symbol, EmLineComponent}
end


abstract type JobState{T <: AbstractRecipe} <: Job{T} end
struct       cJobState{T <: AbstractRecipe} <: JobState{T}
    @copy_fields(cJob)
    pspec::StdSpectrum
    model::GFit.Model
end

function JobState{T}(source::Source, job::Job{T}; id=1) where T <: AbstractRecipe
    pspec = StdSpectrum(job, source, id=id)
    cJobState{T}(getfield.(Ref(job), fieldnames(typeof(job)))..., pspec, Model(pspec.domain))
end


abstract type JobMultiState{T <: AbstractRecipe} <: Job{T} end
struct       cJobMultiState{T <: AbstractRecipe} <: JobMultiState{T}
    @copy_fields(cJob)
    pspecs::Vector{StdSpectrum}
    models::GFit.MultiModel
end

function JobMultiState{T}(source::Source, job::Job{T}) where T <: AbstractRecipe
    pspecs = [StdSpectrum(job, source, id=id) for id in 1:length(source.specs)]
    cJobMultiState{T}(getfield.(Ref(job), fieldnames(typeof(job)))..., pspecs, MultiModel())
end


function run(source::Source, job::Job{T}) where T <: AbstractRecipe
    if length(source.specs) == 1
        return run(JobState{T}(source, job))
    end
    return run(JobMultiState{T}(source, job))
end

abstract type JobResults{T <: AbstractRecipe} <: JobState{T} end
struct       cJobResults{T} <: JobResults{T}
    @copy_fields(cJobState)
    fitres::GFit.FitResult
    elapsed::Float64
    reduced::OrderedDict{Symbol, Any}
end

JobResults(job::JobState{T}, fitres::GFit.FitResult, elapsed::Float64) where T <: AbstractRecipe =
    cJobResults{T}(getfield.(Ref(job), fieldnames(typeof(job)))..., fitres, elapsed, OrderedDict{Symbol, Any}())


include("DefaultRecipe.jl")
include("reduce.jl")
# TODO include("viewer.jl")
# TODO include("gnuplot.jl")
# TODO include("interactive_guess.jl")

end  # module
