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
    id::Int
    domain::GFit.Domain{1}
    data::GFit.Measures{1}
    lcs::OrderedDict{Symbol, EmLineComponent}
end


abstract type JobState{T <: AbstractRecipe} <: Job{T} end
struct       cJobState{T <: AbstractRecipe} <: JobState{T}
    @copy_fields(cJob)
    source::Source
    pspec::StdSpectrum
    model::GFit.Model
end

function JobState{T}(source::Source, job::Job{T}) where T <: AbstractRecipe
    pspec = StdSpectrum(job, source)
    cJobState{T}(getfield.(Ref(job), fieldnames(typeof(job)))..., source, pspec, Model(pspec.domain))
end

# TODO struct JobStateMulti{T <: AbstractRecipe}
# TODO     @copy_fields(Job)
# TODO     source::Source
# TODO     pspecs::Vector{StdSpectrum}
# TODO     models::GFit.MultiModel
# TODO end


abstract type JobResults{T <: AbstractRecipe} <: JobState{T} end
struct       cJobResults{T} <: JobResults{T}
    @copy_fields(cJobState)
    fitres::GFit.FitResult
    EW::OrderedDict{Symbol, Float64}
end

JobResults(job::JobState{T}, fitres::GFit.FitResult) where T <: AbstractRecipe =
    cJobResults{T}(getfield.(Ref(job), fieldnames(typeof(job)))..., fitres, estimate_line_EWs(source, pspec, model))

run(source::Source, job::Job{T}) where T <: AbstractRecipe =
    run(JobState{T}(source, job))


include("DefaultRecipe.jl")
# TODO include("reduce.jl")
# TODO include("viewer.jl")
# TODO include("gnuplot.jl")
# TODO include("interactive_guess.jl")

end  # module
