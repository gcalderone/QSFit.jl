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
using TextParse, DataFrames

include("cosmology.jl")
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

struct Job{T <: AbstractRecipe}
    options::OrderedDict{Symbol, Any}
    cosmo::Cosmology.AbstractCosmology
    logfile::Union{Nothing, String}
    logio::Union{IOStream, Base.TTY}
end


function Job{T}(;
                logfile=nothing,
                cosmo=default_cosmology()) where T <: AbstractRecipe
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
    return Job{T}(Options(T), cosmo,
                  logfile, logio)
end


function close_log(job::Job{T}) where T <: AbstractRecipe
    if !isnothing(job.logfile)
        close(job.logio)
        GFit.showsettings.plain = false
    end
    nothing
end


struct PreparedSpectrum
    id::Int
    domain::GFit.Domain{1}
    data::GFit.Measures{1}
    flux2lum::Float64
    lcs::OrderedDict{Symbol, EmLineComponent}
end


struct JobState{T <: AbstractRecipe}
    @copy_fields(Job)
    source::Source
    pspec::PreparedSpectrum
    model::GFit.Model
end

function JobState{T}(source::Source, job::Job{T}) where T <: AbstractRecipe
    pspec = PreparedSpectrum(job, source)
    JobState{T}(getfield.(Ref(job), fieldnames(typeof(job)))..., source, pspec, Model(pspec.domain))
end

# TODO struct JobStateMulti{T <: AbstractRecipe}
# TODO     @copy_fields(Job)
# TODO     source::Source
# TODO     pspecs::Vector{PreparedSpectrum}
# TODO     models::GFit.MultiModel
# TODO end


include("DefaultRecipe.jl")
# TODO include("reduce.jl")
# TODO include("viewer.jl")
# TODO include("gnuplot.jl")
# TODO include("interactive_guess.jl")

end  # module
