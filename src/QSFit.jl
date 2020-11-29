module QSFit

export QSO, Spectrum, add_spec!, fit!

import GFit: Domain_1D, CompEval,
    Parameter, AbstractComponent, compeval_cdata, compeval_array, evaluate, fit!

using CMPFit, GFit
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, DataFrames, FFTW, Interpolations, QuadGK, Printf, DataStructures
using Unitful, UnitfulAstro

GFit.@with_CMPFit
# GFit.showsettings.showfixed = true

include("utils.jl")
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
include("Spectrum.jl")
#include("spectral_lines.jl")

abstract type AbstractSource end
struct T1_generic <: AbstractSource end

struct QSO{T <: AbstractSource}
    name::String
    z::Float64
    mw_ebv::Float64
    cosmo::Cosmology.AbstractCosmology
    flux2lum::Float64
    log::IO
    domain::Vector{GFit.Domain_1D}
    data::Vector{GFit.Measures_1D}
    broad_lines::Vector{OrderedDict{Symbol, AbstractComponent}}
    narrow_lines::Vector{OrderedDict{Symbol, AbstractComponent}}

    function QSO{T}(name, z; ebv=0., logfile="", cosmo=default_cosmology())  where T <: AbstractSource
        @assert z > 0
        @assert ebv >= 0
        ld = uconvert(u"cm", luminosity_dist(cosmo, float(z)))
        flux2lum = 4pi * ld^2 * unit_flux() / unit_lum()
        log = stdout
        (logfile != "")  &&  (log = open(logfile, "w"))
        return new{T}(string(name), float(z), float(ebv), cosmo, flux2lum, log,
                      Vector{GFit.Domain_1D}(), Vector{GFit.Measures_1D}(),
                      Vector{OrderedDict{Symbol, AbstractComponent}}(),
                      Vector{OrderedDict{Symbol, AbstractComponent}}())
    end
end


function add_spec!(source::QSO, data::Spectrum)
    @assert length(source.data) == 0

    println(source.log, "New spectrum: " * data.label)
    println(source.log, "  good fraction:: ", goodfraction(data))
    if goodfraction(data) < 0.5
        @error "Good fraction < 0.5"
    end
    println(source.log, "  resolution: ~", @sprintf("%.4g", data.resolution), " km / s")

    λ = data.λ ./ (1 + source.z)
    data.good[findall(λ .< options(source)[:wavelength_range][1])] .= false
    data.good[findall(λ .> options(source)[:wavelength_range][2])] .= false

    narrow_lines = OrderedDict{Symbol, AbstractComponent}()
    broad_lines  = OrderedDict{Symbol, AbstractComponent}()

    # Ignore lines on missing data
    println(source.log, "Good samples before line coverage filter: ", length(findall(data.good)))
    for (lname, (ltype, lwave)) in default_known_lines(source)
        if ltype == :Narrow
            comp = default_narrowline(source, lwave)
        elseif ltype == :Broad
            comp = default_broadline(source, lwave)
        else
            @error "Unexpected line type: $ltype"
        end

        (λmin, λmax, coverage) = line_coverage(λ .* data.good, data.resolution, comp.center.val, comp.fwhm.val)
        @info "Line $lname has coverage: $coverage"
        if coverage >= options(source)[:line_minimum_coverage]
            if ltype == :Narrow
                narrow_lines[lname] = comp
            elseif ltype == :Broad
                broad_lines[lname] = comp
            end
        else
            println(source.log, "  neglecting line: ", lname)
            ii = findall(λmin .<= λ .< λmax)
            data.good[ii] .= false
        end
    end
    println(source.log, "Good samples after line coverage filter: ", length(findall(data.good)))

    dered = ccm_unred([1450, 3000, 5100.], source.mw_ebv)
    println(source.log, "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.mw_ebv)

    ii = findall(data.good)
    dom = Domain(data.λ[ii] ./ (1 + source.z))
    lum = Measures(data.flux[ii] .* dered[ii] .* source.flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* source.flux2lum .* (1 + source.z))
    push!(source.domain, dom)
    push!(source.data, lum)
    push!(source.narrow_lines, narrow_lines)
    push!(source.broad_lines , broad_lines)
end

include("default_recipe.jl")

end  # module
