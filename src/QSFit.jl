module QSFit

export QSO, Spectrum, DefaultRecipe, add_spec!, fit, multiepoch_fit

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


abstract type DefaultRecipe end

struct QSO{T <: DefaultRecipe}
    name::String
    z::Float64
    mw_ebv::Float64
    cosmo::Cosmology.AbstractCosmology
    flux2lum::Float64
    log::IO
    domain::Vector{GFit.Domain_1D}
    data::Vector{GFit.Measures_1D}
    line_names::Vector{OrderedDict{Symbol, Symbol}}
    line_comps::Vector{OrderedDict{Symbol, AbstractComponent}}
    options::OrderedDict{Symbol, Any}

    function QSO{T}(name, z; ebv=0., logfile="", cosmo=default_cosmology())  where T <: DefaultRecipe
        @assert z > 0
        @assert ebv >= 0
        ld = uconvert(u"cm", luminosity_dist(cosmo, float(z)))
        flux2lum = 4pi * ld^2 * (scale_flux() * unit_flux()) / (scale_lum() * unit_lum())
        log = stdout
        (logfile != "")  &&  (log = open(logfile, "w"))
        return new{T}(string(name), float(z), float(ebv), cosmo, flux2lum, log,
                      Vector{GFit.Domain_1D}(), Vector{GFit.Measures_1D}(),
                      Vector{OrderedDict{Symbol, AbstractComponent}}(),
                      Vector{OrderedDict{Symbol, AbstractComponent}}(),
                      default_options(T))
    end
end


abstract type AbstractSpectralLine end

struct BroadBaseLine <: AbstractSpectralLine
    name::Symbol
    λ::Float64
end

struct BroadLine <: AbstractSpectralLine
    name::Symbol
    λ::Float64
end

struct NarrowLine <: AbstractSpectralLine
    name::Symbol
    λ::Float64
end

struct CombinedLine <: AbstractSpectralLine
    name::Symbol
    λ::Float64
end

struct AbsorptionLine <: AbstractSpectralLine
    name::Symbol
    λ::Float64
end

struct UnkLine <: AbstractSpectralLine
end


function line_components_and_groups(source::QSO{T}) where T
    comps = OrderedDict{Symbol, AbstractComponent}()
    groups = OrderedDict{Symbol, Symbol}()
    for line in known_spectral_lines(T)
        ltype = string(typeof(line))
        (ltype[1:6] == "QSFit.")  &&  (ltype = ltype[7:end])
        ltype = Symbol(ltype)
        for (lname, lcomp) in line_components(T, line)
            (lname == :Ha_base)       &&  !source.options[:use_broad_Ha_base]  &&  continue
            (lname == :OIII_5007_bw)  &&  !source.options[:use_OIII_5007_bw ]  &&  continue
            comps[lname] = lcomp
            if (ltype == :CombinedLine)
                if string(lname)[1:3] == "br_"
                    groups[lname] = :BroadLine
                elseif string(lname)[1:3] == "na_"
                    groups[lname] = :NarrowLine
                else
                    groups[lname] = Symbol(ltype)
                end
            else
                groups[lname] = Symbol(ltype)
            end
        end
    end
    return (groups, comps)
end


function add_spec!(source::QSO, data::Spectrum)
    println(source.log, "New spectrum: " * data.label)
    println(source.log, "  good fraction:: ", goodfraction(data))
    if goodfraction(data) < 0.5
        error("Good fraction < 0.5")
    end
    println(source.log, "  resolution: ~", @sprintf("%.4g", data.resolution), " km / s")

    λ = data.λ ./ (1 + source.z)
    data.good[findall(λ .< source.options[:wavelength_range][1])] .= false
    data.good[findall(λ .> source.options[:wavelength_range][2])] .= false

    #= Emission line are localized features whose parameter can be
    reliably estimated only if there are sufficient samples to
    constrain the corresponding parameters.  If data coverage is not
    sufficient the component should not be added to the model, and
    corresponding spectral samples should be ignored to avoid
    worsening the fit due to missing model components. =#
    println(source.log, "Good samples before line coverage filter: ", length(findall(data.good)))
    line_names, line_comps = line_components_and_groups(source)
    for (lname, comp) in line_comps
        (λmin, λmax, coverage) = line_coverage(λ .* data.good, data.resolution, comp.center.val, comp.fwhm.val)
        println(source.log, "Line $lname has coverage: $coverage")
        if coverage < source.options[:line_minimum_coverage]
            println(source.log, "  neglecting line: ", lname, "(", λmin, " < λ < ", λmax)
            ii = findall(λmin .<= λ .< λmax)
            data.good[ii] .= false
            delete!(line_names, lname)
            delete!(line_comps, lname)
        end
    end
    println(source.log, "Good samples after line coverage filter: ", length(findall(data.good)))

    # De-reddening
    dered = ccm_unred([1450, 3000, 5100.], source.mw_ebv)
    println(source.log, "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.mw_ebv)

    ii = findall(data.good)
    dom = Domain(data.λ[ii] ./ (1 + source.z))
    lum = Measures(data.flux[ii] .* dered[ii] .* source.flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* source.flux2lum .* (1 + source.z))
    lum.meta[:label] = data.label

    push!(source.domain, dom)
    push!(source.data, lum)
    push!(source.line_names, line_names)
    push!(source.line_comps, line_comps)
end


function populate_metadata!(source, model)
    for id in 1:length(model.preds)
        model.preds[id].meta[:label] = source.name * ", z=" * string(source.z) * ", E(B-V)=" * string(source.mw_ebv)
        model.preds[id].meta[:label_x] = "Rest frame wavelength"
        model.preds[id].meta[:unit_x]  = string(QSFit.unit_λ())
        model.preds[id].meta[:label_y] = "Lum. density"
        model.preds[id].meta[:unit_y]  = string(QSFit.unit_lum_density())
        model.preds[id].meta[:scale_y] = string(QSFit.scale_lum())
    end
end


include("DefaultRecipe.jl")
include("DefaultRecipe_multiepoch.jl")

end  # module
