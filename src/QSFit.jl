module QSFit

export interpol, QSO, Spectrum, add_spec!, fit, multiepoch_fit

import GFit: Domain_1D, CompEval,
    Parameter, AbstractComponent, compeval_cdata, compeval_array, evaluate, fit!

using CMPFit, GFit
using Pkg, Pkg.Artifacts
using Statistics, DelimitedFiles, DataFrames, FFTW, Interpolations, QuadGK, Printf, DataStructures
using Unitful, UnitfulAstro

# GFit.showsettings.showfixed = true

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
include("utils.jl")
include("Spectrum.jl")


abstract type AbstractRecipe end

function default_options(::Type{T}) where T <: AbstractRecipe
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict{Symbol, Float64}(:default => 0.6)
    out[:skip_lines] = Vector{Symbol}()
    return out
end

struct QSO{T <: AbstractRecipe}
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

    function QSO{T}(name, z; ebv=0., logfile="", cosmo=default_cosmology())  where T <: AbstractRecipe
        @assert z > 0
        @assert ebv >= 0
        ld = uconvert(u"cm", luminosity_dist(cosmo, float(z)))
        flux2lum = 4pi * ld^2 * (scale_flux() * unit_flux()) / (scale_lum() * unit_lum())
        if logfile != ""
            log = open(logfile, "w")
            GFit.showsettings.plain = true
        else
            log = stdout
            GFit.showsettings.plain = false
        end
        return new{T}(string(name), float(z), float(ebv), cosmo, flux2lum, log,
                      Vector{GFit.Domain_1D}(), Vector{GFit.Measures_1D}(),
                      Vector{OrderedDict{Symbol, AbstractComponent}}(),
                      Vector{OrderedDict{Symbol, AbstractComponent}}(),
                      default_options(T))
    end
end


function close_log(source::QSO)
    if source.log != stdout
        close(source.log)
        GFit.showsettings.plain = false
        # source.log = stdout
    end
end

abstract type AbstractSpectralLine end

struct BroadBaseLine <: AbstractSpectralLine
    λ::Float64
end

struct BroadLine <: AbstractSpectralLine
    λ::Float64
end

struct NarrowLine <: AbstractSpectralLine
    λ::Float64
end

struct ComboBroadLine <: AbstractSpectralLine
    λ::Float64
end

struct ComboNarrowLine <: AbstractSpectralLine
    λ::Float64
end

struct CombinedLine <: AbstractSpectralLine
    λ::Float64
end

struct AbsorptionLine <: AbstractSpectralLine
    λ::Float64
end

struct UnkLine <: AbstractSpectralLine
    λ::Float64
end


# Note: line_breakdown must always return a Vector{Tuple{Symbol, <: AbstractSpectralLine}}
line_breakdown(::Type{T}, name::Symbol, line::L) where {T <: AbstractRecipe, L <: AbstractSpectralLine} =
    [(name, line)]

line_group_name(::Type{T}, name::Symbol, line::L) where {T <: AbstractRecipe, L <: AbstractSpectralLine} =
    Symbol(L)

function line_components_and_groups(source::QSO{T}) where T
    comps = OrderedDict{Symbol, AbstractComponent}()
    groups = OrderedDict{Symbol, Symbol}()
    for (lname0, line0) in known_spectral_lines(T)
        (lname0 in source.options[:skip_lines])  &&  continue
        for (lname, line) in line_breakdown(T, lname0, line0)
            lcomp = line_component(T, line)
            ltype = string(typeof(lcomp))
            (ltype[1:6] == "QSFit.")  &&  (ltype = ltype[7:end])
            ltype = Symbol(ltype)
            comps[lname] = lcomp
            groups[lname] = line_group_name(T, lname0, line) # Note: here we use original name and specific line type
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
        (λmin, λmax, coverage) = spectral_coverage(λ .* data.good, data.resolution, comp)
        coverage = round(coverage * 1e3) / 1e3  # keep just 3 significant digits...
        threshold = get(source.options[:min_spectral_coverage], lname, source.options[:min_spectral_coverage][:default])
        print(source.log, "Line $lname coverage: $coverage (threshold: $threshold)")
        if coverage < threshold
            print(source.log, "  neglecting range: $λmin < λ <  $λmax")
            ii = findall(λmin .<= λ .< λmax)
            data.good[ii] .= false
            delete!(line_names, lname)
            delete!(line_comps, lname)
        end
        println()
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
        model.meta[id][:label] = source.name * ", z=" * string(source.z) * ", E(B-V)=" * string(source.mw_ebv)
        model.meta[id][:label_x] = "Rest frame wavelength"
        model.meta[id][:unit_x]  = string(QSFit.unit_λ())
        model.meta[id][:label_y] = "Lum. density"
        model.meta[id][:unit_y]  = string(QSFit.unit_lum_density())
        model.meta[id][:scale_y] = string(QSFit.scale_lum())
    end
end


include("DefaultRecipe.jl")
include("DefaultRecipe_multiepoch.jl")

end  # module
