using FITSIO

unit_λ() = u"angstrom"
unit_flux() = u"erg" / u"s" / u"cm"^2
unit_lum() =  u"erg" / u"s"

unit_λ() = u"angstrom"
unit_flux_density() = unit_flux() / unit_λ()
unit_lum_density() = unit_lum() / unit_λ()

scale_λ() = 1.
scale_flux() = 1.e-17
scale_lum() =  1.e42

struct Spectrum
    label::String
    λ::Vector{Float64}
    flux::Vector{Float64}
    err::Vector{Float64}
    good::Vector{Bool}
    resolution::Float64  # km/s
    meta::Dict{Symbol, Any}
    function Spectrum(λ::Vector{T}, flux::Vector{T}, err::Vector{T};
                      good::Union{Nothing, Vector{Bool}}=nothing,
                      label="") where T <: AbstractFloat
        if good == nothing
            good = fill(true, length(λ))
        end
        if length(err) == 0
            @warn "Uncertainty is not given: assuming 10% of flux"
            err = 0.1 .* flux
        end
        @assert length(λ) == length(flux) == length(err) == length(good)

        # Ensure wavelengths are sorted
        ii = sortperm(λ)
        λ = λ[ii]

        # Estimate average spectral resolution in km/s
        resolution = median((λ[2:end] .- λ[1:end-1]) ./ ((λ[2:end] .+ λ[1:end-1]) ./ 2)) * 3e5

        return new(label, λ, flux[ii], err[ii], good, resolution,
                   Dict{Symbol, Any}())
    end
end


Spectrum(λ::Vector{T1}, flux::Vector{T2}; kw...) where {T1 <: Quantity, T2 <: Quantity} = Spectrum(λ, flux, flux[[]]; kw...)
function Spectrum(λ::Vector{T1}, flux::Vector{T2}, err::Vector{T2}; kw...) where {T1 <: Quantity, T2 <: Quantity}
    Spectrum(getproperty.(uconvert.(Ref(unit_λ())              , λ   ), :val) ./ scale_λ(),
             getproperty.(uconvert.(Ref(unit_flux_density()), flux), :val) ./ scale_flux(),
             getproperty.(uconvert.(Ref(unit_flux_density()), err ), :val) ./ scale_flux();
             kw...)
end


function Spectrum(::Val{:SDSS_DR10}, file::AbstractString; ndrop=100)
    f = FITS(file)
    λ = 10 .^read(f[2], "loglam")
    flux = float.(read(f[2], "flux"))
    ivar = float.(read(f[2], "ivar"))
    mask = read(f[2], "and_mask")
    close(f)
    
    ii = sortperm(λ)
    λ    =    λ[ii]
    flux = flux[ii]
    ivar = ivar[ii]
    mask = mask[ii]

    good = convert(Vector{Bool}, ((mask .== 0)  .&
                                  (ivar .> 0)   .&
                                  (flux .> 0)))
    if ndrop > 0
        good[1:ndrop] .= false
        good[end-ndrop+1:end] .= false
    end
    
    out = Spectrum(λ, flux, sqrt.(1 ./ ivar), good=good, label="SDSS: " * file)
    return out
end


function Spectrum(::Val{:ASCII}, file::AbstractString; columns=[1,2,3])
    @assert length(columns) >= 2
    
    λ    = Vector{Float64}()
    flux = Vector{Float64}()
    unc  = Vector{Float64}()
    for l in readlines(file)
        l = strip(strip(l))
        (l[1] == '#')  &&  continue
        s = string.(split(l, keepempty=false))
        @assert length(s) >= maximum(columns)
        push!(λ   , Meta.parse(s[columns[1]]))
        push!(flux, Meta.parse(s[columns[2]]))
        if length(columns) == 3
            push!(unc, Meta.parse(s[columns[3]]))
        end
    end
    good = fill(true, length(λ))
    out = Spectrum(λ, flux, unc, label="ASCII: " * file)
    return out
end


goodfraction(d::Spectrum) = length(findall(d.good)) / length(d.good)
