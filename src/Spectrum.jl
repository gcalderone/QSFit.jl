using FITSIO

unit_λ() = UnitfulAstro.angstrom
unit_flux() = 1.e-17 * UnitfulAstro.erg / UnitfulAstro.s / UnitfulAstro.cm^2
unit_lum() =  1.e42  * UnitfulAstro.erg / UnitfulAstro.s

struct Spectrum
    label::String
    λ::Vector{Float64}
    flux::Vector{Float64}
    err::Vector{Float64}
    good::Vector{Bool}
    resolution::Float64  # km/s
    meta::Dict{Symbol, Any}
    function Spectrum(λ::Vector{T}, flux::Vector{T}, err::Vector{T},
                      good::Union{Nothing, Vector{Bool}}=nothing;
                      label="") where T <: AbstractFloat
        if good == nothing
            good = fill(true, length(λ))
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


Spectrum(λ::Vector{Quantity}, flux::Vector{Quantity}, err::Vector{Quantity}; label="") =
    Spectrum(getproperty.(uconvert.(Ref(unit_λ())   , λ   ), :val),
             getproperty.(uconvert.(Ref(unit_flux()), flux), :val),
             getproperty.(uconvert.(Ref(unit_flux()), err ), :val),
             nothing, label=label)


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
    
    out = Spectrum(λ, flux, sqrt.(1 ./ ivar), good, label="SDSS: " * file)
    return out
end


function Spectrum(::Val{:ASCII}, file::AbstractString; columns=[1,2,3])
    @assert length(columns) >= 2
    if length(columns) < 3
        @warn "Uncertainty is not given: assuming 10% of flux"
    end
    
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
        else
            push!(unc, abs(0.1 * flux[end]))
        end
    end
    good = fill(true, length(λ))
    out = Spectrum(λ, flux, unc, good, label="ASCII: " * file)
    return out
end


goodfraction(d::Spectrum) = length(findall(d.good)) / length(d.good)
