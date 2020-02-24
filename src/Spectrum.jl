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
    meta::Dict{Symbol, Any}
    function Spectrum(λ::Vector{T}, flux::Vector{T}, err::Vector{T},
                       good::Union{Nothing, Vector{Bool}}=nothing; label="") where T <: AbstractFloat
        if good == nothing
            good = fill(true, length(λ))
        end
        @assert length(λ) == length(flux) == length(err) == length(good)
        @assert issorted(λ)
        igood = findall(good)
        new(label, λ, flux, err, good, Dict{Symbol, Any}())
    end
end

Spectrum(wave::Vector{Quantity}, flux::Vector{Quantity}, err::Vector{Quantity},
          good=nothing; label="") =
    Spectrum(getproperty.(uconvert.(Ref(unit_λ())   , wave), :val),
              getproperty.(uconvert.(Ref(unit_flux()), flux), :val),
              getproperty.(uconvert.(Ref(unit_flux()), err ), :val),
              good, label=label)

goodfraction(d::Spectrum) = length(findall(d.good)) / length(d.good)

function Spectrum(source::Symbol, file::AbstractString)
    if source == :SDSS_DR10
        f = FITS(file)
        λ = 10 .^read(f[2], "loglam")
        flux = float.(read(f[2], "flux"))
        ivar = float.(read(f[2], "ivar"))
        mask = read(f[2], "and_mask")
        close(f)
    
        ndrop = 100
        λ    =    λ[ndrop+1:end-ndrop]
        flux = flux[ndrop+1:end-ndrop]
        ivar = ivar[ndrop+1:end-ndrop]
        mask = mask[ndrop+1:end-ndrop]

        ii = sortperm(λ)
        λ    =    λ[ii]
        flux = flux[ii]
        ivar = ivar[ii]
        mask = mask[ii]

        good = convert(Vector{Bool}, ((mask .== 0)  .&
                                      (ivar .> 0)   .&
                                      (flux .> 0)))

        out = Spectrum(λ, flux, sqrt.(1 ./ ivar), good, label="SDSS-DR10: " * file)
        return out
    elseif source == :ASCII
        λ    = Vector{Float64}()
        flux = Vector{Float64}()
        unc  = Vector{Float64}()
        for l in readlines(file)
            l = strip(strip(l))
            (l[1] == '#')  &&  continue
            s = string.(split(l, " ", keepempty=false))
            @assert length(s) >= 3
            push!(λ   , Meta.parse(s[1]))
            push!(flux, Meta.parse(s[2]))
            push!(unc , Meta.parse(s[3]))
        end
        good = fill(true, length(λ))
        out = Spectrum(λ, flux, unc, good, label="ASCII: " * file)
        return out
    else
        error("Unexpected source: $source")
    end
end
