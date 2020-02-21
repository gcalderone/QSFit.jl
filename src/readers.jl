using FITSIO

unit_λ() = UnitfulAstro.angstrom
unit_flux() = 1.e-17 * UnitfulAstro.erg / UnitfulAstro.s / UnitfulAstro.cm^2
unit_lum() =  1.e42  * UnitfulAstro.erg / UnitfulAstro.s

struct QSFitData
    label::String
    λ::Vector{Float64}
    flux::Vector{Float64}
    err::Vector{Float64}
    good::Vector{Bool}
    goodfraction::Float64
    median_flux::Float64
    median_err::Float64
    meta::Dict{Symbol, Any}
    function QSFitData(λ::Vector{T}, flux::Vector{T}, err::Vector{T},
                       good::Union{Nothing, Vector{Bool}}=nothing; label="") where T <: AbstractFloat
        if good == nothing
            good = fill(true, length(λ))
        end
        @assert length(λ) == length(flux) == length(err) == length(good)
        @assert issorted(λ)
        igood = findall(good)
        new(label, λ, flux, err, good, length(igood)/length(λ),
            median(flux[igood]), median(err[igood]), Dict{Symbol, Any}())
    end
end

QSFitData(wave::Vector{Quantity}, flux::Vector{Quantity}, err::Vector{Quantity}, good=nothing; label="") =
    QSFitData(getproperty.(uconvert.(Ref(unit_λ())   , wave), :val),
              getproperty.(uconvert.(Ref(unit_flux()), flux), :val),
              getproperty.(uconvert.(Ref(unit_flux()), err ), :val),
              good, label=label)



function read_sdss_dr10(file::AbstractString)
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

    d = QSFitData(λ, flux, sqrt.(1 ./ ivar), good, label="SDSS-DR10: " * file)
    return d
end
