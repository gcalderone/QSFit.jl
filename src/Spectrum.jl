using FITSIO, DSP

unit_λ() = u"angstrom"
unit_flux() = u"erg" / u"s" / u"cm"^2
unit_lum() =  u"erg" / u"s"

unit_flux_density() = unit_flux() / unit_λ()
unit_lum_density() = unit_lum() / unit_λ()

scale_λ() = 1.
scale_flux() = 1.e-17
scale_lum() = 1.e42

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
                      resolution=NaN,
                      label="") where T <: AbstractFloat
        if good == nothing
            good = fill(true, length(λ))
        end
        if length(err) == 0
            @warn "Uncertainties were not provided: assuming 10% of flux"
            err = 0.1 .* flux
        end
        @assert length(λ) == length(flux) == length(err) == length(good)

        # Ensure wavelengths are sorted
        ii = sortperm(λ)
        λ = λ[ii]

        # Estimate sampling resolution in km/s
        est_res = median((λ[2:end] .- λ[1:end-1]) ./ ((λ[2:end] .+ λ[1:end-1]) ./ 2)) * 3e5
        if isnan(resolution)
            @warn "Resolution is not provided, assuming it is equal to twice the sampling resolution..."
            resolution = 2 * est_res
        end
        @assert est_res < resolution

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

    out = Spectrum(λ, flux, sqrt.(1 ./ ivar), good=good, label=file, resolution=150.)  # TODO: Check resolution is correct
    return out
end


function Spectrum(::Val{:ASCII}, file::AbstractString; columns=[1,2,3], kw...)
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
    out = Spectrum(λ, flux, unc, label=file; kw...)
    return out
end


goodfraction(d::Spectrum) = length(findall(d.good)) / length(d.good)



function instrumental_broadening(λ, flux, σ_kms)
    # In some case the main reducer may evaluate to a single value
    # (rather than a vector)
    (length(flux) == 1)  &&  (return flux)

    #=
    A regular log-λ grid is characterized by:
    log10(λ_i+1) - log10(λ_i)  =  log10(λ_i+1 / λ_i)  =  costant step

    hence:
    λ_i+1 / λ_i  =  cost.
    λ_i+1 / λ_i - 1  =  cost.
    (λ_i+1 - λ_i) / λ_i  =  Δλ / λ  =  1/R   is constant ∀i

    i.e. it has a constant spectral resolution given by:
    R = c / σ_kms

    The step of the grid is:
    log10(σ_kms / c + 1)
    =#
    oversampling = 1
    step = log10(σ_kms / 3e5 + 1) / oversampling

    # Interpolate input flux on the regular log-λ grid with proper resolution
    ll = 10. .^ range(log10(minimum(λ)), log10(maximum(λ)), step=step)
    # check with: extrema([(ll[i+1] - ll[i]) / ll[i] * 3e5 for i in 1:length(ll)-1]) .* oversampling, σ_kms
    ff = Spline1D(λ, flux, k=1)(ll)

    # Instrument response
    gauss(x; μ=0., σ=1.) = exp(-0.5 .* ((x .- μ) ./ σ).^2) ./ sqrt(2pi) ./ σ
    response = gauss.(-5:1. / oversampling:5)

    # Convolve model with instrument response
    i1 = div(length(response), 2) + 1
    i2 = i1 + length(ff) - 1
    out = conv(ff, response)[i1:i2] ./ oversampling

    # Interpolate back to original domain
    out = Spline1D(ll, out, k=1, bc="extrapolate")(λ)

    # Replace data close to the edges
    ee = (1 + 2σ_kms / 3e5) * minimum(λ)
    i = findall(λ  .< ee)
    out[i] = Spline1D(ll, ff, k=1, bc="extrapolate")(λ[i])

    ee = (1 - 2σ_kms / 3e5) * maximum(λ)
    i = findall(λ  .> ee)
    out[i] = Spline1D(ll, ff, k=1, bc="extrapolate")(λ[i])

    # @gp    λ flux "w lp t 'Model'"
    # @gp :- λ out  "w lp t 'Instrumental'"

    return out
end
