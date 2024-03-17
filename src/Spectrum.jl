export Spectrum

abstract type AbstractSpectrum end

mutable struct Spectrum <: AbstractSpectrum
    label::String
    z::Union{Nothing, Float64}
    extlaw::DustExtinction.ExtinctionLaw
    Av::Float64
    unit_x::Unitful.Quantity
    unit_y::Unitful.Quantity
    x::Vector{Float64}
    y::Vector{Float64}
    err::Vector{Float64}
    good::Vector{Bool}
    dered::Vector{Float64}
    resolution::Float64
    meta::Dict{Symbol, Any}
    function Spectrum(x::Vector{T}, y::Vector{T}, err::Vector{T};
                      good::Union{Nothing, Vector{Bool}}=nothing,
                      unit_x=u"angstrom",
                      unit_y=u"erg" / u"s" / u"cm"^2 / u"angstrom",
                      label="", z=nothing, resolution=NaN,
                      extlaw=OD94(Rv=3.1), Av=0.) where {T <: Real}
        @assert dimension(unit_x) == dimension(u"cm")
        @assert dimension(unit_y) == dimension(u"erg" / u"s" / u"cm"^2 / u"angstrom")
        isnothing(good)  &&  (good = fill(true, length(x)))
        if length(err) == 0
            @warn "Uncertainties were not provided: assuming 0.05 * (median(flux) + abs(flux))"
            err = 0.05 .* (median(y) .+ abs.(y))
        end
        @assert length(x) == length(y) == length(err) == length(good)
        @assert minimum(err) > 0 "Uncertainties must be positive!"

        if Av != 0.
            println("De-reddening factors @ 1450, 3000, 5100 AA: ",
                    deredden.(extlaw, [1450, 3000, 5100.], [1., 1., 1.], Av=Av))
            dered = deredden.(extlaw, x .* unit_x,  fill(1., length(x)), Av=Av)
        else
            dered = fill(1., length(x))
        end

        # Ensure wavelengths are sorted
        ii = sortperm(x)

        # Estimate sampling resolution
        l = x[ii]
        Rsampling = median(sampling_resolutions(l))
        if isnan(resolution)
            resolution = Rsampling / 2
            @warn "Resolution is not provided, assuming it is equal to half the sampling resolution: $resolution"
        end
        @assert Rsampling > resolution "Can't have sampling resolution ($Rsampling) greater than actual resolution ($resolution)"

        return new(label, z, extlaw, Av,
                   1. * unit_x, 1. * unit_y,
                   x[ii], y[ii], err[ii], good[ii], dered,
                   resolution, Dict{Symbol, Any}(:sampling_resolution => Rsampling))
    end
end


function convert_units!(spec::AbstractSpectrum, newunit_x::Quantity, newunit_y::Quantity)
    if  ((dimension(spec.unit_x) == dimension(u"cm"))  &&  (dimension(newunit_x) == dimension(u"Hz")))  ||
        ((dimension(spec.unit_x) == dimension(u"Hz"))  &&  (dimension(newunit_x) == dimension(u"cm")))
        c = 3e5 * u"km" / u"s"
        spec.y .*=      spec.x;  spec.err .*= spec.x;  spec.unit_y *=     spec.unit_x
        spec.x  .= 1 ./ spec.x                      ;  spec.unit_x  = c / spec.unit_x
        spec.y ./=      spec.x;  spec.err ./= spec.x;  spec.unit_y /=     spec.unit_x
    end
    spec.x   = ustrip.(uconvert.(unit(newunit_x), spec.x .* spec.unit_x)) ./ ustrip(newunit_x)
    spec.unit_x = newunit_x
    spec.y   = ustrip.(uconvert.(unit(newunit_y), spec.y   .* spec.unit_y)) ./ ustrip(newunit_y)
    spec.err = ustrip.(uconvert.(unit(newunit_y), spec.err .* spec.unit_y)) ./ ustrip(newunit_y)
    spec.unit_y = newunit_y
    return spec
end


function round_unit_scales!(spec::AbstractSpectrum)
    convert_units!(spec,
                   10. ^round(log10(ustrip(spec.unit_x))) * unit(spec.unit_x),
                   10. ^round(log10(ustrip(spec.unit_y))) * unit(spec.unit_y))
end


#=
spec = Spectrum(Val(:SDSS_DR10), "/home/gcalderone/my/work/software/qsfit/data/spec-0752-52251-0323.fits",
                           label="My SDSS source", z=0.3806);
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]

QSFit.convert_units!(spec, 1 * u"cm", 1 * u"erg" / u"s" / u"cm"^2 / u"angstrom");
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]

QSFit.convert_units!(spec, 1 * u"angstrom", 1e-17 * u"erg" / u"s" / u"cm"^2 / u"angstrom");
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]

QSFit.convert_units!(spec, 1 * u"Hz", 1 * u"erg" / u"s" / u"cm"^2 / u"Hz");
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]

QSFit.convert_units!(spec, 1 * u"Hz", 1 * u"Jy");
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]



spec = Spectrum(Val(:SDSS_DR10), "/home/gcalderone/my/work/software/qsfit/data/spec-0752-52251-0323.fits",
                           label="My SDSS source", z=0.3806);
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]

spec = QSFit.RestFrameSpectrum(spec, cosmology(h=0.70, OmegaM=0.3));

QSFit.convert_units!(spec, 1 *u"Hz", 1 * u"erg" / u"s" / u"Hz");
spec.unit_x, spec.x[1], spec.unit_y, spec.y[1]

=#


function show(io::IO, spec::AbstractSpectrum)
    println(io, @sprintf "%s: %s" string(typeof(spec)) spec.label)
    isnothing(spec.z)    ||  println(io, @sprintf "       z: %8.3f" spec.z)
    isnothing(spec.Av)   ||  println(io, @sprintf "      Av: %8.3f  (%s)" spec.Av join(split(string(spec.extlaw), "\n"), ","))
    println(io, @sprintf "  resol.: %8.3f" spec.resolution)
    println(io, @sprintf "    good: %8.3f%%, units: X=%s, Y=%s" count(spec.good) / length(spec.good) * 100 string(spec.unit_x) string(spec.unit_y))
end


# Resolution is ~2000
# (https://www.sdss.org/dr12/spectro/spectro_basics/.) Analysis of an
# HII region (e.g. spec-1176-52791-0591) yield 180 km/s FWHM for the
# [OIII]5007)
function Spectrum(::Val{:SDSS_DR10}, file::AbstractString; ndrop=100, resolution=2000., kws...)
    f = FITS(file)
    wl = 10 .^read(f[2], "loglam")
    flux = float.(read(f[2], "flux"))
    ivar = float.(read(f[2], "ivar"))
    mask = read(f[2], "and_mask")
    close(f)

    good = convert(Vector{Bool}, ((mask .== 0)  .&
                                  (ivar .> 0)   .&
                                  (flux .> 0)))
    if ndrop > 0
        good[1:ndrop] .= false
        good[end-ndrop+1:end] .= false
    end

    out = Spectrum(wl, flux , sqrt.(1 ./ ivar);
                   unit_x = u"angstrom",
                   unit_y = 1e-17 * u"erg" / u"s" / u"cm"^2 / u"angstrom",
                   good=good, resolution=resolution, kws...)
    return out
end


function Spectrum(::Val{:ASCII}, file::AbstractString;
                  columns=[1,2,3],
                  kws...)
    @assert length(columns) >= 2

    wl    = Vector{Float64}()
    flux = Vector{Float64}()
    unc  = Vector{Float64}()
    for l in readlines(file)
        l = strip(strip(l))
        (l[1] == '#')  &&  continue
        s = string.(split(l, keepempty=false))
        @assert length(s) >= maximum(columns)
        push!(wl  , Meta.parse(s[columns[1]]))
        push!(flux, Meta.parse(s[columns[2]]))
        if length(columns) == 3
            push!(unc, Meta.parse(s[columns[3]]))
        end
    end
    out = Spectrum(wl, flux , unc; kws...)
    return out
end


# ====================================================================
mutable struct RestFrameSpectrum <: AbstractSpectrum
    label::String
    z::Float64
    extlaw::DustExtinction.ExtinctionLaw
    Av::Float64
    unit_x::Unitful.Quantity
    unit_y::Unitful.Quantity
    x::Vector{Float64}
    y::Vector{Float64}
    err::Vector{Float64}
    good::Vector{Bool}
    dered::Vector{Float64}
    resolution::Float64
    meta::Dict{Symbol, Any}
end

function RestFrameSpectrum(spec::Spectrum, cosmology::Cosmology.AbstractCosmology)
    @assert !isnothing(spec.z)  &&  (spec.z > 0)
    @assert dimension(spec.unit_x) == dimension(u"cm")
    @assert dimension(spec.unit_y) == dimension(u"erg" / u"s" / u"cm"^2 / u"angstrom")
    ld = uconvert(u"cm", luminosity_dist(cosmology, spec.z))
    flux2lum = 4pi * ld^2
    out = RestFrameSpectrum(spec.label, spec.z, spec.extlaw, spec.Av,
                            spec.unit_x, spec.unit_y * flux2lum,
                            spec.x    ./ (1 + spec.z),
                            spec.y    .* (1 + spec.z),
                            spec.err  .* (1 + spec.z),
                            spec.good, spec.dered, spec.resolution, spec.meta)
    return out
end
