
qsfitversion() = v"0.1.0"
qsfit_data() = artifact"qsfit_data"

gauss(x, μ, σ) = exp.(-0.5 .* ((x .- μ) ./ σ).^2) ./ sqrt(2pi) ./ σ


function planck(λ, T)
    h = 6.6260755  * 1e-27  # Planck's constant [erg s]
    c = 2.99792458 * 1e10   # Vacuum speed of light [cm / s]
    k = 1.380658   * 1e-16  # Boltzmann constant [erg / K]
    b = 2 * h * c^2. ./ (λ.^5.)
    d = h * c ./ (λ .* k .* T)
    return b ./ (exp.(d) .- 1)
end


function convol(v, _k)
    k = reverse(_k)
    nk = length(k)
    @assert nk < length(v)
    @assert mod(nk, 2) == 1
    r = div(nk, 2)
    out = fill(0., length(v))
    for i in r+1:length(v)-r
        out[i-r:i+r] .+= v[i] .* k
    end
    return out ./ sum(abs.(k))
end


function interpol(y, x, newX; allow_extrapolations=false)
    out = fill(0., length(newX))  # TODO: Initialize with NaN, not zero...
    itp = interpolate((x,), y, Gridded(Linear()))
    ii = findall(minimum(x) .<= newX .<= maximum(x))
    out[ii] .= itp(newX[ii])

    if allow_extrapolations
        etp = extrapolate(itp, Interpolations.Line())
        ii = findall(newX .< minimum(x))
        out[ii] .= etp(newX[ii])
        ii = findall(maximum(x) .< newX)
        out[ii] .= etp(newX[ii])
    end
    return out
end
interpol(y, x, newX::Real; kw...) = interpol(y, x, [newX]; kw...)[1]


function estimate_fwhm(λ, f)
    i = argmax(f)
    maxv = f[i]
    i = findall(f .>= (maxv/2))
    (i1, i2) = extrema(i)
    fwhm = λ[i2] - λ[i1]
    return fwhm
end


function estimate_area(λ, f, from=nothing, to=nothing)
    isnothing(from)  &&  (from = λ[1])
    isnothing(to  )  &&  (to   = λ[end])
    itp = interpolate((λ,), f, Gridded(Linear()));
    return quadgk(itp, from, to)
end


function line_coverage(spec_λ, resolution, line_λ, fwhm)
    # Identify min/max wavelengths corresponding to expected FWHM
    λmin = line_λ .* (1 - fwhm / 3.e5 / 2.)
    λmax = line_λ .* (1 + fwhm / 3.e5 / 2.)

    # Calculate step corresponding to the spectral resolution
    δ = resolution / 3e5 * line_λ

    # How many intervals should be considered?
    intervals = Int(ceil((λmax - λmin) / δ))
    bins = range(λmin, λmax, length=intervals+1)

    # How many intervals are sampled?
    good = 0
    for i in 1:intervals
        good += sign(count(bins[i] .<= spec_λ .< bins[i+1]))
    end

    return (λmin, λmax, good / intervals)
end
