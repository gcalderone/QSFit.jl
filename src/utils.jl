
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


function estimate_fwhm(λ, f)
    i = argmax(f)
    maxv = f[i]
    i = findall(f .>= (maxv/2))
    (i1, i2) = extrema(i)
    fwhm = λ[i2] - λ[i1]
    return fwhm
end


function int_tabulated(λ, f; int_k=3)
    @assert issorted(λ)
    itp = Spline1D(λ, f, k=int_k, bc="error")
    return quadgk(itp, λ[1], λ[end])
end


spectral_coverage(spec_λ, resolution, line::Union{SpecLineGauss, SpecLineAsymmGauss, SpecLineLorentz}; kw...) =
    spectral_coverage(spec_λ, resolution, line.center.val, line.center.val * line.fwhm.val / 3.e5; kw...)


function spectral_coverage(spec_λ, resolution, comp::ironuv; kw...)
    rr = extrema(ironuv_read()[1])
    return spectral_coverage(spec_λ, resolution, mean(rr), rr[2]-rr[1]; kw...)
end

function spectral_coverage(spec_λ, resolution, comp::ironopt; kw...)
    rr = extrema(QSFit.ironopt_read(comp.file)[:wavelength])
    return spectral_coverage(spec_λ, resolution, mean(rr), rr[2]-rr[1]; kw...)
end


function spectral_coverage(spec_λ::Vector{Float64}, resolution::Float64,
                           center_λ::Float64, span_λ::Float64; min_steps::Int=5)

    # Identify min/max wavelengths
    λmin = center_λ - span_λ / 2.
    λmax = center_λ + span_λ / 2.

    # Calculate step corresponding to the spectral resolution
    δ = resolution / 3e5 * center_λ

    # How many steps should be considered?
    steps = Int(ceil((λmax - λmin) / δ))
    (steps < min_steps)  &&  (steps = min_steps)
    bin_edges = range(λmin, λmax, length=steps+1)

    # How many intervals are sampled?
    good = 0
    for i in 1:steps
        good += sign(count(bin_edges[i] .<= spec_λ .< bin_edges[i+1]))
    end

    return (λmin, λmax, good / steps)
end


function invert_dictionary(d::OrderedDict{K, V}) where {K, V}
    out = OrderedDict{V, Vector{K}}()
    kk = collect(keys(d))
    vv = collect(values(d))
    for v in unique(vv)
        out[v] = kk[findall(v .== vv)]
    end
    return out
end
