# The following macro is taken from ReusePatterns.jl
macro copy_fields(T)
    out = Expr(:block)
    for name in fieldnames(__module__.eval(T))
        e = Expr(Symbol("::"))
        push!(e.args, name)
        push!(e.args, fieldtype(__module__.eval(T), name))
        push!(out.args, e)
    end
    return esc(out)
end


gauss(x, μ, σ) = exp.(-0.5 .* ((x .- μ) ./ σ).^2) ./ sqrt(2pi) ./ σ


function int_tabulated(x, y)
    @assert issorted(x)
    @assert all(isfinite.(x))
    @assert all(isfinite.(y))
    b = x[2:end] .- x[1:end-1]
    h = y[2:end] .+ y[1:end-1]
    return sum(b .* h) / 2
end


function planck(λ, T)
    h = 6.6260755  * 1e-27  # Planck's constant [erg s]
    c = 2.99792458 * 1e10   # Vacuum speed of light [cm / s]
    k = 1.380658   * 1e-16  # Boltzmann constant [erg / K]
    b = 2 * h * c^2. ./ (λ.^5.)
    d = h * c ./ (λ .* k .* T)
    return b ./ (exp.(d) .- 1)
end



function estimate_fwhm(λ, f; plot=false)
    imax = argmax(f)
    half_max = f[imax] / 2

    (i1, i2) = extrema(findall(f .>= half_max))

    j1 = maximum(findall(f[1:imax]   .<  half_max))
    j2 = minimum(findall(f[imax:end] .<  half_max)) + imax - 1
    @assert j1 < i1
    @assert i2 < j2
    x = collect(range(λ[j1], λ[i1], length=1000))
    y = Dierckx.Spline1D(λ, f)(x)
    λ1 = x[argmin(abs.(y .- half_max))]
    x = collect(range(λ[i2], λ[j2], length=1000))
    y = Dierckx.Spline1D(λ, f)(x)
    λ2 = x[argmin(abs.(y .- half_max))]
    fwhm = λ2 - λ1

    if plot
        yr = [extrema(f)...]
        @gp    :estfwhm λ f "w l notit"
        @gp :- :estfwhm λ1 .* [1,1] yr "w l notit lc rgb 'gray'"
        @gp :- :estfwhm λ2 .* [1,1] yr "w l notit lc rgb 'gray'"
        @gp :- :estfwhm [λ1, λ2] half_max .* [1,1]  "w l notit lc rgb 'gray'"
    end
    return fwhm
end


# Wavelength span is assumed to be equal to the initial FWHM plus the
# spectral resolution (in quadrature)
spectral_coverage(spec_λ, resolution, line::Union{SpecLineGauss, SpecLineLorentz, SpecLineVoigt}; kw...) =
    spectral_coverage(spec_λ, resolution, line.center.val, line.center.val * sqrt(line.fwhm.val^2 + resolution^2) / 3.e5; kw...)

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
