gauss(x, μ, σ) = exp.(-0.5 .* ((x .- μ) ./ σ).^2) ./ sqrt(2pi) ./ σ

# Calculate sampling resolutions of a wavelength vector
function sampling_resolutions(x::AbstractVector{<: Real})
    if !issorted(x)
        return sampling_resolution(sort(x))
    end
    return 1 ./ (x[2:end] ./ x[1:end-1] .- 1)
end

# Calculate number of logsteps corresponding to a spectral resolution Rspec on a
# logregular grid with sampling resolution Rsampling.
nlogsteps(Rsampling, Rspec) = logstep(Rspec)  /  logstep(Rsampling)

# Generate a gaussian kernel with spectral resolution Rspec on a logregular grid
# with resolution Rspec, extending in the range +/-nsigma.
function gauss_kernel(Rsampling, Rspec, nsigma=5)
    N = nlogsteps(Rsampling, Rspec)
    grid = -ceil(N * nsigma):ceil(N * nsigma)
    return gauss(grid, 0., N)
end

# Interpolate on a logregular grid with constant resolution R
function interpolate_on_logregular_grid(x, y, R)
    out_x = 10 .^ range(log10(minimum(x)), log10(maximum(x)), step=logstep(R))
    out_y = Dierckx.Spline1D(x, y, k=1)(out_x)
    return (out_x, out_y)
end


function gauss_broadening(x, y, σ_kms)
    Rsampling = 2 * maximum(sampling_resolutions(x))
    grid_x, grid_y = interpolate_on_logregular_grid(x, y, Rsampling)
    kernel = gauss_kernel(Rsampling, 3e5/σ_kms)
    @assert isodd(length(kernel))
    convolved = direct_conv1d(grid_y, kernel)
    # Interpolate back to original domain
    return Dierckx.Spline1D(grid_x, convolved, k=1, bc="extrapolate")(x)
end


# This provides results similar to https://docs.sciml.ai/Integrals/stable/basics/SampledIntegralProblem/, while being significantly simpler.
function int_tabulated(x, y)
    @assert all(isfinite.(x))
    @assert all(isfinite.(y))
    if !issorted(x)
        i = sortperm(x)
        return int_tabulated(x[i], y[i])
    end
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


# Enlarge a wavelength range to account for the instrumental resolution
function broaden_range(λmin, λmax, resolution)
    λcenter = (λmin + λmax) / 2
    λspan   =  λmax - λmin
    λspan   = sqrt(λspan^2 + (λcenter / resolution)^2)
    λmin = λcenter - λspan / 2
    λmax = λcenter + λspan / 2
    return (λmin, λmax)
end


# Wavelength range spanned by components
wavelength_span(line::AbstractSpecLineComp) = (line.center.val * (1 + line.fwhm.val / 3.e5 / 2),
                                               line.center.val * (1 - line.fwhm.val / 3.e5 / 2))

wavelength_span(comp::ironuv)  = extrema(ironuv_read()[1])
wavelength_span(comp::ironopt) = extrema(ironopt_read(comp.file)[:wavelength])


function spectral_coverage(spec_λ::Vector{Float64}, resolution::Float64,
                           comp::AbstractComponent; comp_Npoints::Int=5)
    # Get wavelength range spanned by the component and enlarge it to
    # account for the instrumental resolution
    comp_λmin, comp_λmax = wavelength_span(comp)
    comp_λmin, comp_λmax = broaden_range(comp_λmin, comp_λmax, resolution)

    # # Calculate step corresponding to the spectral resolution
    # Δλ = resolution / 3e5 * (comp_λmin + comp_λmax) / 2
    # res_Npoints = Int(ceil((comp_λmax - comp_λmin) / Δλ))

    # Calculate Npoints within the component range at given instrument resolution
    #=
    ((l1 + l0) / 2) / (l1 - l0) = resolution
    ((l1 + l0) / 2) = (l1 - l0) * resolution
    l0/2 + l0 * resolution = l1 * resolution - l1/2
    l0 * (1/2 + resolution) = l1 * (resolution - 1/2)
    l1 = l0 * (1/2 + resolution) / (resolution - 1/2)
    =#
    res_Npoints = 0
    l = comp_λmin
    while l <= comp_λmax
        l *= (1/2 + resolution) / (resolution - 1/2)
        res_Npoints += 1
    end

    # Consider Npoints as the largest among res_Npoints and com_Npoints
    Npoints = max(res_Npoints, comp_Npoints)

    # Calculate the grid of potentially sampled intervals
    grid = 10. .^range(log10(comp_λmin), log10(comp_λmax), length=Npoints+1)

    # How many intervals are actually sampled in the spectrum?
    good = 0
    for i in 1:Npoints
        good += sign(count(grid[i] .<= spec_λ .< grid[i+1]))
    end

    return (comp_λmin, comp_λmax, good / Npoints)
end


macro batch_when_threaded(args...)
    if Threads.nthreads() > 1
        out = Expr(:macrocall, Symbol("@batch"), LineNumberNode(1, nothing))
        append!(out.args, args)
        return esc(out)
    else
        for arg in args
            if isa(arg, Expr)
                if arg.head == :for
                    return esc(arg)
                end
            end
        end
    end
end
