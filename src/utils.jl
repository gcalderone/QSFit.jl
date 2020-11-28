using Pkg, Pkg.Artifacts
using Interpolations, QuadGK



qsfitversion() = v"0.1.0"
qsfit_data() = artifact"qsfit_data"

function gauss(x, μ, σ)
    return exp.(-0.5 .* ((x .- μ) ./ σ).^2) ./ sqrt(2pi) ./ σ
end

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

function interpol(y, x, X)
    out = fill(0., length(X))
    itp = interpolate((x,), y, Gridded(Linear()))
    ii = findall(minimum(x) .<= X .<= maximum(x))
    out[ii] .= itp(X[ii])
    return out
end
interpol(y, x, X::Number) = interpol(y, x, [X])[1]

function interpol1(y, x, X)
    out = fill(0., length(X))
    itp = interpolate((x,), y, Gridded(Linear()))
    ii = findall(minimum(x) .<= X .<= maximum(x))
    out[ii] .= itp(X[ii])
    etp = extrapolate(itp, Interpolations.Line())
    ii = findall(minimum(x) .> X)
    out[ii] .= etp(X[ii])
    ii = findall(maximum(x) .< X)
    out[ii] .= etp(X[ii])
    return out
end
interpol1(y, x, X::Number) = interpol1(y, x, [X])[1]


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
