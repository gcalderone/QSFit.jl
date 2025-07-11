#=
Notes on spectral resolution
A regular log-λ grid is characterized by a constant step:
logstep = log10(λ_i+1) - log10(λ_i)  =  log10(λ_i+1 / λ_i)

hence:
λ_i+1 / λ_i  =  const.

by subtracting a constant 1 to both sides
λ_i+1 / λ_i - 1  =  (λ_i+1 - λ_i) / λ_i  =  Δλ / λ  = const.

We interpret the constant as the reciprocal of the resolution:
Δλ / λ  =  1/R

i.e.
λ_i+1 / λ_i = 1/R + 1
logstep = log10(1/R + 1)

The resolution can also be expressed in km/s:
R = c / σ_kms
=#

# Calculate logstep corresponding to resolution R = λ / Δλ
logstep(R) = log10(1/R + 1)

function logregular_grid(x, R)
    ee = extrema(x)
    ee = (ee[1] - ee[1]/R, ee[2] + ee[2]/R)
    return 10. .^collect(log10(ee[1]):logstep(R):log10(ee[2]))
end



# ____________________________________________________________________
# GaussConv
#
# Convolution of a spectrum sampled on a log-regular grid with a Gaussian kernel
#
mutable struct GaussConv <: AbstractComponent
    depcname::Symbol
    R::Float64 # λ / Δλ
    kernel::Vector{Float64}
    buffer::Vector{Float64}
    depdomain::Domain{1}
    GaussConv(depcname, R) = new(depcname, R, Vector{Float64}(), Vector{Float64}(), Domain(1))
end

dependencies(comp::GaussConv) = [comp.depcname]
dependency_domain(comp::GaussConv, ::Domain) = comp.depdomain

function prepare!(comp::GaussConv, domain::Domain{1})
    comp.depdomain = Domain(logregular_grid(coords(domain), 2 * comp.R))  # sampling is twice the resolution
    comp.buffer = fill(0., length(comp.depdomain))
    nsigma = 5
    grid = -ceil(2 * nsigma):ceil(2 * nsigma)
    comp.kernel = gauss(grid, 0., 2.)
end

function evaluate!(comp::GaussConv, domain::Domain{1}, output::Vector, deps)
    direct_conv1d!(comp.buffer, deps[1], comp.kernel, Val(:edge_mirror))
    output .= Dierckx.Spline1D(coords(comp.depdomain), comp.buffer, k=1)(coords(domain))
end



#=
function mydelta(x)
    out = x .* 0.
    mm = div(length(out), 6)
    out[argmin(abs.(x .- 1250))] = 1.;
    out[argmin(abs.(x .- 1750))] = 1.;
    out[argmin(abs.(x .- 2600))] = 1.;
    return out
end

x = 10 .^(3:0.0001:3.5);

model = Model(:orig => @fd(x -> mydelta(x)), :convol => QSFit.GaussConv(:orig, 2000))
meval = GModelFit.ModelEval(model, Domain(x))
GModelFit.scan_model!(meval)
c = GModelFit.update_eval!(meval)

ox = coords(meval.cevals[:orig].domain)
oy = meval.cevals[:orig].tpar.buffer
QSFit.int_tabulated(ox, oy), QSFit.int_tabulated(x, c)  # <-- these should be equal
@gp xlog=true xr=extrema(x) ox oy "w lp" x c "w lp notit"

i = argmin(abs.(x .- 1250))
QSFit.estimate_fwhm(x[(i-20):(i+20)], c[(i-20):(i+20)]) / x[i] * 3e5 / 2.355  # should be ~150 km/s, i.e. 3e5 / 2000
=#
