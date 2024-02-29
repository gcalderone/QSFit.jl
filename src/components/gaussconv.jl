# ____________________________________________________________________
# GaussConv
#
# Convolution of a spectrum sampled on a log-regular grid with a Gaussian kernel
#
mutable struct GaussConv <: AbstractComponent
    input::Symbol
    resolution::Float64 # λ / Δλ
    n::Int
    kernel::Vector{Float64}
    GaussConv(input, resolution) = new(input, resolution, 5, Vector{Float64}())
end

dependencies(comp::GaussConv) = [comp.input]

function prepare!(comp::GaussConv, domain::Domain{1})
    x = coords(domain)
    # Check sampling resolution in input grid
    logsteps = log10.(x[2:end]) .- log10.(x[1:end-1])
    # Ensure resolution is approximately constant
    @assert 0 <= ((maximum(logsteps) - minimum(logsteps)) / mean(logsteps)) < 1.e-2 "Grid is not log-regular"
    logstep = mean(logsteps)
    # Calculate width of kernel in units of steps
    sigma = log10(1/comp.resolution + 1) / logstep
    @assert sigma > 1 "Sampling resolution is coarser than instrumental resolution"
    grid_kernel = -ceil(comp.n * sigma):ceil(comp.n * sigma)
    empty!(comp.kernel)
    append!(comp.kernel, gauss(grid_kernel, 0, sigma))
    return fill(NaN, length(domain))
end

function evaluate!(buffer::Vector{Float64}, comp::GaussConv, x::Domain{1}, arg)
    direct_conv1d!(buffer, arg[1], comp.kernel)
end

#=
x = 10 .^(3:0.0001:3.5);
y = x .* 0.;
mm = div(length(y), 6);
y[ mm] = 1.;
y[3mm] = 1.;
y[5mm] = 1.;
comp = QSFit.GaussConv(:a, 2000)
c = comp(Domain(x), y)
QSFit.int_tabulated(x, y), QSFit.int_tabulated(x, c)  # <-- these should be equal
@gp xlog=true xr=extrema(x) x y "w l notit" x c "w l notit"
QSFit.estimate_fwhm(x[2mm:4mm], c[2mm:4mm]) / x[3mm] * 3e5 / 2.355  # should be ~150 km/s, i.e. 3e5 / 2000
=#
