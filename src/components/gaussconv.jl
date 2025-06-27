# ____________________________________________________________________
# GaussConv
#
# Convolution of a spectrum sampled on a log-regular grid with a Gaussian kernel
#
mutable struct GaussConv <: AbstractComponent
    input::Symbol
    center::Float64
    R::Float64 # λ / Δλ
    kernel::Vector{Float64}
    GaussConv(input, center, R) = new(input, center, R, Vector{Float64}())
end

dependencies(comp::GaussConv) = [comp.input]

function prepare!(comp::GaussConv, domain::Domain{1})
    x = coords(domain)
    # Check sampling resolution around central point
    i = argmin(abs.(comp.center .- x))
    i1 = i - 10;  (i1 < 1)          &&  (i1 = 1)
    i2 = i + 10;  (i2 > length(x))  &&  (i2 = length(x))
    Rsampling = sampling_resolutions(x[i1:i2])
    # Ensure sampling resolution is approximately constant
    if !(0 <= ((maximum(Rsampling) - minimum(Rsampling)) / mean(Rsampling)) < 5.e-2)
        @warn "Grid may not be log-regular for component $(comp.input), hence broadening may not be reliable"
    end
    Rsampling = mean(Rsampling)
    (Rsampling < comp.R)  &&  (@warn  "Sampling resolution for $(comp.input) ($(Rsampling)) is too small w.r.t spectral resolution ($(comp.R)) for the required broadening")
    comp.kernel = gauss_kernel(Rsampling, comp.R)
end

function evaluate!(comp::GaussConv, domain::Domain{1}, output::Vector, deps)
    direct_conv1d!(output, deps[1], comp.kernel, Val(:edge_mirror))
end

#=
x = 10 .^(3:0.0001:3.5);
y = x .* 0.;
mm = div(length(y), 6);
y[ mm] = 1.;
y[3mm] = 1.;
y[5mm] = 1.;
comp = QSFit.GaussConv(:a, 2000)
ceval = GModelFit.CompEval(comp, Domain(x))
push!(ceval.tpar.deps, GModelFit.CompEvalT{Float64}(y))
c = GModelFit.update_eval!(ceval, Float64[])

QSFit.int_tabulated(x, y), QSFit.int_tabulated(x, c)  # <-- these should be equal
@gp xlog=true xr=extrema(x) x y "w l notit" x c "w l notit"
QSFit.estimate_fwhm(x[2mm:4mm], c[2mm:4mm]) / x[3mm] * 3e5 / 2.355  # should be ~150 km/s, i.e. 3e5 / 2000
=#
