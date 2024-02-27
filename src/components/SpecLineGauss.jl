# ____________________________________________________________________
# SpecLineGauss
#
mutable struct SpecLineGauss <: AbstractSpecLineComp
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    resolution::Float64

    function SpecLineGauss(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  0.)
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.center.fixed = true
        return out
    end
end

function evaluate!(buffer::Vector{Float64}, comp::SpecLineGauss, x::Domain{1},
                   norm, center, fwhm, voff)
    x0 = center - (voff / 3.e5) * center
    σ_res = comp.resolution / 2.355 / 3.e5 * center
    σ     = fwhm            / 2.355 / 3.e5 * center
    σ = sqrt(σ^2 + σ_res^2)
    X = (coords(x) .- x0) ./ σ

    function profile(x)
        (abs(x) > 4)  &&  (return 0.)
        return norm * exp(-x^2 / 2) / sqrt(2pi) / σ
    end
    map!(profile, buffer, X)
end
