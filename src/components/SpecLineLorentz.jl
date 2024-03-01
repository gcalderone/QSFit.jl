# ____________________________________________________________________
# SpecLineLorentz
#
mutable struct SpecLineLorentz <: AbstractSpecLineComp
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter

    function SpecLineLorentz(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0))
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.center.fixed = true
        return out
    end
end

function evaluate!(buffer::Vector{Float64}, comp::SpecLineLorentz, x::Domain{1},
                   norm, center, fwhm, voff)
    x0 = center - (voff / 3.e5) * center
    γ     = fwhm            / 2     / 3.e5 * center
    X = coords(x) .- x0
    buffer .= norm .* γ ./ (pi .* (X.^2. .+ γ^2.))
end
