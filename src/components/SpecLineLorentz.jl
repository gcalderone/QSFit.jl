# ____________________________________________________________________
# SpecLineLorentz
#
mutable struct SpecLineLorentz <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    index::Vector{Int}  # optimization
    norm_integrated::Bool

    function SpecLineLorentz(center::Number)
        out = new(Parameter(1),
                  Parameter(center),
                  Parameter(3000),
                  Parameter(0),
                  Vector{Int}(), true)
        @assert center > 0
        out.norm.low = 0
        out.center.low = 0
        out.fwhm.low = 0
        out.voff.low = 0
        out.center.fixed = true
        return out
    end
end

function prepare!(comp::SpecLineLorentz, domain::Domain{1})
    comp.index = collect(1:length(domain))
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::SpecLineLorentz, x::Domain{1},
                   norm, center, fwhm, voff)
    buffer[comp.index] .= 0.
    empty!(comp.index)

    x0 = center - (voff / 3.e5) * center
    hwhm = fwhm / 3.e5 * center / 2  # Note: this is in `center` units

    X = (x .- x0) ./ hwhm
    i = findall(abs.(X) .< 20) # optimization
    append!(comp.index, i)
    buffer[i] .= norm ./ (1 .+ X[i].^2.)
    if comp.norm_integrated
        buffer[i] ./= pi * hwhm
    end
end


#=
    x = Domain(500:1:1500.)
    comp = QSFit.SpecLineLorentz(1000.)
    comp.fwhm.val = 3e4
    ceval = GFit.CompEval(comp, x)
    GFit.evaluate_cached(ceval)
    @gp x ceval.buffer ./ maximum(ceval.buffer) "w l"
=#
