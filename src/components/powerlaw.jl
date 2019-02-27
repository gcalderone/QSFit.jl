# ____________________________________________________________________
# powerlaw
#
mutable struct powerlaw <: AbstractComponent
    norm::Parameter
    x0::Parameter
    alpha::Parameter

    function powerlaw(x0::Number)
        out = new(
            Parameter(1),
            Parameter(x0),
            Parameter(-1))
        out.norm.low = 0
        out.x0.low = 0
        out.x0.fixed = true
        out.alpha.low = -5; out.alpha.high = 5
        return out
    end
end

struct powerlaw_cdata <: AbstractComponentData; end
cdata(comp::powerlaw, domain::AbstractDomain) = powerlaw_cdata()

function evaluate!(cdata::powerlaw_cdata, output::Vector{Float64}, domain::Domain_1D,
                   norm, x0, alpha)
    x = domain[1]
    output .= norm .* (x ./ x0).^alpha
    return output
end
