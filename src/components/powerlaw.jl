# ____________________________________________________________________
# powerlaw
#
mutable struct powerlaw <: AbstractComponent
    norm::Parameter
    x0::Parameter
    alpha::Parameter

    function powerlaw(x0::Number)
        out = new(Parameter(1),
                  Parameter(x0),
                  Parameter(-1))
        out.norm.low = 0
        out.x0.low = 0
        out.x0.fixed = true
        out.alpha.low = -5; out.alpha.high = 5
        return out
    end
end

compeval_cdata(comp::powerlaw, domain::Domain_1D) = nothing
compeval_array(comp::powerlaw, domain::Domain_1D) = fill(NaN, length(domain))

function evaluate(c::CompEval{powerlaw, Domain_1D},
                   norm, x0, alpha)
    x = c.domain[1]
    c.buffer .= norm .* (x ./ x0).^alpha
end
