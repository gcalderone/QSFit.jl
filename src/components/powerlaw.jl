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

ceval_data(domain::Domain_1D, comp::powerlaw) = (nothing, length(domain))

function evaluate(c::CompEval{Domain_1D, powerlaw},
                   norm, x0, alpha)
    x = c.domain[1]
    c.eval .= norm .* (x ./ x0).^alpha
end
