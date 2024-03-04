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

function evaluate!(ceval::CompEval{powerlaw, Domain{1}},
                   norm, x0, alpha)
    ceval.buffer .= norm .* (coords(ceval.domain) ./ x0).^alpha
end
