# ____________________________________________________________________
# Cutoff powerlaw
#
mutable struct cutoff_powerlaw <: AbstractComponent
    norm::Parameter
    x0::Parameter
    alpha::Parameter
    beta::Parameter

    function cutoff_powerlaw(x0::Number)
        out = new(Parameter(1),
                  Parameter(x0),
                  Parameter(-1),
                  Parameter(1))
        out.norm.low = 0
        out.x0.low = 0
        out.alpha.low = -5
        out.alpha.high = 5
        out.beta.low = 0.1
        out.beta.high = 10.
        
        return out
    end
end

compeval_cdata(comp::cutoff_powerlaw, domain::Domain_1D) = nothing
compeval_array(comp::cutoff_powerlaw, domain::Domain_1D) = fill(NaN, length(domain))

function evaluate(c::CompEval{cutoff_powerlaw, Domain_1D},
                   norm, x0, alpha, beta)
    x = c.domain[1]
    if alpha * beta < 0
        @warn "alpha and beta should have the same sign"
    end        
    c.buffer .= norm .* (x ./ x0).^alpha .* exp.(1 .- ((x ./ x0) .^ beta))
end
