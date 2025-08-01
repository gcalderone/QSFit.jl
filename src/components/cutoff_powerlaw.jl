# ____________________________________________________________________
# Cutoff powerlaw
#
mutable struct cutoff_powerlaw <: AbstractComponent
    norm::Parameter
    x0::Parameter
    alpha::Parameter
    beta::Parameter # when beta is positive the cutoff affects x > x0 (and viceversa)

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


function evaluate!(::cutoff_powerlaw, domain::Domain{1}, output::Vector,
                   norm, x0, alpha, beta)
    x = coords(domain)
    output .= norm .* (x ./ x0).^alpha .* exp.(1 .- ((x ./ x0) .^beta))
end
