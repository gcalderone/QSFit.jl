# ____________________________________________________________________
# Smoothly Broken Power-law
#
mutable struct sbpl <: AbstractComponent
    norm::Parameter
    x0::Parameter
    alpha1::Parameter
    alpha2::Parameter
    delta::Parameter

    function sbpl(x0::Number)
        out = new(Parameter(1),
                  Parameter(x0),
                  Parameter(-1),
                  Parameter(-1),
                  Parameter(0.1))
        out.norm.low = 0
        out.x0.low = 0
        out.x0.fixed = true
        out.alpha1.low = -5; out.alpha1.high = 5
        out.alpha2.low = -5; out.alpha2.high = 5
        out.delta.low = 1e-3; out.delta.high = 1
        return out
    end
end

function evaluate!(buffer, comp::sbpl, x::Domain{1},
                   norm, x0, alpha1, alpha2, delta)
    xx = x ./ x0

    # The quantity `t = (x / x_b)^(1 / delta)` can become quite large.
    # To avoid overflow errors we will start by calculating its
    # natural logarithm:
    logt = log.(xx) ./ delta
    
    threshold = 30  # corresponding to exp(30) ~ 1e13
    i = findall(logt .> threshold)
    if length(i) > 0
        buffer[i] .= norm .* xx[i].^alpha2 .*
            (0.5 .^((alpha2 - alpha1) * delta))
    end

    i = findall(logt .< -threshold)
    if length(i) > 0
        buffer[i] .= norm .* xx[i].^alpha1 .*
            (0.5 .^((alpha2 - alpha1) * delta))
    end
    
    i = findall(abs.(logt) .< threshold)
    if length(i) > 0
        buffer[i] .= norm .* xx[i].^alpha1 .*
            ((0.5 .* (1 .+ xx[i].^(1/delta))) .^((alpha2 - alpha1) * delta))
    end
end


#=

x = Domain(0.7:0.01:500)

comp = QSFit.sbpl(20)
comp.alpha1.val =  2
comp.alpha2.val = -2
comp.delta.val = 0.5
ceval = GFit.CompEval(comp, x)

ceval.comp.delta.val = 0.5
evaluate!(ceval);
@gp x[:] ceval.buffer "w l" xlog=true ylog=true

ceval.comp.delta.val = 0.3
evaluate!(ceval);
@gp :- x[:] ceval.buffer "w l" xlog=true ylog=true

ceval.comp.delta.val = 0.1
evaluate!(ceval);
@gp :- x[:] ceval.buffer "w l" xlog=true ylog=true

=#
