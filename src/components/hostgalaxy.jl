# ____________________________________________________________________
# hostgalaxy
#
mutable struct hostgalaxy <: AbstractComponent
    norm::Parameter
    template::String
    function hostgalaxy(template::String)
        out = new(Parameter(1), template)
        out.norm.val = 1
        out.norm.low = 0
        return out
    end
end

mutable struct hostgalaxy_cdata
    template::String
    base::Vector{Float64}
end

function ceval_data(domain::Domain_1D, comp::hostgalaxy)
    d = readdlm(qsfitpath() * "/data/swire/" * comp.template * "_template_norm.sed")
    @assert typeof(d) == Matrix{Float64}
    itp = interpolate((d[:,1],), d[:,2], Gridded(Linear()))
    base = collect(itp(domain[1]))
    base ./= itp(5500.)
    return (hostgalaxy_cdata(comp.template, base), length(domain))
end

function evaluate(c::CompEval{Domain_1D, hostgalaxy},
                  norm)
    c.eval .= norm .* c.cdata.base
end
