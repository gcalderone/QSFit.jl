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

function compeval_cdata(comp::hostgalaxy, domain::Domain_1D)
    d = readdlm(qsfit_data() * "/swire/" * comp.template * "_template_norm.sed")
    @assert typeof(d) == Matrix{Float64}
    itp = interpolate((d[:,1],), d[:,2], Gridded(Linear()))
    base = collect(itp(domain[1]))
    base ./= itp(5500.)
    return hostgalaxy_cdata(comp.template, base)
end
compeval_array(comp::hostgalaxy, domain::Domain_1D) = fill(NaN, length(domain))

function evaluate(c::CompEval{hostgalaxy, Domain_1D},
                  norm)
    c.buffer .= norm .* c.cdata.base
end
