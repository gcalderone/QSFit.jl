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

mutable struct hostgalaxy_cdata <: AbstractComponentData
    template::String
    base::Vector{Float64}
end

function cdata(comp::hostgalaxy, domain::AbstractDomain)
    d = readdlm(qsfitpath() * "/components/swire/" * comp.template * "_template_norm.sed")
    @assert typeof(d) == Matrix{Float64}
    itp = interpolate((d[:,1],), d[:,2], Gridded(Linear()))
    base = collect(itp(domain[1]))
    base ./= itp(5500.)
    return hostgalaxy_cdata(comp.template, base)
end

function evaluate!(cdata::hostgalaxy_cdata, output::Vector{Float64}, domain::Domain_1D,
                   norm)
    output .= norm .* cdata.base
    return output
end
