# ____________________________________________________________________
# hostgalaxy
#
mutable struct hostgalaxy <: AbstractComponent
    norm::Parameter
    template::String
    base::Vector{Float64}

    function hostgalaxy(template::String)
        out = new(Parameter(1), template, Vector{Float64}())
        out.norm.val = 1
        out.norm.low = 0
        return out
    end
end

function prepare!(comp::hostgalaxy, domain::Domain{1})
    d = readdlm(qsfit_data() * "/swire/" * comp.template * "_template_norm.sed")
    @assert typeof(d) == Matrix{Float64}
    itp = interpolate((d[:,1],), d[:,2], Gridded(Linear()))
    comp.base = collect(itp(domain[:]))
    comp.base ./= itp(5500.)
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::hostgalaxy, domain::Domain{1},
                   norm)
    buffer .= norm .* comp.base
end
