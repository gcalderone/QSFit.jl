# ____________________________________________________________________
# hostgalaxy
#
mutable struct hostgalaxy <: AbstractComponent
    norm::Parameter
    template::String
    base::Vector{Float64}

    function hostgalaxy(template::String)
        filename = qsfit_data() * "/swire/" * template * "_template_norm.sed"
        if !isfile(filename)
            filename = template
        end
        @assert isfile(filename) "File $filename do not exists"
        out = new(Parameter(1), filename, Vector{Float64}())
        out.norm.val = 1
        out.norm.low = 0
        return out
    end
end

function prepare!(comp::hostgalaxy, domain::Domain{1})
    d = readdlm(comp.template, comments=true)
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
