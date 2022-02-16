# ____________________________________________________________________
# hostgalaxy
#
mutable struct hostgalaxy <: AbstractComponent
    norm::Parameter
    library::String
    template::String
    filename::String
    base::Vector{Float64}

    hostgalaxy(d::Dict) = hostgalaxy(d[:template], library=d[:library])

    function hostgalaxy(template::String; library="swire")
        if library == "swire"
            filename = template * "_template_norm.sed"
        end
        if library == "ILBERT2009"
            filename = template * ".sed"
        end
        filename = joinpath(qsfit_data(), library, filename)
        @assert isfile(filename) "File $filename do not exists"
        out = new(Parameter(1), library, template, filename, Vector{Float64}())
        out.norm.val = 1
        out.norm.low = 0
        return out
    end
end

function prepare!(comp::hostgalaxy, domain::Domain{1})
    # Read template
    d = readdlm(comp.filename, comments=true)
    @assert typeof(d) == Matrix{Float64}
    x = d[:,1]
    y = d[:,2]
    # Ensure wavelength values are unique
    jj = sortmerge(x, x)
    i = findall(countmatch(jj, 1) .== 1)
    x = x[i]
    y = y[i]
    # Ensure wavelength values are increasing monotonically
    i = sortperm(x)
    x = x[i]
    y = y[i]
    itp = Spline1D(x, y, k=1, bc="error")
    comp.base = itp(domain[:])
    comp.base ./= itp(5500.)
    return fill(NaN, length(domain))
end

function evaluate!(buffer, comp::hostgalaxy, domain::Domain{1},
                   norm)
    buffer .= norm .* comp.base
end
