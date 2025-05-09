# ____________________________________________________________________
# hostgalaxy
#
function list_hostgalaxy_templates()
    out = Vector{Dict}()
    for file in readdir(joinpath(qsfit_data(), "swire"))
        m = match(r"(.*)_template_norm.sed", file)
        isnothing(m)  &&  continue
        d = Dict(:library=>"swire", :template=>string(m.captures[1]))
        push!(out, d)
    end
    for file in readdir(joinpath(qsfit_data(), "ILBERT2009"))
        m = match(r"(.*).sed$", file)
        isnothing(m)  &&  continue
        d = Dict(:library=>"ILBERT2009", :template=>string(m.captures[1]))
        push!(out, d)
    end
    return out
end


mutable struct hostgalaxy <: AbstractComponent
    norm::Parameter
    refwl::Float64
    library::String
    template_name::String
    filename::String
    template::Vector{Float64}

    function hostgalaxy(template_name::String; library="swire", refwl=5500.)
        if library == "swire"
            filename = template_name * "_template_norm.sed"
        end
        if library == "ILBERT2009"
            filename = template_name * ".sed"
        end
        filename = joinpath(qsfit_data(), library, filename)
        @assert isfile(filename) "File $filename do not exists"
        out = new(Parameter(1), refwl, library, template_name, filename, Vector{Float64}())
        out.norm.val = 1
        out.norm.low = 0
        out.refwl = refwl
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
    itp = Dierckx.Spline1D(x, y, k=1, bc="error")
    comp.template = fill(0., length(domain))
    i = findall(minimum(x) .< coords(domain) .< maximum(x))
    comp.template[i] = itp(coords(domain)[i])
    comp.template ./= itp(comp.refwl)
end

function evaluate!(ceval::CompEval{hostgalaxy, Domain{1}},
                   norm)
    ceval.buffer .= norm .* ceval.comp.template
end
