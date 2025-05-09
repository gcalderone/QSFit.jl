# ____________________________________________________________________
# ironopt
#
function ironopt_read(file)
    a = readdlm(file, '\t', comments=true)
    df = OrderedDict()
    df[:line] = string.(a[:,1])
    df[:transition] = string.(a[:,2])
    df[:wavelength] = fill(NaN, length(df[:line]));
    df[:ul]  = fill(NaN, length(df[:line]));
    df[:aat] = fill(NaN, length(df[:line]));
    df[:wht] = fill(NaN, length(df[:line]));

    ii = findall(typeof.(a[ :,3]) .== Float64)
    df[:ul][ii] = float.(a[ii,3])

    ii = findall(typeof.(        a[ :,4]) .== Float64)
    df[:wavelength][ii] = float.(a[ii,4])

    ii = findall(typeof.( a[ :,5]) .== Float64)
    df[:aat][ii] = float.(a[ii,5])

    ii = findall(typeof.( a[ :,6]) .== Float64)
    df[:wht][ii] = float.(a[ii,6])

    ii = findall(isfinite.(df[:wavelength])  .&  isfinite.(df[:wht]))
    for (key, val) in df
        df[key] = val[ii]
    end

    return df
end


mutable struct ironopt <: AbstractComponent
    file::String
    fwhm::Float64
    L::Vector{Float64}
    norm::Parameter

    function ironopt(file::String, fwhm::Number)
        out = new(file, float(fwhm), Vector{Float64}(), Parameter(1))
        out.norm.low = 0
        return out
    end
end

function prepare!(comp::ironopt, domain::Domain{1})
    df = ironopt_read(comp.file)
    # Drop Balmer lines (they are accounted for in the main QSFit code)
    ii = findall(getindex.(df[:line], Ref(1:2)) .!= "H\$")
    for (key, val) in df
        df[key] = val[ii]
    end

    x = coords(domain)
    comp.L = x .* 0.
    for i in 1:length(df[:wavelength])
        λ = df[:wavelength][i]
        comp.L .+= df[:wht][i] .* gauss(x, λ, comp.fwhm / 2.355 / 3e5 * λ)
    end
    comp.L ./= int_tabulated(x, comp.L)
end


function evaluate!(ceval::CompEval{ironopt, Domain{1}},
                   norm)
    ceval.buffer .= norm .* ceval.comp.L
end

ironopt_broad( fwhm) = ironopt(joinpath(qsfit_data(), "VC2004", "TabA1"), fwhm)
ironopt_narrow(fwhm) = ironopt(joinpath(qsfit_data(), "VC2004", "TabA2"), fwhm)
