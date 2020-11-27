# ____________________________________________________________________
# ironopt
#
function ironopt_read(file)
    a = readdlm(file, '\t', comments=true)
    df = DataFrame()
    df[!, :line] = string.(a[:,1])
    df[!, :transition] = string.(a[:,2])
    ii = findall(typeof.(a[:,3]) .== Float64)
    df[!, :ul] .= NaN; df[ii, :ul] = float.(a[ii,3])
    ii = findall(typeof.(a[:,4]) .== Float64)
    df[!, :wavelength] .= NaN; df[ii, :wavelength] = float.(a[ii,4])
    ii = findall(typeof.(a[:,5]) .== Float64)
    df[!, :aat] .= NaN; df[ii, :aat] = float.(a[ii,5])
    ii = findall(typeof.(a[:,6]) .== Float64)
    df[!, :wht] .= NaN; df[ii, :wht] = float.(a[ii,6])
            
    ii = findall(isfinite.(df[:, :wavelength])  .&  isfinite.(df[:, :wht]))
    df = df[ii,:]
    return df
end


mutable struct ironopt_cdata
    L::Vector{Float64}
    fwhm::Float64
end


mutable struct ironopt <: AbstractComponent
    norm::Parameter
    λ::Vector{Float64}
    A::Vector{Float64}
    fwhm::Float64
    
    function ironopt(file::String, fwhm::Number)
        df = ironopt_read(file)
        # Drop Balmer lines (they are accounted for in the main QSFit code)
        ii = findall(getindex.(df[:, :line], Ref(1:2)) .!= "H\$")
        df = df[ii,:]
        out = new(Parameter(1), df[:, :wavelength], df[:, :wht], float(fwhm))
        out.norm.low = 0
        return out
    end
end

function compeval_cdata(comp::ironopt, domain::Domain_1D)
    σ0 = comp.fwhm / 2.35 / 3.e5
    lmin, lmax = extrema(comp.λ)
    lmin -= 3 * σ0 * lmin
    lmax += 3 * σ0 * lmax
    λ = collect(lmin:(σ0*lmin/5):lmax)
    L = fill(0., length(λ))
    for ii in 1:length(comp.A)
        L .+= comp.A[ii] .* gauss(λ, comp.λ[ii], σ0 * comp.λ[ii])
    end
    L ./= sum(comp.A)
    out = interpol(L, λ, domain[1])
    return ironopt_cdata(out, comp.fwhm)
end
compeval_array(comp::ironopt, domain::Domain_1D) = fill(NaN, length(domain))

function evaluate(c::CompEval{ironopt, Domain_1D},
                  norm)
    c.buffer .= norm .* c.cdata.L
end

ironopt_broad( fwhm) = ironopt(qsfit_data() * "/VC2004/TabA1", fwhm)
ironopt_narrow(fwhm) = ironopt(qsfit_data() * "/VC2004/TabA2", fwhm)
