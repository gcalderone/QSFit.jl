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
    ii = findall(getindex.(df[:, :line], Ref(1:2)) .!= "H\$")
    df = df[ii,:]
    λ0 = df[:, :wavelength]
    A0 = df[:, :wht]

    σ0 = comp.fwhm / 2.35 / 3.e5
    lmin, lmax = extrema(λ0)
    lmin -= 3 * σ0 * lmin
    lmax += 3 * σ0 * lmax
    λ = collect(lmin:(σ0*lmin/5):lmax)
    L = fill(0., length(λ))
    for ii in 1:length(A0)
        L .+= A0[ii] .* gauss(λ, λ0[ii], σ0 * λ0[ii])
    end
    L ./= sum(A0)
    comp.L = Spline1D(λ, L, k=1, bc="zero")(domain[:])
    return fill(NaN, length(domain))
end


function evaluate!(buffer, comp::ironopt, domain::Domain{1},
                   norm)
    buffer .= norm .* comp.L
end

ironopt_broad( fwhm) = ironopt(qsfit_data() * "/VC2004/TabA1", fwhm)
ironopt_narrow(fwhm) = ironopt(qsfit_data() * "/VC2004/TabA2", fwhm)
