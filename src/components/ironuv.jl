# ____________________________________________________________________
# ironuv
#
function ironuv_read()
    ia   = readdlm(qsfit_data() * "/VW2001/Fe_UVtemplt_A.asc")
    ib   = readdlm(qsfit_data() * "/VW2001/Fe_UVtemplt_B.asc")
    i34  = readdlm(qsfit_data() * "/VW2001/Fe3UV34modelB2.asc")
    i47  = readdlm(qsfit_data() * "/VW2001/Fe3_UV47.asc")
    i191 = readdlm(qsfit_data() * "/VW2001/Fe2_UV191.asc")

    lmin, lmax = extrema([ia[:,1]; ib[:,1]; i34[:,1]; i47[:,1]; i191[:,1]])

    λ = lmin:(ia[2,1]-ia[1,1]):lmax
    ia   = interpol(  ia[:,2],   ia[:,1], λ)
    ib   = interpol(  ib[:,2],   ib[:,1], λ)
    i34  = interpol( i34[:,2],  i34[:,1], λ)
    i47  = interpol( i47[:,2],  i47[:,1], λ)
    i191 = interpol(i191[:,2], i191[:,1], λ)

    L = ib .+ i47 .+ i191
    L ./= sum(L .* (λ[2] - λ[1]))
    return (λ, L)
end


mutable struct ironuv_cdata
    L::Vector{Float64}
    fwhm::Float64
end


mutable struct ironuv <: AbstractComponent
    norm::Parameter
    λ::Vector{Float64}
    L::Vector{Float64}
    fwhm::Float64

    function ironuv(fwhm::Number)
        λ, L = ironuv_read()
        out = new(Parameter(1), λ, L, float(fwhm))
        out.norm.low = 0
        return out
    end
end


function compeval_cdata(comp::ironuv, domain::Domain{1})
    @assert comp.fwhm > 900
    σ0 = comp.fwhm / 2.35 / 3.e5
    lmin, lmax = extrema(comp.λ)
    lmin -= 3 * σ0 * lmin
    lmax += 3 * σ0 * lmax
    lmin -= 100
    lmax += 100

    logλ = collect(log10(lmin):log10(lmax/(lmax-σ0*lmax)):log10(lmax))
    logL = interpol(comp.L, comp.λ, 10 .^logλ)
    σ = sqrt(comp.fwhm^2. - 900^2) / 2.35 / 3.e5
    kernel = gauss(logλ, mean(logλ), σ)

    conv = real.(ifft(fft(logL) .* fft(kernel)))
    conv = [conv[div(length(conv), 2):end]; conv[1:div(length(conv), 2)-1]]

    conv = interpol(conv, 10 .^logλ, domain[1])
    return ironuv_cdata(conv, comp.fwhm)
end


function evaluate(buffer, comp::ironuv, domain::Domain{1}, cdata,
                  norm)
    buffer .= norm .* cdata.L
end

