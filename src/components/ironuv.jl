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
    ia   = Spline1D(  ia[:,1],   ia[:,2], k=1, bc="zero")(λ)
    ib   = Spline1D(  ib[:,1],   ib[:,2], k=1, bc="zero")(λ)
    i34  = Spline1D( i34[:,1],  i34[:,2], k=1, bc="zero")(λ)
    i47  = Spline1D( i47[:,1],  i47[:,2], k=1, bc="zero")(λ)
    i191 = Spline1D(i191[:,1], i191[:,2], k=1, bc="zero")(λ)

    L = ib .+ i47 .+ i191
    L ./= sum(L .* (λ[2] - λ[1]))
    return (λ, L)
end


mutable struct ironuv <: AbstractComponent
    norm::Parameter
    fwhm::Float64
    L::Vector{Float64}  # pre-computed

    function ironuv(fwhm::Number)
        out = new(Parameter(1), float(fwhm), Vector{Float64}())
        out.norm.low = 0
        return out
    end
end


function prepare!(comp::ironuv, domain::Domain{1})
    @assert comp.fwhm > 900
    λ, L = ironuv_read()
    L = QSFit.conv_gauss(λ, L, sqrt(comp.fwhm^2. - 900^2) / 2.35)
    L ./= int_tabulated(λ, L)[1]
    comp.L = Spline1D(λ, L, k=1, bc="zero")(domain[:])
    return fill(NaN, length(domain))
end


function evaluate!(buffer, comp::ironuv, domain::Domain{1},
                   norm)
    buffer .= norm .* comp.L
end
