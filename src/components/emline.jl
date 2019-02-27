# ____________________________________________________________________
# emline
#
mutable struct emline <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter
    voff::Parameter
    profile::Symbol
    function emline(center::Number)
        out = new(
            Parameter(1),
            Parameter(center),
            Parameter(3000),
            Parameter(0),
            :Gaussian)
        out.norm.low = 0
        out.center.low = 0        
        out.fwhm.low = 0        
        out.voff.low = 0        
        out.center.fixed = true
        return out
    end
end

mutable struct emline_cdata <: AbstractComponentData
    profile::Symbol
end
cdata(comp::emline, domain::AbstractDomain) = emline_cdata(comp.profile)

maxvalue(w::DataFitting.UI{DataFitting.WComponent}) = maxvalue(DataFitting.wrappee(w).comp)
function maxvalue(comp::emline)
    d = Domain([comp.center.val])
    cd = cdata(comp, d)
    out = [0.]
    evaluate!(cd, out, d, comp.norm.val, comp.center.val, comp.fwhm.val, comp.voff.val)
    return out[1]
end

function evaluate!(cdata::emline_cdata, output::Vector{Float64}, domain::Domain_1D,
                   norm, center, fwhm, voff)

    x = domain[1]

    if cdata.profile == :Lorentzian
        x0 = center - (voff / 3.e5) * center
        xx = (x .- x0) ./ (fwhm / 3.e5 * center / 2.)
        output .= norm ./ (1 .+ xx.^2.)
        return output
    end

    @assert cdata.profile == :Gaussian
    x0 = center - (voff / 3.e5) * center
    sigma = (fwhm  / 3.e5) * center / 2.35
    #improve performance (SQRT(10*2.) = 4.4721360)
    i = findall(abs.(x .- x0) .< 4.4721360 * sigma)
    ee = ((x[i] .- x0) ./ sigma).^2. ./ 2.
    line = norm * exp.(-ee) ./ 2.50663 ./ sigma #SQRT(2*!PI) = 2.50663
    output .= 0.
    output[i] .= line
    return output
end
