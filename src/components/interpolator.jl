struct Interpolator <: GModelFit.AbstractComponent
    orig_x::Vector{Float64}
    orig_y::Vector{Float64}
    interp_y::Vector{Float64}
    norm::GModelFit.Parameter

    function Interpolator(orig_x, orig_y)
        norm = GModelFit.Parameter(1)
        interp_y = Vector{Float64}()
        return new(orig_x, orig_y, interp_y, norm)
    end
end

function prepare!(comp::Interpolator, domain::AbstractDomain{1})
    itp = Dierckx.Spline1D(comp.orig_x, comp.orig_y, bc="error")
    append!(comp.interp_y, itp(coords(domain)))
end

result_length(comp::Interpolator, domain::AbstractDomain{1}) = length(comp.interp_y)

# Component evaluation (apply scaling factor)
function evaluate!(comp::Interpolator, domain::AbstractDomain{1}, output::Vector,
                   norm)
    output .= norm .* comp.interp_y
end
