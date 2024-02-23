import Gnuplot
import Gnuplot.recipe

function Gnuplot.recipe(spec::T) where T <: AbstractSpectrum
    i = findall(spec.good)
    return Gnuplot.parseSpecs("set bars 0", title=spec.label * (isnothing(spec.z) ? "" : ", z=$(spec.z)"),
                              xlabel="[x" * string(spec.unit_x) * "]", ylabel="[x" * string(spec.unit_y) * "]",
                              spec.x   , spec.y   , spec.err   , "with yerr notit pt 0      lc rgb 'gray'",
                              spec.x[i], spec.y[i], spec.err[i], "with p    notit pt 1 ps 1 lc rgb 'black'")
end

function Gnuplot.recipe(spec::RestFrameSpectrum)
    out = @invoke Gnuplot.recipe(spec::AbstractSpectrum)
    append!(out, Gnuplot.parseSpecs(xlabel="Rest frame [x" * string(spec.unit_x) * "]"))
    return out
end


function Gnuplot.recipe(res::Results)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    i = findall(isnothing.(match.(r"SpecLine", ctypes)))
    keep = keys(res.bestfit)[i]
    out = [Gnuplot.recipe(res.spec)..., Gnuplot.recipe(res.bestfit, keep=keep)...]
    return out
end


function residuals(res::Results)
    resid = (values(res.pspec.data) .- res.bestfit()) ./ uncerts(res.pspec.data)
    return Gnuplot.parseSpecs(
        "set grid", "set key outside horizontal",
        coords(res.pspec.data.domain), resid, "with p t 'Residuals' lc rgb 'red'",
        Gnuplot.line(extrema(coords(res.pspec.data.domain)), 0., "with l notit dt 2 lc rgb 'black'"),
        coords(res.pspec.data.domain), cumsum(resid.^2), "with l t 'Cumulative residuals^2' lc rgb 'blue' axes x1y2")
end
