import Gnuplot
import Gnuplot.recipe

function Gnuplot.recipe(spec::T) where T <: Spectrum
    i = findall(spec.good)
    xlabel = (isnan(spec.localtorestfactor)  ?  ""  :  "Rest frame ") * "wavelenth "
    return Gnuplot.parseSpecs("set bars 0", title=spec.label,
                              xlabel=xlabel * "[x" * string(spec.unit_x) * "]",
                              ylabel="[x" * string(spec.unit_y) * "]",
                              spec.x   , spec.y   , spec.err   , "with yerr notit pt 0      lc rgb 'gray'")
end


function Gnuplot.recipe(res::Results)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    i = findall(isnothing.(match.(r"SpecLine", ctypes)))
    keep = keys(res.bestfit)[i]
    out = [Gnuplot.recipe(res.spec)..., Gnuplot.recipe(res.bestfit, keep=keep)...]
    return out
end


function residuals(bestfit::GModelFit.ModelSnapshot, data::Measures)
    resid = (values(data) .- bestfit()) ./ uncerts(data)
    return Gnuplot.parseSpecs(
        "set grid", "set key outside horizontal",
        coords(data.domain), resid, "with p t 'Residuals' lc rgb 'red'",
        Gnuplot.line(extrema(coords(data.domain)), 0., "with l notit dt 2 lc rgb 'black'"),
        coords(data.domain), cumsum(resid.^2), "with l t 'Cumulative residuals^2' lc rgb 'blue' axes x1y2")
end
residuals(res::Results) = residuals(res.bestfit, res.pspec.data)
