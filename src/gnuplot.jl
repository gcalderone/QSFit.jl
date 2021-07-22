import Gnuplot
import Gnuplot.recipe

export residuals

function Gnuplot.recipe(res::QSFitResults)
    out = [Gnuplot.recipe(res.model)..., Gnuplot.recipe(res.pspec.data)]
    return reverse(out)
end


function residuals(res::QSFitResults)
    out1 = Gnuplot.PlotElement(
        cmds=["set grid"],
        data = Gnuplot.DatasetBin(res.pspec.data.domain[:],
                                  (res.pspec.data.val .- res.model()) ./ res.pspec.data.val),
        plot="with p t 'Residuals' lc rgb 'red'")
    out2 = Gnuplot.PlotElement(
        data = Gnuplot.DatasetText([extrema(res.pspec.data.domain[:])...], [0., 0.]),
        plot="with l notit dt 2 lc rgb 'black'")
    return [out1, out2]
end
