import Gnuplot
import Gnuplot.recipe

function Gnuplot.recipe(res::QSFitResults)
    out = [Gnuplot.recipe(res.model)..., Gnuplot.recipe(res.pspec.data)]
    return reverse(out)
end
