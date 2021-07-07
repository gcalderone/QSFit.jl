import Gnuplot
import Gnuplot.recipe

function Gnuplot.recipe(res::QSFitResults)
    out = [Gnuplot.recipe(res.pspec.data), Gnuplot.recipe(res.model)...]
    return reverse(out)
end
