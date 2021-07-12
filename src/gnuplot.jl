import Gnuplot
import Gnuplot.recipe

function Gnuplot.recipe(res::QSFitResults)
    out = [Gnuplot.recipe(res.model)..., Gnuplot.recipe((domain(res.model), res.pspec.data))]
    return reverse(out)
end
