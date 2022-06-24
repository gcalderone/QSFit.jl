import Gnuplot
import Gnuplot.recipe

export residuals

function Gnuplot.recipe(spec::Spectrum)
    i = findall(spec.good)
    return [Gnuplot.PlotElement(cmds=["set bars 0"], title=spec.label),
            Gnuplot.PlotElement(data=Gnuplot.DatasetBin(spec.λ[i], spec.flux[i], spec.err[i]), plot="with yerr notit ps 0 lc rgb 'black'"),
            Gnuplot.PlotElement(data=Gnuplot.DatasetBin(spec.λ   , spec.flux   , spec.err   ), plot="with yerr notit ps 0 lc rgb 'gray'")
            ]

end

function Gnuplot.recipe(res::JobResults)
    out = [Gnuplot.recipe(res.pspec.data), Gnuplot.recipe(res.model)...]
    return reverse(out)
end


function residuals(res::JobResults)
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
