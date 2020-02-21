function plot(qsfit, model)
    @gp "set bars 0" "set grid" "set key horizontal" xr=extrema(domain(model)) :-
    @gp :- title=qsfit.name * ", z=" * string(qsfit.z) * ", E(B-V)=" * string(qsfit.ebv) :-
    @gp :- xlabel="Rest frame wavelength [A]" ylabel="Lum. density [10^{42} erg s^{-1} A^{-1}]" :-
    @gp :- qsfit.domain[1][1] qsfit.data[1].val qsfit.data[1].unc "w yerr t 'Data' pt 0 lc rgb 'black'" :-
    @gp :- domain(model) model() "w l t 'Model' lw 2 lc rgb 'orange'" :-
    @gp :- domain(model) model(:continuum) "w l t 'Cont' dt 2 lc rgb 'red'" :-
    qsfit.options.use_galaxy  &&  (@gp :- domain(model) model(:continuum) .+ model(:galaxy) "w l t 'Cont+Host' lc rgb 'red'" :-)
    qsfit.options.use_balmer  &&  (@gp :- domain(model) model(:balmer) "w l t 'Balmer' dt 4 lc rgb 'dark-green'" :-)
    qsfit.options.use_ironuv  &&  (@gp :- domain(model) model(:ironuv) "w l t 'IronUV' lc rgb 'dark-green'" :-)
    qsfit.options.use_ironopt &&  (@gp :- domain(model) model(:ironoptbr) "w l t 'IronOptB' lc rgb 'dark-green'" :-)
    qsfit.options.use_ironopt &&  (@gp :- domain(model) model(:ironoptna) "w l t 'IronOptN' lc rgb 'dark-green'" :-)
    @gp :- domain(model) model(:broad_lines) "w l t 'Broad' lw 2 lc rgb 'blue'" :-
    @gp :- domain(model) model(:narrow_lines) "w l t 'Narrow' lw 2 lc rgb 'brown'" :-
    @gp :- domain(model) model(:unknown) "w l t 'Unknown' lc rgb 'purple'"


    @gp    :resid "set bars 0" "set grid" xr=extrema(domain(model)) :-
    @gp :- :resid title=qsfit.name * ", z=" * string(qsfit.z) * ", E(B-V)=" * string(qsfit.ebv) :-
    @gp :- :resid xlabel="Rest frame wavelength [A]" ylabel="Residuals [{/Symbol s}]" :-
    resid = (qsfit.data[1].val .- model()) ./ qsfit.data[1].unc
    nn = length(resid)
    @gp :- :resid qsfit.domain[1][1] resid fill(1., nn) "w p t 'Data' lc rgb 'black'" :-
    @gp :- :resid [extrema(qsfit.domain[1][1])...] [0,0] "w line notitle dt 2 lw 2 lt rgb 'orange'" :-

    params = collect(values(DataFitting.getparams(DataFitting.wrappee(model))))
    ifree = findall(.! getfield.(params, :cfixed))
    dof = nn - length(ifree)
    fs = cumsum(resid.^2) / dof
    @gp :- :resid "set y2label 'Cumulative {/Symbol c}^2_{red}'" :-
    @gp :- :resid "set ytics nomirror" :-
    @gp :- :resid "set y2tics" :-
    @gp :- :resid "set format y2 '%.1f'" :-
    @gp :- :resid "set y2range [" * string(minimum(fs)) * ":" * string(maximum(fs)) * "]" :-
    @gp :- :resid "set key bottom right" :-
    @gp :- :resid qsfit.domain[1][1] fs "w l title 'Reduced cumul. {/Symbol c}^2' ls 1 lw 2 lt rgb 'red' axes x1y2"
end
