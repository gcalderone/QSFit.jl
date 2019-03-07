function plot(qsfit, model)
    @gp "set bars 0" "set grid" xr=extrema(model(:domain)) :-
    @gp :- title=qsfit.name * ", z=" * string(qsfit.z) * ", E(B-V)=" * string(qsfit.ebv) :-
    @gp :- xlabel="Rest frame wavelength [A]" ylabel="Lum. density [10^{42} erg s^{-1} A^{-1}]" :-
    @gp :- qsfit.domain[1][1] qsfit.data[1].val qsfit.data[1].unc "w yerr t 'Data' pt 0 lc rgb 'black'" :-
    @gp :- model(:domain) model() "w l t 'Model' lw 2 lc rgb 'orange'" :-
    @gp :- model(:domain) model(:continuum) "w l t 'Cont' dt 2 lc rgb 'red'" :-
    qsfit.options.use_galaxy  &&  (@gp :- model(:domain) model(:continuum) .+ model(:galaxy) "w l t 'Cont+Host' lc rgb 'red'" :-)
    qsfit.options.use_balmer  &&  (@gp :- model(:domain) model(:balmer) "w l t 'Balmer' dt 4 lc rgb 'dark-green'" :-)
    qsfit.options.use_ironuv  &&  (@gp :- model(:domain) model(:ironuv) "w l t 'IronUV' lc rgb 'dark-green'" :-)
    qsfit.options.use_ironopt &&  (@gp :- model(:domain) model(:ironoptbr) "w l t 'IronOptB' lc rgb 'dark-green'" :-)
    qsfit.options.use_ironopt &&  (@gp :- model(:domain) model(:ironoptna) "w l t 'IronOptN' lc rgb 'dark-green'" :-)
    @gp :- model(:domain) model(:broad_lines) "w l t 'Broad' lw 2 lc rgb 'blue'" :-
    @gp :- model(:domain) model(:narrow_lines) "w l t 'Narrow' lw 2 lc rgb 'brown'" :-
    @gp :- model(:domain) model(:unknown) "w l t 'Unknown' lc rgb 'purple'"
end
