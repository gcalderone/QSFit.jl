function plot(qsfit, model)
    @gp "set bars 0" :-
    @gp :- qsfit.domain[1][1] qsfit.data[1].val qsfit.data[1].unc "w yerr t 'Data' pt 0 lc rgb 'black'" :-
    @gp :- model(:domain) model() "w l t 'Model' lw 2 lc rgb 'red'" :-
    @gp :- model(:domain) model(:continuum) "w l t 'Cont'" :-
    qsfit.options.use_galaxy  &&  (@gp :- model(:domain) model(:galaxy) "w l t 'Host'" :-)
    qsfit.options.use_balmer  &&  (@gp :- model(:domain) model(:balmer) "w l t 'Balmer'" :-)
    qsfit.options.use_ironuv  &&  (@gp :- model(:domain) model(:ironuv) "w l t 'IronUV'" :-)
    qsfit.options.use_ironopt &&  (@gp :- model(:domain) model(:ironoptbr) "w l t 'IronOptB'" :-)
    qsfit.options.use_ironopt &&  (@gp :- model(:domain) model(:ironoptna) "w l t 'IronOptN'" :-)
    @gp :- model(:domain) model(:narrow_lines) "w l t 'Narrow'" :-
    @gp :- model(:domain) model(:broad_lines) "w l t 'Broad'" :-
    @gp :- model(:domain) model(:unknown) "w l t 'Unknown'"
end
