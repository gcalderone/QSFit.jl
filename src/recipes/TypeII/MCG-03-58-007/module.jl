import DataFitting: fit!

const showstep = true

mutable struct MCG0358007_Options
    # The wavelength range used to for fitting.  Wavelengths outside
    # the range are ignored.
    λ_range::NTuple{2, Float64}

    # Fraction of negative residuals after continuum re-normalization
    cont_negative_fraction::Float64

    # Name of host galaxy template
    galaxy_template::String

    # If true use Lorentzian (rather than Gaussian) profiles for
    # emission lines
    lorentzian::Bool

    function MCG0358007_Options()
        return new(
            (1210., 1e5), #λ_range
            0.9,    #cont_negative_fraction
            "Ell5", #       galaxy_template
            false)  #            lorentzian
    end
end



@quasiabstract struct MCG0358007 <: Source
    options::MCG0358007_Options
    function MCG0358007(args...; kw...)
        tmp = Source(args...; kw...)
        return new(getfield.(Ref(tmp), fieldnames(typeof(tmp)))..., MCG0358007_Options())
    end
end


function add_spec!(source::MCG0358007, data::Spectrum)
    empty!(source.lines)
    l = NarrowLine(        "na_NeIII_3869"  , 3869.81 );  push!(source.lines, l)
    l = NarrowLine(        "na_Hb"          , 4862.68 );  push!(source.lines, l)
    l = NarrowLine(        "na_OIII_4959"   , 4960.295);  push!(source.lines, l)
    l = NarrowLine(        "na_OIII_5007"   , 5008.240);  push!(source.lines, l)
    l = NarrowLine(        "na_OIII_5007_bw", 5008.240);  l.fwhm = 500; l.fwhm_limits = (1e2, 1e3); l.voff_limits = (0, 2e3); push!(source.lines, l)
    l = NarrowLine(        "na_NII_6549"    , 6549.86 );  push!(source.lines, l)
    l = NarrowLine(        "na_Ha"          , 6564.61 );  push!(source.lines, l)
    l = NarrowLine(        "na_NII_6583"    , 6585.27 );  push!(source.lines, l)
    l = NarrowLine(        "na_SII_6716"    , 6718.29 );  push!(source.lines, l)
    l = NarrowLine(        "na_SII_6731"    , 6732.67 );  push!(source.lines, l)
    println(source.log, "New data: " * data.label)
    println(source.log, "  good fraction:: " * string(goodfraction(data)))
    if goodfraction(data) < 0.5
        @error "Good fraction < 0.5"
    end

    λ = data.λ ./ (1 + source.z)
    ii = findall(λ .<= source.options.λ_range[1])
    data.good[ii] .= false

    dered = ccm_unred([1450, 3000, 5100.], source.ebv)
    println(source.log, "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.ebv)

    igood = findall(data.good)
    dom = Domain(data.λ[igood] ./ (1 + source.z))
    lum = Measures(data.flux[igood] .* dered[igood] .* source.flux2lum .* (1 + source.z),
                   data.err[ igood] .* dered[igood] .* source.flux2lum .* (1 + source.z))
    push!(source.domain, dom)
    push!(source.data, lum)
end


function fit!(source::MCG0358007)
    elapsed = time()
    mzer = cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # First dataset have a special meaning
    λ = source.domain[1][1]

    # Prepare model
    model = Model()
    for ii in 1:length(source.domain)
        add_dom!(model, source.domain[ii])
    end

    # List of component names
    cnames = Dict{Symbol,   Vector{Symbol}}()
    cnames[:broadband]    = Vector{Symbol}()
    cnames[:narrow_lines] = Vector{Symbol}()

    # AGN continuum
    add_comp!(model, :continuum => powerlaw((minimum(λ) + maximum(λ)) / 2.))
    push!(cnames[:broadband], :continuum)

    # Host galaxy
    add_comp!(model, :galaxy => hostgalaxy(source.options.galaxy_template))
    push!(cnames[:broadband], :galaxy)
    add_expr!(model, :broadband, Meta.parse(join(cnames[:broadband], ".+")), cmp=false)
    #add_expr!(model, :broadband, :(continuum .+ galaxy))

    # Emission lines
    for line in source.lines
        if line.enabled
            if isa(line, NarrowLine)  ||  isa(line, CombinedNarrowLine)
                cname = addEmLine(model, line)
                push!(cnames[:narrow_lines], cname)
                model[cname].norm.val = 0.
                model[cname].fixed = true
                if cname == :na_OIII_4959
                    model[cname].voff.expr = "na_OIII_5007_voff"
                    model[cname].voff.fixed = true
                end
                if cname == :na_OIII_5007_bw
                    model[cname].fwhm.expr = "na_OIII_5007_fwhm + na_OIII_5007_bw_fwhm"
                    model[cname].voff.expr = "na_OIII_5007_voff + na_OIII_5007_bw_voff"
                end
            end
        end
    end
    add_expr!(model, :narrow_lines, Meta.parse(join(cnames[:narrow_lines], ".+")), cmp=false)
    add_expr!(model, :known, :(broadband .+ narrow_lines), cmp=false)
    add_expr!(model, :all, :(known .* 1.))

    # AGN continuum
    let
        L = source.data[1].val
        model.continuum.norm.val = interpol(L, λ, model.continuum.x0.val)
        model.continuum.alpha.val = -1.7
        model.continuum.alpha.low  = -3
        model.continuum.alpha.high = 1 #--> -3, 1 in frequency
        model.continuum.alpha.fixed = true # avoid degeneracy with galaxy template
    end

    # Host galaxy
    let
        L = source.data[1].val
        model.galaxy.norm.val = interpol(L, λ, 5500)
        model.galaxy.norm.val = interpol(L, λ, maximum(λ))
        (model.galaxy.norm.val < 1e-4)  &&  (model.galaxy.norm.val = 1e-4)
    end
    show(model)
    bestfit = fit!(model, source.data, minimizer=mzer)
    showstep  &&   (show(model); show(bestfit); plot(source, model); readline())

    # Continuum renormalization
    let
        println(source.log, "Cont. norm. (before): ", model.continuum.norm.val)
        check_fraction = -1.
        last_fraction = check_fraction
        yy = source.data[1].val
        ee = source.data[1].unc
        while true
            mm = model(1)

            residuals = (mm .- yy) ./ ee
            check_fraction = length(findall(residuals .< 0)) / length(yy)
            (last_fraction == check_fraction)  &&  break
            last_fraction = check_fraction
            (check_fraction > source.options.cont_negative_fraction)  &&  break

            model.continuum.norm.val *= 0.99
            evaluate!(model)
        end
        println(source.log, "Cont. norm. (after) : ", model.continuum.norm.val)
    end
    evaluate!(model)
    showstep  &&   (show(model); show(bestfit); plot(source, model); readline())
    model.continuum.fixed = true
    model.galaxy.fixed = true

    # Emission lines
    let
        x = domain(model)
        y = source.data[1].val - model(:broadband)
        for line in source.lines
            if line.enabled
                if isa(line, NarrowLine)  ||  isa(line, CombinedNarrowLine)
                    yatline = interpol(y, x, line.λ)
                    cname = Symbol(line.label)
                    model[cname].norm.val = 1.
                    model[cname].norm.val = abs(yatline) / maxvalue(model[cname])
                    model[cname].fixed = false
                end
            end
        end

        bestfit = fit!(model, source.data, minimizer=mzer)
        showstep  &&   (show(model); show(bestfit); plot(source, model); readline())
        for line in source.lines
            if isa(line, NarrowLine)  ||  isa(line, CombinedNarrowLine)
                cname = Symbol(line.label)
                line.enabled  &&  (model[cname].fixed = true)
            end
        end
    end
    showstep  &&   (show(model); show(bestfit); plot(source, model); readline())

    # Last run with all parameters free
    model.continuum.fixed = false
    model.galaxy.fixed = false
    for cname in cnames[:narrow_lines];  model[cname].fixed = false;  end
    @info "last run"
    bestfit = fit!(model, source.data, minimizer=mzer)
    showstep  &&   (show(model); show(bestfit); plot(source, model); readline())

    elapsed = time() - elapsed
    @info elapsed
    showstep  &&   (show(model); show(bestfit); plot(source, model); readline())
    return (model, bestfit)
    # catch err
    # finally
    #     if source.log != stdout
    #         close(source.log)
    #     end
    #     return nothing
    # end
end


function plot(source::MCG0358007, model)
    @gp "set bars 0" "set grid" "set key horizontal" xr=extrema(domain(model)) :-
    @gp :- title=source.name * ", z=" * string(source.z) * ", E(B-V)=" * string(source.ebv) :-
    @gp :- xlabel="Rest frame wavelength [A]" ylabel="Lum. density [10^{42} erg s^{-1} A^{-1}]" :-
    @gp :- source.domain[1][1] source.data[1].val source.data[1].unc "w yerr t 'Data' pt 0 lc rgb 'black'" :-
    @gp :- domain(model) model() "w l t 'Model' lw 2 lc rgb 'orange'" :-
    @gp :- domain(model) model(:continuum) "w l t 'Cont' dt 2 lc rgb 'red'" :-
    @gp :- domain(model) model(:continuum) .+ model(:galaxy) "w l t 'Cont+Host' lc rgb 'red'" :-
    @gp :- domain(model) model(:narrow_lines) "w l t 'Narrow' lw 2 lc rgb 'brown'"


    @gp    :resid "set bars 0" "set grid" xr=extrema(domain(model)) :-
    @gp :- :resid title=source.name * ", z=" * string(source.z) * ", E(B-V)=" * string(source.ebv) :-
    @gp :- :resid xlabel="Rest frame wavelength [A]" ylabel="Residuals [{/Symbol s}]" :-
    resid = (source.data[1].val .- model()) ./ source.data[1].unc
    nn = length(resid)
    @gp :- :resid source.domain[1][1] resid fill(1., nn) "w p t 'Data' lc rgb 'black'" :-
    @gp :- :resid [extrema(source.domain[1][1])...] [0,0] "w line notitle dt 2 lw 2 lt rgb 'orange'" :-

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
    @gp :- :resid source.domain[1][1] fs "w l title 'Reduced cumul. {/Symbol c}^2' ls 1 lw 2 lt rgb 'red' axes x1y2"
end
