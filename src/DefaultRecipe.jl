# ====================================================================
# Default recipe

function options(::QSO)
    out = OrderedDict{Symbol, Any}()
    out[:host_template] = "Ell5"
    out[:wavelength_range] = [1215, 7.3e3]
    out[:line_minimum_coverage] = 0.6
    return out
end


function line_components(::QSO, line::BroadBaseLine)
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return [line.name => comp]
end

function line_components(::QSO, line::BroadLine)
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3
    return [line.name => comp]
end

function line_components(::QSO, line::NarrowLine)
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = 2e3
    comp.voff.low  = -1e3
    comp.voff.high =  1e3
    return [line.name => comp]
end

function line_components(::QSO, line::CombinedLine)
    br = SpecLineGauss(line.λ)
    br.fwhm.val  = 5e3
    br.fwhm.low  = 900
    br.fwhm.high = 1.5e4
    br.voff.low  = -3e3
    br.voff.high =  3e3

    na = SpecLineGauss(line.λ)
    na.fwhm.val  = 5e2
    na.fwhm.low  = 100
    na.fwhm.high = 1e3
    na.voff.low  = -1e3
    na.voff.high =  1e3

    return [Symbol(:br_, line.name) => br
            Symbol(:na_, line.name) => na]
end

function line_components(::QSO, line::UnkLine)
    comp = SpecLineGauss(5e3)
    comp.center.fixed = false
    comp.center.low = 0
    comp.center.high = Inf
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 600
    comp.fwhm.high = 1e4
    comp.voff.fixed = true
    return [:Unk => comp]
end


function known_spectral_lines(::QSO)
    list = Vector{AbstractSpectralLine}()
    #push!(list,CombinedLine(:Lyb          , 1026.0  ))
    push!(list, CombinedLine( :Lya          , 1215.24 ))
    push!(list, NarrowLine(   :NV_1241      , 1240.81 ))
    push!(list, BroadLine(    :OI_1306      , 1305.53 ))
    push!(list, BroadLine(    :CII_1335     , 1335.31 ))
    push!(list, BroadLine(    :SiIV_1400    , 1399.8  ))
    push!(list, CombinedLine( :CIV_1549     , 1549.48 ))
    #push!(list,BroadLine(    :HeII         , 1640.4  ))
    #push!(list,BroadLine(    :OIII         , 1665.85 ))
    #push!(list,BroadLine(    :AlIII        , 1857.4  ))
    push!(list, BroadLine(    :CIII_1909    , 1908.734))
    push!(list, BroadLine(    :CII          , 2326.0  ))
    push!(list, BroadLine(    :F2420        , 2420.0  ))
    push!(list, CombinedLine( :MgII_2798    , 2799.117))
    #push!(list,NarrowLine(   :NeVN         , 3346.79 ))
    push!(list, NarrowLine(   :NeVI_3426    , 3426.85 ))
    push!(list, NarrowLine(   :OII_3727     , 3729.875))
    push!(list, NarrowLine(   :NeIII_3869   , 3869.81 ))
    push!(list, BroadLine(    :Hd           , 4102.89 ))
    push!(list, BroadLine(    :Hg           , 4341.68 ))
    #push!(list,NarrowLine(   :OIII_4363    , 4363.00 ))  # TODO: Check wavelength is correct
    push!(list, BroadLine(    :HeII         , 4686.   ))
    push!(list, CombinedLine( :Hb           , 4862.68 ))
    push!(list, NarrowLine(   :OIII_4959    , 4960.295))
    push!(list, NarrowLine(   :OIII_5007    , 5008.240))
    push!(list, NarrowLine(   :OIII_5007_bw , 5008.240))
    push!(list, BroadLine(    :HeI_5876     , 5877.30 ))
    push!(list, NarrowLine(   :OI_6300      , 6300.00 ))  # TODO: Check wavelength is correct
    push!(list, NarrowLine(   :OI_6364      , 6364.00 ))  # TODO: Check wavelength is correct
    push!(list, NarrowLine(   :NII_6549     , 6549.86 ))
    push!(list, CombinedLine( :Ha           , 6564.61 ))
    push!(list, BroadBaseLine(:Ha_base      , 6564.61 ))
    push!(list, NarrowLine(   :NII_6583     , 6585.27 ))
    push!(list, NarrowLine(   :SII_6716     , 6718.29 ))
    push!(list, NarrowLine(   :SII_6731     , 6732.67 ))
    return list
end


function GFit.fit!(source::QSO)
    @assert length(source.data) == 1
    elapsed = time()
    mzer = cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # Initialize components and guess initial values
    λ = source.domain[1][1]
    model = Model(source.domain[1],
                  :Broadband => Reducer(sum, [:continuum, :galaxy, :balmer]),
                  :galaxy    => QSFit.hostgalaxy(options(source)[:host_template]),
                  :continuum => QSFit.cutoff_powerlaw(1216),
                  :balmer    => QSFit.balmercont(0.1, 0.5))

    # Galaxy
    model[:galaxy].norm.val = interpol(source.data[1].val, λ, 5500)

    # Continuum
    c = model[:continuum]
    c.norm.val = interpol(source.data[1].val, λ, c.x0.val)
    c.alpha.val = -1.8
    c.beta.low = -50
    c.beta.val = -5

    # Balmer continuum
    c = model[:balmer]
    c.norm.val  = 0.1
    c.norm.fixed = false
    c.norm.high = 1.5
    c.ratio.val = 0.5
    c.ratio.fixed = false
    c.ratio.low  = 0.3
    c.ratio.high = 1
    patch!(model) do m
        m[:balmer].norm *= m[:continuum].norm
    end
    bestfit = fit!(model, source.data, minimizer=mzer);  show(bestfit)

    # Continuum renormalization
    freeze(model, :continuum)
    c = model[:continuum]
    println(source.log, "Cont. norm. (before): ", c.norm.val)
    check_fraction = -1.
    last_fraction = check_fraction
    yy = source.data[1].val;
    ee = source.data[1].unc;
    while true
        #global last_fraction
        mm = model()
        residuals = (mm .- yy) ./ ee
        check_fraction = count(residuals .< 0) / length(residuals)
        (last_fraction == check_fraction)  &&  break
        last_fraction = check_fraction
        (check_fraction > 0.9)  &&  break
        c.norm.val *= 0.99
        evaluate(model)
    end
    println(source.log, "Cont. norm. (after) : ", c.norm.val)

    freeze(model, :galaxy)
    freeze(model, :continuum)
    freeze(model, :balmer)
    evaluate(model)

    # Fit iron template
    add!(model, :Iron => Reducer(sum, [:ironuv, :ironoptbr, :ironoptna]),
         :ironuv    => QSFit.ironuv(3000),
         :ironoptbr => QSFit.ironopt_broad(3000),
         :ironoptna => QSFit.ironopt_narrow(500))
    add!(model, :main  => Reducer(sum, [:Broadband, :Iron]))
    evaluate(model)
    model[:ironuv].norm.val = 0.5
    model[:ironoptbr].norm.val = 0.5
    model[:ironoptna].norm.val = 0.0
    freeze(model, :ironoptna)  # will be freed during last run
    evaluate(model)
    bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)

    freeze(model, :ironuv)
    freeze(model, :ironoptbr)
    freeze(model, :ironoptna)

    # Add emission lines
    line_names = collect(keys(source.line_names[1]))
    line_groups = unique(collect(values(source.line_names[1])))
    add!(model, source.line_comps[1])
    for (group, lnames) in invert_dictionary(source.line_names[1])
        add!(model, group  => Reducer(sum, lnames))
    end
    add!(model, :main => Reducer(sum, [:Broadband, :Iron, line_groups...]))

    if haskey(model.comps, :MgII_2798)
        model[:MgII_2798].voff.low  = -1000
        model[:MgII_2798].voff.high =  1000
    end
    if haskey(model.comps, :OIII_5007_bw)
        model[:OIII_5007_bw].fwhm.val  = 500
        model[:OIII_5007_bw].fwhm.low  = 1e2
        model[:OIII_5007_bw].fwhm.high = 1e3
        model[:OIII_5007_bw].voff.low  = 0
        model[:OIII_5007_bw].voff.high = 2e3
    end

    # Guess values
    evaluate(model)
    y = source.data[1].val - model()
    for cname in line_names
        c = model[cname]
        yatline = interpol(y, λ, c.center.val)
        c.norm.val = 1.
        c.norm.val = abs(yatline) / QSFit.maxvalue(model[cname])
    end

    model[:OIII_4959].voff.fixed = true
    patch!(model) do m
        m[:OIII_4959].voff = m[:OIII_5007].voff
        m[:OIII_5007_bw].voff += m[:OIII_5007].voff
        m[:OIII_5007_bw].fwhm += m[:OIII_5007].fwhm
    end
    bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)
    for lname in line_names
        freeze(model, lname)
    end

    # Add unknown lines
    nunk = 10
    tmp = OrderedDict{Symbol, Any}()
    for j in 1:nunk
        tmp[Symbol(:unk, j)] = line_components(source, UnkLine())[1][2]
        tmp[Symbol(:unk, j)].norm.val = 0
    end
    add!(model, :UnkLines => Reducer(sum, collect(keys(tmp))), tmp)
    add!(model, :main => Reducer(sum, [:Broadband, :Iron, line_groups..., :UnkLines]))
    evaluate(model)
    for j in 1:nunk
        freeze(model, Symbol(:unk, j))
    end
    evaluate(model)

    # Set "unknown" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λunk = Vector{Float64}()
    while true
        (length(λunk) >= nunk)  &&  break
        evaluate(model)
        Δ = (source.data[1].val - model()) ./ source.data[1].unc

        # Avoid considering again the same region (within 1A)
        for l in λunk
            Δ[findall(abs.(l .- λ) .< 1)] .= 0.
        end

        # Avoidance regions
        Δ[abs.(λ .- 4863) .<  50] .= 0.
        Δ[abs.(λ .- 6565) .< 150] .= 0.

        # Do not add lines close to from the edges since these may
        # affect continuum fitting
        Δ[findall((λ .< minimum(λ)*1.02)  .|
                  (λ .> maximum(λ)*0.98))] .= 0.
        iadd = argmax(Δ)
        (Δ[iadd] <= 0)  &&  break  # No residual is greater than 0, skip further residuals....
        push!(λunk, λ[iadd])

        cname = Symbol(:unk, length(λunk))
        model[cname].norm.val = 1.
        model[cname].center.val  = λ[iadd]
        model[cname].center.low  = λ[iadd] - λ[iadd]/10. # allow to move 10%
        model[cname].center.high = λ[iadd] + λ[iadd]/10.

        thaw(model, cname)
        bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)
        freeze(model, cname)
    end
    evaluate(model)

    # Last run with all parameters free
    thaw(model, :galaxy)
    thaw(model, :continuum)
    thaw(model, :balmer)
    thaw(model, :ironuv)
    thaw(model, :ironoptbr)
    thaw(model, :ironoptna)

    for lname in line_names
        thaw(model, lname)
    end
    for j in 1:nunk
        cname = Symbol(:unk, j)
        if model[cname].norm.val > 0
            thaw(model, cname)
        else
            freeze(model, cname)
        end
    end

    bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)
    elapsed = time() - elapsed
    println("Elapsed time: $elapsed s")
    return (model, bestfit)
end
