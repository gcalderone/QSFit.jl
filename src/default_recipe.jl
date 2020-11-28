# ====================================================================
# Default recipe

function default_broadline(::QSO, λ)
    out = SpecLineGauss(λ)

    out.fwhm.val  = 5e3
    out.fwhm.low  = 900
    out.fwhm.high = 1.5e4

    out.voff.low  = -3e3
    out.voff.high =  3e3
    return out
end


function default_narrowline(::QSO, λ)
    out = SpecLineGauss(λ)
    out.fwhm.val  = 5e2
    out.fwhm.low  = 100
    out.fwhm.high = 2e3

    out.voff.low  = -1e3
    out.voff.high =  1e3
    return out
end


function default_combinedlines(::QSO, λ)
    br = SpecLineGauss(λ)
    br.fwhm.val  = 5e3
    br.fwhm.low  = 900
    br.fwhm.high = 1.5e4

    br.voff.low  = -3e3
    br.voff.high =  3e3

    na = SpecLineGauss(λ)
    na.fwhm.val  = 5e2
    na.fwhm.low  = 100
    na.fwhm.high = 1e3

    na.voff.low  = -1e3
    na.voff.high =  1e3

    return [br, na]
end


function default_unknownline(::QSO)
    out = SpecLineGauss(5e3)
    out.center.fixed = false
    out.center.low = 0
    out.center.high = Inf

    out.fwhm.val  = 5e3
    out.fwhm.low  = 600
    out.fwhm.high = 1e4
    out.voff.fixed = true
    return out
end


function default_known_lines(::QSO)
    list = OrderedDict{Symbol, Tuple{Symbol, Float64}}()
    #list[:Lyb         ] = (:Combined, 1026.0  )
    list[:Lya          ] = (:Combined, 1215.24 )
    list[:NV_1241      ] = (:Narrow,   1240.81 )
    list[:OI_1306      ] = (:Broad,    1305.53 )
    list[:CII_1335     ] = (:Broad,    1335.31 )
    list[:SiIV_1400    ] = (:Broad,    1399.8  )
    list[:CIV_1549     ] = (:Combined, 1549.48 )
    # list[:HeII       ] = (:Broad,    1640.4  )
    # list[:OIII       ] = (:Broad,    1665.85 )
    # list[:AlIII      ] = (:Broad,    1857.4  )
    list[:CIII_1909    ] = (:Broad,    1908.734)
    list[:CII          ] = (:Broad,    2326.0  )
    list[:F2420        ] = (:Broad,    2420.0  )
    list[:MgII_2798    ] = (:Combined, 2799.117)
    # list[:NeVN       ] = (:Narrow,   3346.79 )
    list[:NeVI_3426    ] = (:Narrow,   3426.85 )
    list[:OII_3727     ] = (:Narrow,   3729.875)
    list[:NeIII_3869   ] = (:Narrow,   3869.81 )
    list[:Hd           ] = (:Broad,    4102.89 )
    list[:Hg           ] = (:Broad,    4341.68 )
    #list[:OIII_4363   ] = (:Narrow,    4363.00)  # TODO: Check wavelength is correct
    list[:HeII         ] = (:Broad,    4686.   )
    list[:Hb           ] = (:Combined, 4862.68 )
    list[:OIII_4959    ] = (:Narrow,   4960.295)
    list[:OIII_5007    ] = (:Narrow,   5008.240)
    list[:OIII_5007_bw ] = (:Narrow,   5008.240)
    list[:HeI_5876     ] = (:Broad,    5877.30 )
    list[:OI_6300      ] = (:Narrow,   6300.00 )  # TODO: Check wavelength is correct
    list[:OI_6364      ] = (:Narrow,   6364.00 )  # TODO: Check wavelength is correct
    list[:NII_6549     ] = (:Narrow,   6549.86 )
    list[:Ha           ] = (:Combined, 6564.61 )
    list[:Ha_base      ] = (:Broad,    6564.61 )
    list[:NII_6583     ] = (:Narrow,   6585.27 )
    list[:SII_6716     ] = (:Narrow,   6718.29 )
    list[:SII_6731     ] = (:Narrow,   6732.67 )

    # Expand combined lines into broad and narrow components
    out = OrderedDict{Symbol, Tuple{Symbol, Float64}}()
    for (lname, (ltype, lwave)) in list
        @assert ltype in [:Narrow, :Broad, :Combined]
        if ltype == :Combined
            out[Symbol(:na_, lname)] = (:Narrow, lwave)
            out[Symbol(:br_, lname)] = (:Broad , lwave)
        else
            out[lname] = (ltype, lwave)
        end
    end
    return out
end


function GFit.fit!(source::QSO)
    @assert length(source.data) == 1
    elapsed = time()
    mzer = cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # Initialize components and predictions
    λ = source.domain[1][1]
    model = Model(source.domain[1],
                  :Broadband => Reducer(sum, [:continuum, :galaxy, :balmer]),
                  :galaxy    => QSFit.hostgalaxy("Ell5"),
                  :continuum => QSFit.cutoff_powerlaw(1216),
                  :balmer    => QSFit.balmercont(0.1, 0.5))

    # Guess initial values

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
    d = OrderedDict(zip(keys(source.broad_lines[1]), deepcopy(values(source.broad_lines[1]))))
    add!(model, :BroadLines  => Reducer(sum, collect(keys(d))) , d)
    d = OrderedDict(zip(keys(source.narrow_lines[1]), deepcopy(values(source.narrow_lines[1]))))
    add!(model, :NarrowLines => Reducer(sum, collect(keys(d))) , d)
    add!(model, :main => Reducer(sum, [:Broadband, :Iron, :BroadLines, :NarrowLines]))
    line_names = [collect(keys(source.broad_lines[1])); collect(keys(source.narrow_lines[1]))]

    if haskey(model.comps, :Ha_base)
        model[:Ha_base].fwhm.val  = 2e4
        model[:Ha_base].fwhm.low  = 1e4
        model[:Ha_base].fwhm.high = 3e4
        model[:Ha_base].voff.val = 0
        model[:Ha_base].voff.fixed = true
    end
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
        tmp[Symbol(:unk, j)] = default_unknownline(source)
        tmp[Symbol(:unk, j)].norm.val = 0
    end
    add!(model, :UnkLines=>Reducer(sum, collect(keys(tmp))), tmp)
    add!(model, :main => Reducer(sum, [:Broadband, :Iron, :BroadLines, :NarrowLines, :UnkLines]))
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
