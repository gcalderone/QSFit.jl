# Uncomment the following to implement a custom recipe:
# import QSFit: options, line_components, known_spectral_lines, fit!

function default_options(::Type{T}) where T <: DefaultRecipe
    out = OrderedDict{Symbol, Any}()
    out[:host_template] = "Ell5"
    out[:wavelength_range] = [1215, 7.3e3]
    out[:line_minimum_coverage] = 0.6
    out[:use_host_template] = true
    out[:use_balmer] = true
    out[:use_ironuv] = true
    out[:use_ironopt] = true
    out[:use_OIII_5007_bw] = false
    out[:n_unk] = 10
    out[:unk_avoid] = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]
    return out
end


function line_components(::QSO{T}, line::BroadBaseLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return [line.name => comp]
end

function line_components(::QSO{T}, line::BroadLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3
    return [line.name => comp]
end

function line_components(::QSO{T}, line::NarrowLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = 2e3
    comp.voff.low  = -1e3
    comp.voff.high =  1e3
    return [line.name => comp]
end

function line_components(::QSO{T}, line::CombinedLine) where T <: DefaultRecipe
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

function line_components(::QSO{T}, line::UnkLine) where T <: DefaultRecipe
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


function known_spectral_lines(source::QSO{T}) where T <: DefaultRecipe
    list = Vector{AbstractSpectralLine}()
    #push!(list,CombinedLine(  :Lyb          , 1026.0  ))
    push!(list, CombinedLine(  :Lya          , 1215.24 ))
    push!(list, NarrowLine(    :NV_1241      , 1240.81 ))
    push!(list, BroadLine(     :OI_1306      , 1305.53 ))
    push!(list, BroadLine(     :CII_1335     , 1335.31 ))
    push!(list, BroadLine(     :SiIV_1400    , 1399.8  ))
    push!(list, CombinedLine(  :CIV_1549     , 1549.48 ))
    #push!(list,BroadLine(     :HeII         , 1640.4  ))
    #push!(list,BroadLine(     :OIII         , 1665.85 ))
    #push!(list,BroadLine(     :AlIII        , 1857.4  ))
    push!(list, BroadLine(     :CIII_1909    , 1908.734))
    push!(list, BroadLine(     :CII          , 2326.0  ))
    push!(list, BroadLine(     :F2420        , 2420.0  ))
    push!(list, CombinedLine(  :MgII_2798    , 2799.117))
    #push!(list,NarrowLine(    :NeVN         , 3346.79 ))
    push!(list, NarrowLine(    :NeVI_3426    , 3426.85 ))
    push!(list, NarrowLine(    :OII_3727     , 3729.875))
    push!(list, NarrowLine(    :NeIII_3869   , 3869.81 ))
    push!(list, BroadLine(     :Hd           , 4102.89 ))
    push!(list, BroadLine(     :Hg           , 4341.68 ))
    #push!(list,NarrowLine(    :OIII_4363    , 4363.00 ))  # TODO: Check wavelength is correct
    push!(list, BroadLine(     :HeII         , 4686.   ))
    push!(list, CombinedLine(  :Hb           , 4862.68 ))
    push!(list, NarrowLine(    :OIII_4959    , 4960.295))
    push!(list, NarrowLine(    :OIII_5007    , 5008.240))
    if source.options[:use_OIII_5007_bw]
        push!(list, NarrowLine( :OIII_5007_bw, 5008.240))
    end
    push!(list, BroadLine(     :HeI_5876     , 5877.30 ))
    push!(list, NarrowLine(    :OI_6300      , 6300.00 ))  # TODO: Check wavelength is correct
    push!(list, NarrowLine(    :OI_6364      , 6364.00 ))  # TODO: Check wavelength is correct
    push!(list, NarrowLine(    :NII_6549     , 6549.86 ))
    push!(list, CombinedLine(  :Ha           , 6564.61 ))
    push!(list, BroadBaseLine( :Ha_base      , 6564.61 ))
    push!(list, NarrowLine(    :NII_6583     , 6585.27 ))
    push!(list, NarrowLine(    :SII_6716     , 6718.29 ))
    push!(list, NarrowLine(    :SII_6731     , 6732.67 ))
    return list
end


function fit(source::QSO{T}; dataset=1) where T <: DefaultRecipe
    elapsed = time()
    mzer = cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # Initialize components and guess initial values
    λ = source.domain[dataset][1]
    cont_components = [:qso_cont]
    model = Model(source.domain[dataset], :Continuum => Reducer(sum, cont_components),
                  :qso_cont => QSFit.powerlaw(3000))
    c = model[:qso_cont]
    c.norm.val = interpol(source.data[dataset].val, λ, c.x0.val)
    c.alpha.val = -1.8

    # Host galaxy template
    if source.options[:use_host_template]
        push!(cont_components, :galaxy)
        add!(model, :Continuum => Reducer(sum, cont_components),
             :galaxy => QSFit.hostgalaxy(source.options[:host_template]))
        model[:galaxy].norm.val = interpol(source.data[dataset].val, λ, 5500)
    end

    # Balmer continuum and pseudo-continuum
    if source.options[:use_balmer]
        push!(cont_components, :balmer)
        add!(model, :Continuum => Reducer(sum, cont_components),
             :balmer => QSFit.balmercont(0.1, 0.5))
        c = model[:balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 1.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.3
        c.ratio.high = 1
        patch!(model) do m
            m[:balmer].norm *= m[:qso_cont].norm
        end
    end

    bestfit = fit!(model, source.data, minimizer=mzer);  show(bestfit)

    # QSO continuum renormalization
    freeze(model, :qso_cont)
    c = model[:qso_cont]
    println(source.log, "Cont. norm. (before): ", c.norm.val)
    check_fraction = -1.
    last_fraction = check_fraction
    yy = source.data[dataset].val;
    ee = source.data[dataset].unc;
    while true
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

    freeze(model, :qso_cont)
    source.options[:use_host_template]  &&  freeze(model, :galaxy)
    source.options[:use_balmer]         &&  freeze(model, :balmer)
    evaluate(model)

    # Fit iron templates
    iron_components = Vector{Symbol}()
    if source.options[:use_ironuv]
        add!(model, :ironuv => QSFit.ironuv(3000))
        model[:ironuv].norm.val = 0.5
        push!(iron_components, :ironuv)
    end
    if source.options[:use_ironopt]
        add!(model,
             :ironoptbr => QSFit.ironopt_broad(3000),
             :ironoptna => QSFit.ironopt_narrow(500))
        model[:ironoptbr].norm.val = 0.5
        model[:ironoptna].norm.val = 0.0
        freeze(model, :ironoptna)  # will be freed during last run
        push!(iron_components, :ironoptbr, :ironoptna)
    end
    if length(iron_components) > 0
        add!(model, :Iron => Reducer(sum, iron_components))
        add!(model, :main => Reducer(sum, [:Continuum, :Iron]))
        evaluate(model)
        bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)
    else
        add!(model, :Iron => Reducer(() -> [0.], Symbol[]))
        add!(model, :main => Reducer(sum, [:Continuum, :Iron]))
    end
    source.options[:use_ironuv]   &&  freeze(model, :ironuv)
    source.options[:use_ironopt]  &&  freeze(model, :ironoptbr)
    source.options[:use_ironopt]  &&  freeze(model, :ironoptna)
    evaluate(model)

    # Add emission lines
    line_names = collect(keys(source.line_names[dataset]))
    line_groups = unique(collect(values(source.line_names[dataset])))
    add!(model, source.line_comps[dataset])
    for (group, lnames) in invert_dictionary(source.line_names[dataset])
        add!(model, group  => Reducer(sum, lnames))
    end
    add!(model, :main => Reducer(sum, [:Continuum, :Iron, line_groups...]))

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
    y = source.data[dataset].val - model()
    for cname in line_names
        c = model[cname]
        yatline = interpol(y, λ, c.center.val)
        c.norm.val = 1.
        c.norm.val = abs(yatline) / QSFit.maxvalue(model[cname])
    end

    # Patch parameters
    if  haskey(model.comps, :OIII_4959)  &&
        haskey(model.comps, :OIII_5007)
        model[:OIII_4959].voff.fixed = true
        patch!(model) do m
            m[:OIII_4959].voff = m[:OIII_5007].voff
        end
    end
    if  haskey(model.comps, :NII_6549)  &&
        haskey(model.comps, :NII_6583)
        model[:NII_6549].voff.fixed = true
        patch!(model) do m
            m[:NII_6549].voff = m[:NII_6583].voff
        end
    end
    if  haskey(model.comps, :OIII_5007_bw)  &&
        haskey(model.comps, :OIII_5007)
        patch!(model) do m
            m[:OIII_5007_bw].voff += m[:OIII_5007].voff
            m[:OIII_5007_bw].fwhm += m[:OIII_5007].fwhm
        end
    end
    bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)
    for lname in line_names
        freeze(model, lname)
    end

    # Add unknown lines
    if source.options[:n_unk] > 0
        tmp = OrderedDict{Symbol, Any}()
        for j in 1:source.options[:n_unk]
            tmp[Symbol(:unk, j)] = line_components(source, UnkLine())[1][2]
            tmp[Symbol(:unk, j)].norm.val = 0
        end
        add!(model, :UnkLines => Reducer(sum, collect(keys(tmp))), tmp)
        add!(model, :main => Reducer(sum, [:Continuum, :Iron, line_groups..., :UnkLines]))
        evaluate(model)
        for j in 1:source.options[:n_unk]
            freeze(model, Symbol(:unk, j))
        end
        evaluate(model)

        # Set "unknown" line center wavelength where there is a maximum in
        # the fit residuals, and re-run a fit.
        λunk = Vector{Float64}()
        while true
            (length(λunk) >= source.options[:n_unk])  &&  break
            evaluate(model)
            Δ = (source.data[dataset].val - model()) ./ source.data[dataset].unc

            # Avoid considering again the same region (within 1A)
            for l in λunk
                Δ[findall(abs.(l .- λ) .< 1)] .= 0.
            end

            # Avoidance regions
            for rr in source.options[:unk_avoid]
                Δ[findall(rr[1] .< λ .< rr[2])] .= 0.
            end

            # Do not add lines close to from the edges since these may
            # affect qso_cont fitting
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
    end
    evaluate(model)

    # Last run with all parameters free
    thaw(model, :qso_cont)
    source.options[:use_host_template]  &&  thaw(model, :galaxy)
    source.options[:use_balmer]         &&  thaw(model, :balmer)
    source.options[:use_ironuv]         &&  thaw(model, :ironuv)
    source.options[:use_ironopt]        &&  thaw(model, :ironoptbr)
    source.options[:use_ironopt]        &&  thaw(model, :ironoptna)

    for lname in line_names
        thaw(model, lname)
    end
    for j in 1:source.options[:n_unk]
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

    populate_metadata!(source, model)
    return (model, bestfit)
end
