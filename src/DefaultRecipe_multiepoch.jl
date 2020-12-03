T(i, args...) = Symbol(:T, i, :_, args...)
macro T(i, a)
    esc(:(T($i, $a)))
end
macro T(i, a, b)
    esc(:(T($i, $a, $b)))
end

calibsum(calib, args...) = calib .* (.+(args...))

function multiepoch_fit(source::QSO{TRecipe}; ref_id=1) where TRecipe <: DefaultRecipe
    Nspec = length(source.domain)

    elapsed = time()
    mzer = cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # Arrays containing best fit values to be constrained across epochs
    galaxy_best = Vector{Float64}()
    galaxy_unc  = Vector{Float64}()
    OIII_best = Vector{Float64}()
    OIII_unc  = Vector{Float64}()

    # Initialize components and guess initial values
    preds = Vector{Prediction}()
    for id in 1:Nspec
        λ = source.domain[id][1]
        pred = Prediction(source.domain[id], :Continuum => Reducer(sum, T.(id, [:qso_cont])),
                          T(id, :qso_cont) => qso_cont_component(TRecipe))
        push!(preds, pred)
        c = pred[@T id :qso_cont]
        c.norm.val = interpol(source.data[id].val, λ, c.x0.val)
    end
    model = Model(preds)

    for id in 1:Nspec
        λ = source.domain[id][1]

        # Host galaxy template
        if source.options[:use_host_template]
            add!(model, id=id, :Continuum => Reducer(sum, T.(id, [:qso_cont, :galaxy])),
                 T(id, :galaxy) => QSFit.hostgalaxy(source.options[:host_template]))
            model[@T id :galaxy].norm.val = interpol(source.data[id].val, λ, 5500)
        end

        # Balmer continuum and pseudo-continuum
        if source.options[:use_balmer]
            tmp = [:qso_cont, :balmer]
            source.options[:use_host_template]  &&  push!(tmp, :galaxy)
            add!(model, id=id, :Continuum => Reducer(sum, T.(id, tmp)),
                 T(id, :balmer) => QSFit.balmercont(0.1, 0.5))
            c = model[@T id :balmer]
            c.norm.val  = 0.1
            c.norm.fixed = false
            c.norm.high = 0.5
            c.ratio.val = 0.5
            c.ratio.fixed = false
            c.ratio.low  = 0.3
            c.ratio.high = 1
            patch!(model) do m
                m[@T id :balmer].norm *= m[@T id :qso_cont].norm
            end
        end
    end
    for id in 1:Nspec
        bestfit = fit!(model, id=id, source.data, minimizer=mzer);  show(bestfit)
        push!(galaxy_best, bestfit[@T id :galaxy].norm.val)
        push!(galaxy_unc , bestfit[@T id :galaxy].norm.unc)
    end

    # QSO continuum renormalization
    for id in 1:Nspec
        freeze(model, @T id :qso_cont)
        c = model[@T id :qso_cont]
        println(source.log, "Cont. norm. (before): ", c.norm.val)
        while true
            residuals = (model(id=id) - source.data[id].val) ./ source.data[id].unc
            (count(residuals .< 0) / length(residuals) > 0.9)  &&  break
            c.norm.val *= 0.99
            evaluate(model)
        end
        println(source.log, "Cont. norm. (after) : ", c.norm.val)

        freeze(model, @T id :qso_cont)
        source.options[:use_host_template]  &&  freeze(model, @T id :galaxy)
        source.options[:use_balmer]         &&  freeze(model, @T id :balmer)
    end
    evaluate(model)

    # Fit iron template
    for id in 1:Nspec
        iron_components = Vector{Symbol}()
        if source.options[:use_ironuv]
            add!(model, id=id, T(id, :ironuv) => QSFit.ironuv(3000))
            model[@T id :ironuv].norm.val = 0.5
            push!(iron_components, :ironuv)
        end
        if source.options[:use_ironopt]
            add!(model, id=id,
                 T(id, :ironoptbr) => QSFit.ironopt_broad(3000),
                 T(id, :ironoptna) => QSFit.ironopt_narrow(500))
            model[@T id :ironoptbr].norm.val = 0.5
            model[@T id :ironoptna].norm.val = 0.0
            freeze(model, @T id :ironoptna)  # will be freed during last run
            push!(iron_components, :ironoptbr, :ironoptna)
        end
        if length(iron_components) > 0
            add!(model, id=id, :Iron => Reducer(sum, T.(id, iron_components)))
            add!(model, id=id, :main => Reducer(sum, [:Continuum, :Iron]))
            evaluate(model)
            bestfit = fit!(model, id=id, source.data, minimizer=mzer); show(bestfit)
        else
            add!(model, id=id, :Iron => Reducer(() -> [0.], Symbol[]))
            add!(model, id=id, :main => Reducer(sum, [:Continuum, :Iron]))
        end
        source.options[:use_ironuv]   &&  freeze(model, @T id :ironuv)
        source.options[:use_ironopt]  &&  freeze(model, @T id :ironoptbr)
        source.options[:use_ironopt]  &&  freeze(model, @T id :ironoptna)
    end
    evaluate(model)

    # Add emission lines
    for id in 1:Nspec
        λ = source.domain[id][1]

        line_names = collect(keys(source.line_names[id]))
        line_groups = unique(collect(values(source.line_names[id])))
        for (lname, lcomp) in source.line_comps[id]
            add!(model, id=id, T(id, lname) => lcomp)
        end
        for (group, lnames) in invert_dictionary(source.line_names[id])
            add!(model, id=id, group => Reducer(sum, T.(id, lnames)))
        end
        add!(model, id=id, :main => Reducer(sum, [:Continuum, :Iron, line_groups...]))

        if haskey(model.comps, @T id :MgII_2798)
            model[@T id :MgII_2798].voff.low  = -1000
            model[@T id :MgII_2798].voff.high =  1000
        end
        if haskey(model.comps, @T id :OIII_5007_bw)
            model[@T id :OIII_5007_bw].fwhm.val  = 500
            model[@T id :OIII_5007_bw].fwhm.low  = 1e2
            model[@T id :OIII_5007_bw].fwhm.high = 1e3
            model[@T id :OIII_5007_bw].voff.low  = 0
            model[@T id :OIII_5007_bw].voff.high = 2e3
        end

        # Guess values
        evaluate(model)
        y = source.data[id].val - model(id=id)
        for cname in line_names
            c = model[@T id cname]
            yatline = interpol(y, λ, c.center.val)
            c.norm.val = 1.
            c.norm.val = abs(yatline) / QSFit.maxvalue(model[@T id cname])
        end

        # Patch parameters
        if  haskey(model.comps, @T id :OIII_4959)  &&
            haskey(model.comps, @T id :OIII_5007)
            model[@T id :OIII_4959].voff.fixed = true
            patch!(model) do m
                m[@T id :OIII_4959].voff = m[@T id :OIII_5007].voff
            end
        end
        if  haskey(model.comps, @T id :NII_6549)  &&
            haskey(model.comps, @T id :NII_6583)
            model[@T id :NII_6549].voff.fixed = true
            patch!(model) do m
                m[@T id :NII_6549].voff = m[@T id :NII_6583].voff
            end
        end
        if  haskey(model.comps, @T id :OIII_5007_bw)  &&
            haskey(model.comps, @T id :OIII_5007)
            patch!(model) do m
                m[@T id :OIII_5007_bw].voff += m[@T id :OIII_5007].voff
                m[@T id :OIII_5007_bw].fwhm += m[@T id :OIII_5007].fwhm
            end
        end
        bestfit = fit!(model, id=id, source.data, minimizer=mzer); show(bestfit)
        push!(OIII_best, bestfit[@T id :OIII_5007].norm.val)
        push!(OIII_unc , bestfit[@T id :OIII_5007].norm.unc)

        for lname in line_names
            freeze(model, @T id lname)
        end
    end

    # Add unknown lines
    if source.options[:n_unk] > 0
        for id in 1:Nspec
            tmp = OrderedDict([@T(id, :unk, j) => line_components(TRecipe, UnkLine())[1][2] for j in 1:source.options[:n_unk]])
            add!(model, id=id, :UnkLines => Reducer(sum, collect(keys(tmp))), tmp)
            line_groups = unique(collect(values(source.line_names[id])))
            add!(model, id=id, :main => Reducer(sum, [:Continuum, :Iron, line_groups..., :UnkLines]))
            evaluate(model)
            for j in 1:source.options[:n_unk]
                freeze(model, @T id :unk j)
            end
        end
    else
        for id in 1:Nspec
            add!(model, id=id, :UnkLines => Reducer(() -> [0.], Symbol[]))
            line_groups = unique(collect(values(source.line_names[id])))
            add!(model, id=id, :main => Reducer(sum, [:Continuum, :Iron, line_groups..., :UnkLines]))
        end
    end
    evaluate(model)


    # Set "unknown" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    for id in 1:Nspec
        λ = source.domain[id][1]
        λunk = Vector{Float64}()
        while true
            (length(λunk) >= source.options[:n_unk])  &&  break
            evaluate(model)
            Δ = (source.data[id].val - model(id=id)) ./ source.data[id].unc

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

            cname = @T id :unk length(λunk)
            model[cname].norm.val = 1.
            model[cname].center.val  = λ[iadd]
            model[cname].center.low  = λ[iadd] - λ[iadd]/10. # allow to shift 10%
            model[cname].center.high = λ[iadd] + λ[iadd]/10.

            thaw(model, cname)
            bestfit = fit!(model, id=id, source.data, minimizer=mzer); show(bestfit)
            freeze(model, cname)
        end
    end
    evaluate(model)

    # ----------------------------------------------------------------
    # Constrain component normalization across epochs.  Note:
    # reference spectrum must have reliable estimation of all common
    # components
    rval = [galaxy_best[ref_id], OIII_best[ref_id]]
    runc = [galaxy_unc[ ref_id], OIII_unc[ ref_id]]
    @assert count((rval .!= 0)  .&  (runc .!= 0)) == 2

    # Estimate calibration in all epochs w.r.t. reference epoch
    for id in 1:Nspec
        if id != ref_id
            val = [galaxy_best[id], OIII_best[id]]
            unc = [galaxy_unc[ id], OIII_unc[ id]]
            j = findall((val .!= 0)  .&  (unc .!= 0))

            R = val[j] ./ rval[j]
            R_unc = (unc[j] ./ val[j] .+ runc[j] ./ rval[j]) .* R
            ratio = sum(R ./ R_unc) ./ sum(1 ./ R_unc)
            for cname in keys(preds[id].cevals)
                if :norm in propertynames(model[cname])
                    model[cname].norm.val /= ratio
                end
            end
            model[@T id :galaxy].norm.fixed = true
            model[@T id :OIII_5007].norm.fixed = true
        else
            ratio = 1.
        end
        line_groups = unique(collect(values(source.line_names[id])))
        add!(model, id=id,
             :main => Reducer(calibsum, [@T(id, :calib), :Continuum, :Iron, line_groups..., :UnkLines]),
             T(id, :calib) => ratio)
    end
    evaluate(model)
    model[@T ref_id :calib].par.fixed = true

    for id in 1:Nspec
        if id != ref_id
            patch!(model) do m
                m[@T id    :galaxy].norm = m[@T ref_id    :galaxy].norm
                m[@T id :OIII_5007].norm = m[@T ref_id :OIII_5007].norm
            end
        end
    end
    evaluate(model)

    # Last run with all parameters free
    for id in 1:Nspec
        source.options[:use_host_template]  &&  thaw(model, @T id :qso_cont)
        source.options[:use_balmer]         &&  thaw(model, @T id :galaxy)
        source.options[:use_ironuv]         &&  thaw(model, @T id :balmer)
        source.options[:use_ironopt]        &&  thaw(model, @T id :ironoptbr)
        source.options[:use_ironopt]        &&  thaw(model, @T id :ironoptna)

        line_names = collect(keys(source.line_names[id]))
        for lname in line_names
            thaw(model, @T id lname)
        end
        for j in 1:source.options[:n_unk]
            cname = @T id :unk j
            if model[cname].norm.val > 0
                thaw(model, cname)
            else
                freeze(model, cname)
            end
        end
        if id != ref_id
            thaw(model, @T id :calib)  # parameter is fixed in preds[1]
        end
    end
    bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)

    # Disable "unknown" lines whose normalization uncertainty is larger
    # than 3 times the normalization
    needs_fitting = false
    for id in 1:Nspec
        for ii in 1:source.options[:n_unk]
            cname = @T id :unk ii
            model.cfixed[cname]  &&  continue
            if bestfit[cname].norm.val == 0.
                freeze(model, cname)
                needs_fitting = true
                @info "Disabling $cname (norm. = 0)"
            elseif bestfit[cname].norm.unc / bestfit[cname].norm.val > 3
                model[cname].norm.val = 0.
                freeze(model, cname)
                needs_fitting = true
                @info "Disabling $cname (unc. / norm. > 3)"
            end
        end
    end
    if needs_fitting
        bestfit = fit!(model, source.data, minimizer=mzer); show(bestfit)
    end

    elapsed = time() - elapsed
    println("Elapsed time: $elapsed s")
    populate_metadata!(source, model)
    return (model, bestfit)
end
