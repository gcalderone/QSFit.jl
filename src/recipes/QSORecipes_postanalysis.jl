function postanalysis(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe

    outm = Vector{OrderedDict{Symbol, Any}}()

    for ith in 1:nmodels(fp)
        out = OrderedDict{Symbol, Any}()
        model = getmodel(fp, ith)

        # Data stats
        dict = OrderedDict{Symbol, Float64}()
        data = getdata(fp, ith)
        dict[:min_wavelength] = minimum(coords(getdomain(fp ,ith)))
        dict[:max_wavelength] = maximum(coords(getdomain(fp ,ith)))
        dict[:min]     = minimum( values(data))
        dict[:max]     = maximum( values(data))
        dict[:mean]    = mean(    values(data))
        dict[:stdev]   = std(     values(data))
        dict[:median]  = median(  values(data))
        dict[:nmad]    = mad(     values(data))
        dict[:q0p01]   = quantile(values(data), 0.01)
        dict[:q0p05]   = quantile(values(data), 0.05)
        dict[:q0p16]   = quantile(values(data), 0.16)
        dict[:q0p84]   = quantile(values(data), 0.84)
        dict[:q0p95]   = quantile(values(data), 0.95)
        dict[:q0p99]   = quantile(values(data), 0.99)
        dict[:min_unc]    = minimum(uncerts(data))
        dict[:max_unc]    = maximum(uncerts(data))
        dict[:mean_unc]   = mean(   uncerts(data))
        dict[:median_unc] = median( uncerts(data))
        dict[:SNR] = median(abs.(values(data) ./ uncerts(data)))
        out[:Data_stats] = dict

        # Issues
        dict = OrderedDict{Symbol, String}()
        for cname in reverse(getfield.(GModelFit.flatten(GModelFit.deptree(model)), :cname))
            comp = model[cname]
            comp_issues = String[]
            for (pname, par) in GModelFit.getparams(comp)
                par_issues = String[]
                if !par.fixed
                    isnan(par.unc)                                                   &&  push!(par_issues, "uncertainty is NaN")
                    (par.val in [par.low, par.high])                                 &&  push!(par_issues, "value hit a limit")
                    if (par.low < 0)  &&  (par.high > 0)
                        ((par.unc / abs(par.val)) > recipe.qflag_relunc_threshold)   &&  push!(par_issues, "relative uncertainty exceeds threshold")
                    end
                end
                if length(par_issues) > 0
                    push!(comp_issues, string(pname) * " (" * join(par_issues, ", ") * ")")
                end
            end

            dep_issues = String[]
            for dep in GModelFit.dependencies(model, cname)
                if dep in keys(dict)
                    push!(dep_issues, string(dep))
                end
            end
            if length(dep_issues) > 0
                push!(comp_issues, "Issue in depedencies: " * join(dep_issues, ", "))
            end

            if length(comp_issues) > 0
                dict[cname] = join(comp_issues, ", ")
            end
        end
        kk = reverse(collect(keys(dict)))
        out[:Issues] = OrderedDict([k => dict[k] for k in kk])

        # Continuum
        dict = OrderedDict{Symbol, Float64}()
        if :QSOcont in keys(out[:Issues])
            dict[:l1450] = NaN
            dict[:l3000] = NaN
            dict[:l5100] = NaN
        else
            comp = deepcopy(model[:QSOcont])
            dict[:l1450] = comp(Domain([1450.]))[1]
            dict[:l3000] = comp(Domain([3000.]))[1]
            dict[:l5100] = comp(Domain([5100.]))[1]
        end
        out[:Continuum_luminosity] = dict

        # Equivalent widths of emission lines
        cont = deepcopy(geteval(fp, ith, :Continuum))
        (:Iron in keys(model))           &&  (cont .+= geteval(fp, ith, :Iron))
        (:NuisanceLines in keys(model))  &&  (cont .+= geteval(fp, ith, :NuisanceLines))
        @assert all(cont .> 0) "Continuum model is zero or negative"

        dict = OrderedDict{Symbol, Float64}()
        for cname in keys(recipe.specs[ith].ctx[:lines])
            haskey(model, cname) || continue
            if cname in keys(out[:Issues])
                dict[Symbol(cname)] = NaN
            else
                dict[Symbol(cname)] = QSFit.int_tabulated(coords(getdomain(fp ,ith)),
                                                          geteval(fp, ith, cname) ./ cont)[1]
            end
        end
        out[:Equivalent_widths] = dict

        dict = OrderedDict{Symbol, Any}()
        if  (:Ha_br in keys(model))  &&
            (:Ha_bb in keys(model))
            if  !(:Ha_br in keys(out[:Issues]))  &&
                !(:Ha_bb in keys(out[:Issues]))
                dict[:Ha_ncomp] = 2
                dict[:Ha_comps] = "Ha_br Ha_bb"
                dict[:Ha_norm] = model[:Ha_br].norm.val + model[:Ha_bb].norm.val
                fwhm, voff = QSFit.estimate_fwhm_voff(coords(getdomain(fp ,ith)),
                                                      geteval(fp, ith, :Ha_br) +
                                                      geteval(fp, ith, :Ha_bb),
                                                      model[:Ha_br].center.val)
                dict[:Ha_fwhm] = fwhm * 3e5
                dict[:Ha_voff] = voff * 3e5
            end
        end
        out[:Line_associations] = dict

        push!(outm, out)
    end

    return outm
end
