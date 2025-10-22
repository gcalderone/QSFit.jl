function postanalysis(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe

    outm = Vector{OrderedDict{Symbol, Any}}()

    for ith in 1:nmodels(fp)
        out = OrderedDict{Symbol, Any}()
        model = getmodel(fp, ith)

        # Continuum
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

        # Continuum
        dict = OrderedDict{Symbol, Float64}()
        comp = deepcopy(model[:QSOcont])
        dict[:l1450] = comp(Domain([1450.]))[1]
        dict[:l3000] = comp(Domain([3000.]))[1]
        dict[:l5100] = comp(Domain([5100.]))[1]
        out[:Continuum_luminosity] = dict

        # Equivalent widths of emission lines
        cont = deepcopy(geteval(fp, ith, :Continuum))
        (:Iron in keys(model))           &&  (cont .+= geteval(fp, ith, :Iron))
        (:NuisanceLines in keys(model))  &&  (cont .+= geteval(fp, ith, :NuisanceLines))
        @assert all(cont .> 0) "Continuum model is zero or negative"

        dict = OrderedDict{Symbol, Float64}()
        for cname in keys(recipe.specs[ith].meta[:lines])
            haskey(model, cname) || continue
            dict[Symbol(cname)] = QSFit.int_tabulated(coords(getdomain(fp ,ith)),
                                                      geteval(fp, ith, cname) ./ cont)[1]
        end
        out[:Equivalent_widths] = dict

        # Component and parameter quality flags
        dict = OrderedDict{Symbol, Int64}()
        for (cname, comp) in model
            dict[cname] = 0
            for (pname, par) in GModelFit.getparams(comp)
                dict[Symbol(cname, :_, pname)] = 0
                if !par.fixed
                    qflag = 0
                    isnan(par.unc)                                                   &&  (qflag += 2^0)
                    (par.val in [par.low, par.high])                                 &&  (qflag += 2^1)
                    if (par.low < 0)  &&  (par.high > 0)
                        ((par.unc / abs(par.val)) > recipe.qflag_relunc_threshold)   &&  (qflag += 2^2)
                    end
                    dict[Symbol(cname, :_, pname)] = qflag
                    (qflag != 0)  &&  (dict[cname] = 2^0)  # set also component flag
                end
            end
        end
        out[:Quality_flags] = dict

        push!(outm, out)
    end

    return outm
end
