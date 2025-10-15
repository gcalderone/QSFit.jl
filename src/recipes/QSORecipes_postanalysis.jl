function postanalysis(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe

    outm = Vector{OrderedDict{Symbol, Any}}()

    for ith in 1:nmodels(fp)
        out = OrderedDict{Symbol, Any}()
        model = getmodel(fp, ith)

        # Continuum
        dict = OrderedDict{Symbol, Float64}()
        comp = deepcopy(model[:QSOcont])
        dict[:l1450] = comp(Domain([1450.]))[1]
        dict[:l3000] = comp(Domain([3000.]))[1]
        dict[:l5100] = comp(Domain([5100.]))[1]
        out[:Continuum_luminosity] = dict

        # Equivalent widths
        cont = deepcopy(geteval(fp, ith, :Continuum))
        (:Iron in keys(model))           &&  (cont .+= geteval(fp, ith, :Iron))
        (:NuisanceLines in keys(model))  &&  (cont .+= geteval(fp, ith, :NuisanceLines))
        @assert all(cont .> 0) "Continuum model is zero or negative"

        cnames = Symbol[]
        for i in 1:length(recipe.specs)
            append!(cnames, keys(recipe.specs[i].meta[:lines]))
        end
        cnames = sort(unique(cnames))

        dict = OrderedDict{Symbol, Float64}()
        for cname in cnames
            haskey(model, cname) || continue
            dict[Symbol(cname)] = QSFit.int_tabulated(coords(getdomain(fp ,ith)),
                                                      geteval(fp, ith, cname) ./ cont)[1]
        end
        out[:Equivalent_widths] = dict

        # Parameter quality flags
        dict = OrderedDict{Symbol, Int64}()
        for (cname, comp) in model
            for (pname, par) in GModelFit.getparams(comp)
                if !par.fixed
                    qflag = 0
                    isnan(par.unc)                                                &&  (qflag += 2^1)
                    ((par.unc / abs(par.val)) > recipe.qflag_relunc_threshold)    &&  (qflag += 2^2)
                    (par.val in [par.low, par.high])                              &&  (qflag += 2^3)
                    (qflag != 0)  &&  (dict[Symbol(cname, :_, pname)] = qflag)
                end
            end
        end
        out[:Quality_flags] = dict

        push!(outm, out)
    end

    return outm
end
