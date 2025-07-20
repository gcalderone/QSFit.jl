function postanalysis(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe

    outm = Vector{OrderedDict{Symbol, Any}}()

    for ith in 1:nmodels(fp)
        out = OrderedDict{Symbol, Any}()
        cnames = Symbol[]
        for i in 1:length(recipe.specs)
            append!(cnames, keys(recipe.specs[i].meta[:lines]))
        end
        cnames = sort(unique(cnames))
        model = getmodel(fp, ith)

        # Continuum
        comp = deepcopy(model[:QSOcont])
        out[:Cont_1450] = comp(Domain([1450.]))[1]
        out[:Cont_3000] = comp(Domain([3000.]))[1]
        out[:Cont_5100] = comp(Domain([5100.]))[1]

        cont = deepcopy(geteval(fp, ith, :Continuum))
        (:Iron in keys(model))           &&  (cont .+= geteval(fp, ith, :Iron))
        (:NuisanceLines in keys(model))  &&  (cont .+= geteval(fp, ith, :NuisanceLines))
        @assert all(cont .> 0) "Continuum model is zero or negative"

        for cname in cnames
            haskey(model, cname) || continue
            out[Symbol(:EW_, cname)] = QSFit.int_tabulated(coords(getdomain(fp ,ith)),
                                                           geteval(fp, ith, cname) ./ cont)[1]
        end
        push!(outm, out)
    end

    return outm
end
