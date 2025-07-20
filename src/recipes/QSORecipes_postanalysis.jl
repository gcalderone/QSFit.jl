function postanalysis(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe

    out = Vector{OrderedDict{Symbol, Any}}()

    for ith in 1:nmodels(fp)
        cnames = Symbol[]
        for i in 1:length(recipe.specs)
            append!(cnames, keys(recipe.specs[i].meta[:lines]))
        end
        cnames = sort(unique(cnames))

        model = getmodel(fp, ith)
        cont = deepcopy(geteval(fp, ith, :Continuum))
        (:Iron in keys(model))           &&  (cont .+= geteval(fp, ith, :Iron))
        (:NuisanceLines in keys(model))  &&  (cont .+= geteval(fp, ith, :NuisanceLines))
        @assert all(cont .> 0) "Continuum model is zero or negative"

        EW = OrderedDict{Symbol, Float64}()
        for cname in cnames
            haskey(model, cname) || continue
            EW[cname] = QSFit.int_tabulated(coords(getdomain(fp ,ith)),
                                            geteval(fp, ith, cname) ./ cont)[1]
        end
        push!(out, OrderedDict{Symbol, Any}(:EW => EW))
    end

    return out
end
