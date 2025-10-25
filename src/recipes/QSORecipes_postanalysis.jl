function measure_line_profile(x, y, center, comps::Vector{T}) where T <: QSFit.AbstractSpecLineComp
    fwhm, voff = QSFit.estimate_fwhm_voff(x, y, center)
    fwhm *= 3e5
    voff *= 3e5
    norms = [comp.norm.val for comp in comps]
    weights = norms ./ sum(norms)
    out = OrderedDict{Symbol, Any}()
    out[:ncomp]    = length(norms)
    out[:norm]     = sum(norms)
    out[:norm_unc] = sum(weights .* [comp.norm.unc                 for comp in comps])
    out[:fwhm]     = fwhm
    out[:fwhm_unc] = sum(weights .* [comp.fwhm.unc / comp.fwhm.val for comp in comps]) * fwhm
    out[:voff]     = voff
    out[:fwhm_unc] = sum(weights .* [comp.fwhm.unc                 for comp in comps])
    return out
end

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
                    !isnothing(par.patch)                                                  &&  push!(par_issues, "patched value")
                    isnan(par.unc)                                                         &&  push!(par_issues, "uncertainty is NaN")
                    (par.val in [par.low, par.high])                                       &&  push!(par_issues, "value hit a limit")
                    if (par.low < 0)  &&  (par.high > 0)
                        ((par.unc / abs(par.val)) > recipe.reliability_relunc_threshold)   &&  push!(par_issues, "relative uncertainty exceeds threshold")
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
        λ = coords(getdomain(fp ,ith))
        comp = deepcopy(model[:QSOcont])
        for wl in [1450, 3000, 5100]
            dict[Symbol(:l, wl)] = NaN
            if  !(:QSOcont in keys(out[:Issues]))  &&
                (minimum(λ) <= wl <= maximum(λ))
                dict[Symbol(:l, wl)] = comp(Domain([float(wl)]))[1] * wl # lL_l [10^42 erg/s]
            end
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

        # Line associations
        nuisance_cnames = (:NuisanceLines in keys(model)  ?  GModelFit.dependencies(model, :NuisanceLines)  :  Symbol[])
        for cname in keys(recipe.check_for_line_assoc)
            if cname in keys(model)
                assoc = [cname]
                assoc_lines_with_issues = Symbol[]

                for d in recipe.check_for_line_assoc[cname]
                    if d in keys(model)
                        # Avoid checking if d in keys(out[:Issues])
                        # since that component is typically patched to
                        # the main one and one or more of its
                        # parameters may rise issues
                        println("Line association: $cname - $d")
                        push!(assoc, d)
                    end
                end

                center = model[cname].center.val * (1 - (model[cname].voff.val / 3e5))
                i = 1
                while i <= length(nuisance_cnames)
                    d = nuisance_cnames[i]
                    limit = model[d].fwhm.val * model[d].center.val / 3e5 / 2
                    if abs(center - model[d].center.val) < limit
                        println("Line association: $cname - $d")
                        deleteat!(nuisance_cnames, i)
                        if d in keys(out[:Issues])
                            push!(assoc_lines_with_issues, d)
                        else
                            push!(assoc, d)
                        end
                        continue
                    end
                    i += 1
                end

                if length(assoc_lines_with_issues) > 0
                    if cname in keys(out[:Issues])
                        out[:Issues][cname] *= ", "
                    else
                        out[:Issues][cname]  = ""
                    end
                    out[:Issues][cname] *= "Associated lines with issues (" * join(string.(assoc_lines_with_issues), ", ") * ")"
                else
                    if length(assoc) > 1
                        out[Symbol(cname, :_assoc)] = measure_line_profile(coords(getdomain(fp ,ith)),
                                                                           sum([geteval(fp, ith, assoc[i]) for i in 1:length(assoc)]),
                                                                           model[cname].center.val,
                                                                           [model[k] for k in assoc])
                        out[Symbol(cname, :_assoc)][:comps] = join(string.(assoc), ", ")
                    end
                end
            end
        end

        push!(outm, out)
    end

    return outm
end
