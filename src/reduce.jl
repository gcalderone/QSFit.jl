
function estimate_line_EWs(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: AbstractRecipe
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(model())
    for (lname, lc) in pspec.lcs
        haskey(model, lname) || continue
        cont .-= model(lname)
    end
    for (lname, lc) in pspec.lcs
        haskey(model, lname) || continue
        EW[lname] = int_tabulated(domain(model)[:],
                                  model(lname) ./ cont)[1]
    end
    return EW
end


struct QSFitResults{T}
    source::QSO{T}
    pspec::PreparedSpectrum
    model::Model
    bestfit::GFit.BestFitResult
    EW::OrderedDict{Symbol, Float64}
end

QSFitResults(source::QSO{T}, pspec::PreparedSpectrum, model::Model, bestfit::GFit.BestFitResult) where T <: AbstractRecipe =
    QSFitResults{T}(source, pspec, model, bestfit, estimate_line_EWs(source, pspec, model))

struct QSFitMultiResults{T}
    source::QSO{T}
    multi::MultiModel
    bestfit::GFit.BestFitMultiResult
    EW::Vector{OrderedDict{Symbol, Float64}}
end

QSFitMultiResults(source::QSO{T}, multi::MultiModel, bestfit::GFit.BestFitMultiResult) where T <: AbstractRecipe =
    QSFitMultiResults{T}(source, multi, bestfit, [estimate_line_EWs(source, multi[id]) for id in 1:length(multi)])
