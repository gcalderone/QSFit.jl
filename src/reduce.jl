
function estimate_line_EWs(source::QSO{T}, model::Model) where T <: AbstractRecipe
    EW = OrderedDict{Symbol, Float64}()

    for (lname0, line0) in known_spectral_lines(T)
        cont = deepcopy(model())
        for (lname, line) in line_breakdown(T, lname0, line0)
            haskey(model, lname) || continue
            cont .-= model(lname)
        end

        for (lname, line) in line_breakdown(T, lname0, line0)
            haskey(model, lname) || continue
            EW[lname] = int_tabulated(domain(model)[:],
                                      model(lname) ./ cont)[1]
        end
    end
    return EW
end


struct QSFitResults{T}
    source::QSO{T}
    model::Model
    bestfit::GFit.BestFitResult
    EW::OrderedDict{Symbol, Float64}
end

QSFitResults(source::QSO{T}, model::Model, bestfit::GFit.BestFitResult) where T <: AbstractRecipe =
    QSFitResults{T}(source, model, bestfit, estimate_line_EWs(source, model))

struct QSFitMultiResults{T}
    source::QSO{T}
    multi::MultiModel
    bestfit::GFit.BestFitMultiResult
    EW::Vector{OrderedDict{Symbol, Float64}}
end

QSFitMultiResults(source::QSO{T}, multi::MultiModel, bestfit::GFit.BestFitMultiResult) where T <: AbstractRecipe =
    QSFitMultiResults{T}(source, multi, bestfit, [estimate_line_EWs(source, multi[id]) for id in 1:length(multi)])
