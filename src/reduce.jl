
struct QSFitResults{T}
    source::QSO{T}
    model::Model
    bestfit::GFit.BestFitResult
    EW::Vector{OrderedDict{Symbol, Float64}}
end


function reduce(source::QSO{T}, model::Model, bestfit::GFit.BestFitResult) where T <: AbstractRecipe
    # Estimate emission line EWs
    EW = Vector{OrderedDict{Symbol, Float64}}()

    for id in 1:length(model.preds)
        push!(EW, OrderedDict{Symbol, Float64}())
        
        for (lname0, line0) in known_spectral_lines(T)
            cont = deepcopy(model[id]())
            for (lname, line) in line_breakdown(T, lname0, line0)
                haskey(model[id], lname) || continue
                cont .-= model[id](lname)
            end
            
            for (lname, line) in line_breakdown(T, lname0, line0)
                haskey(model[id], lname) || continue
                EW[id][lname] = int_tabulated(domain(model[id])[:],
                                              model[id](lname) ./ cont)[1]
            end
        end
    end
    return QSFitResults{T}(source, model, bestfit, EW)
end
