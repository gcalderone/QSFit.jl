export QSFitResults, QSFitMultiResults

function estimate_line_EWs(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: AbstractRecipe
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(model())
    for (lname, lc) in pspec.lcs
        haskey(model, lname) || continue
        cont .-= model(lname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
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
    fitres::GFit.FitResult
    EW::OrderedDict{Symbol, Float64}
end

QSFitResults(source::QSO{T}, pspec::PreparedSpectrum, model::Model, fitres::GFit.FitResult) where T <: AbstractRecipe =
    QSFitResults{T}(source, pspec, model, fitres, estimate_line_EWs(source, pspec, model))

struct QSFitMultiResults{T}
    source::QSO{T}
    pspecs::Vector{PreparedSpectrum}
    multi::MultiModel
    fitres::GFit.FitResult
    EW::Vector{OrderedDict{Symbol, Float64}}
end

QSFitMultiResults(source::QSO{T}, pspecs::Vector{PreparedSpectrum}, multi::MultiModel, fitres::GFit.FitResult) where T <: AbstractRecipe =
    QSFitMultiResults{T}(source, pspecs, multi, fitres, [estimate_line_EWs(source, pspecs[id], multi[id]) for id in 1:length(multi)])


function qsfit(res::QSFitResults{T}) where T
    fitres = fit!(res.source, res.model, res.pspec)
    return QSFitResults(res.source, res.pspec, res.model, fitres)
end

function qsfit(res::QSFitMultiResults{T}) where T
    fitres = fit!(res.source, res.model, res.pspecs)
    return QSFitMultiResults(res.source, res.pspecs, res.model, fitres)
end
