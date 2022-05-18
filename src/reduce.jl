function estimate_line_EWs(job::JobState{T}) where T <: AbstractRecipe
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(job.model())
    for (lname, lc) in job.pspec.lcs
        haskey(job.model, lname) || continue
        cont .-= job.model(lname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for (lname, lc) in job.pspec.lcs
        haskey(job.model, lname) || continue
        EW[lname] = int_tabulated(domain(job.model)[:],
                                  job.model(lname) ./ cont)[1]
    end
    return EW
end


struct QSFitMultiResults{T}
    source::QSO{T}
    pspecs::Vector{PreparedSpectrum}
    multi::MultiModel
    fitres::GFit.FitResult
    EW::Vector{OrderedDict{Symbol, Float64}}
end

QSFitMultiResults(source::QSO{T}, pspecs::Vector{PreparedSpectrum}, multi::MultiModel, fitres::GFit.FitResult) where T <: AbstractRecipe =
    QSFitMultiResults{T}(source, pspecs, multi, fitres, [estimate_line_EWs(source, pspecs[id], multi[id]) for id in 1:length(multi)])


function qsfit(res::JobResults{T}) where T
    fitres = fit!(res.source, res.model, res.pspec)
    return JobResults(res.source, res.pspec, res.model, fitres)
end

function qsfit(res::QSFitMultiResults{T}) where T
    fitres = fit!(res.source, res.model, res.pspecs)
    return QSFitMultiResults(res.source, res.pspecs, res.model, fitres)
end
