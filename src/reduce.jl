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
