function reduce(recipe::RRef{T}, state::State) where T <: AbstractRecipe
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(state.model())
    for (lname, lc) in state.pspec.lcs
        haskey(state.model, lname) || continue
        cont .-= state.model(lname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for (lname, lc) in state.pspec.lcs
        haskey(state.model, lname) || continue
        EW[lname] = int_tabulated(coords(domain(state.model)),
                                  state.model(lname) ./ cont)[1]
    end
    return OrderedDict{Symbol, Any}(:EW => EW)
end
