using GModelFitViewer
import GModelFitViewer: Meta, gmfv_lower

function Meta(res::Results; kws...)
    # ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    return GModelFitViewer.Meta(; title=res.spec[:label],
                                xlabel="Wavelength",
                                xunit=res.spec[:xunit],
                                xscale=res.spec[:xscale],
                                ylabel="Lum. density",
                                yunit=res.spec[:yunit],
                                yscale=res.spec[:yscale],
                                kws...)
end


function gmfv_lower(res::Results; kws...)
    extra = OrderedDict{Symbol, Any}()
    for (tab, dict) in res.post
        @assert isa(dict, AbstractDict)

        if tab == :Data_stats
            extra[tab] = OrderedDict(:label => "Data stats", :fields => OrderedDict{Symbol, Any}())
            extra[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Quantity"),
                                                   :data => string.(collect(keys(dict))))
            extra[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "Value"),
                                                   :data => collect(values(dict)))
        elseif tab == :Issues
            extra[tab] = OrderedDict(:label => "Issues", :fields => OrderedDict{Symbol, Any}())
            extra[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Component"),
                                                   :data => Vector{String}())
            extra[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "Parameter / dependency"),
                                                   :data => Vector{String}())
            extra[tab][:fields][:c3] = OrderedDict(:meta => Dict(:desc => "Issue"),
                                                   :data => Vector{String}())
            for (cname, subdict) in dict
                for (par_or_dep, issue) in subdict
                    push!(extra[tab][:fields][:c1][:data], string(cname))
                    push!(extra[tab][:fields][:c2][:data], string(par_or_dep))
                    push!(extra[tab][:fields][:c3][:data], issue)
                end
            end
        elseif tab == :Continuum_luminosity
            extra[tab] = OrderedDict(:label => "Continuum luminosities", :fields => OrderedDict{Symbol, Any}())
            extra[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Wavelength [A]"),
                                                   :data => string.(collect(keys(dict))))
            extra[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "lambda L_lambda [10^42 erg/s]"),
                                                   :data => string.(collect(values(dict))))
        elseif tab == :Equivalent_widths
            extra[tab] = OrderedDict(:label => "Equivalent widths", :fields => OrderedDict{Symbol, Any}())
            extra[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Emission line"),
                                                   :data => string.(collect(keys(dict))))
            extra[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "[A]"),
                                                   :data => string.(collect(values(dict))))
        elseif isa(dict, OrderedDict)
            vv = string.(collect(values(dict)))
            extra[tab] = OrderedDict(:label => tab, :fields => OrderedDict{Symbol, Any}())
            extra[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Key"),
                                                   :data => string.(collect(keys(dict))))
            extra[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "Value"),
                                                   :data => string.(collect(values(dict))))
        else
            @warn "Unexpected entry in post-analysis dictionary: $tab"
        end
    end
    out = GModelFitViewer.gmfv_lower(res.bestfit, res.fsumm, res.data, Meta(res; kws...))
    out.dict[:extra] = TypedJSON.lower(extra)
    return out
end


#=
function ViewerData(multires::MultiResults; kws...)
    N = length(multires.bestfit)
    res = [Results(multires.timestamp,
                   multires.elapsed,
                   multires.spec[i],
                   multires.data[i],
                   multires.bestfit[i],
                   multires.fsumm::GModelFit.FitSummary,
                   multires.post[i]) for i in 1:N]

    meta = [ViewerData_meta(res[i]; kws...) for i in 1:N]
    out = GModelFitViewer.ViewerData(multires.bestfit, multires.fsumm, multires.data, meta=meta)
    for i in 1:N
        out.data[1][i]["extra"] = ViewerData_extra(res[i])
    end
    return out
end
=#
