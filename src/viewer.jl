using GModelFitViewer
import GModelFitViewer: ViewerData

function ViewerData_meta(res::Results; kws...)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    return GModelFitViewer.Meta(; title=res.spec[:label],
                                xlabel="Wavelength",
                                xunit=res.spec[:xunit],
                                xscale=res.spec[:xscale],
                                ylabel="Lum. density",
                                yunit=res.spec[:yunit],
                                yscale=res.spec[:yscale],
                                kws...)
end


function ViewerData_extra(res::Results)
    out = OrderedDict{Symbol, Any}()
    for (majorkey, dict) in res.post
        out[majorkey] = OrderedDict{Symbol, Any}()
        out[majorkey][:label] = replace(string(majorkey), "_" => " ")
        out[majorkey][:fields] = OrderedDict{Symbol, Any}()
        out[majorkey][:fields][:Key] = OrderedDict{Symbol, Any}()
        out[majorkey][:fields][:Key][:meta] = OrderedDict{Symbol, Any}()
        out[majorkey][:fields][:Key][:meta][:desc] = "Key"
        out[majorkey][:fields][:Key][:data] = collect(keys(dict))

        out[majorkey][:fields][:Value] = OrderedDict{Symbol, Any}()
        out[majorkey][:fields][:Value][:meta] = OrderedDict{Symbol, Any}()
        out[majorkey][:fields][:Value][:meta][:desc] = "Value"
        vv = collect(values(dict))
        if isa(vv, Vector)
            out[majorkey][:fields][:Value][:data] = vv
        elseif isa(vv, Vector{NTuple{2, Float64}})
            out[majorkey][:fields][:Value][:data] = getindex.(vv, 1)
            out[majorkey][:fields][:Uncert] = OrderedDict{Symbol, Any}()
            out[majorkey][:fields][:Uncert][:meta] = OrderedDict{Symbol, Any}()
            out[majorkey][:fields][:Uncert][:meta][:desc] = "Uncertainty"
            out[majorkey][:fields][:Uncert][:data] = getindex.(vv, 2)
        end
    end
    return [out]
end


function ViewerData(res::Results; kws...)
    meta = ViewerData_meta(res; kws...)
    out = GModelFitViewer.ViewerData(res.bestfit, res.fsumm, res.data, meta=meta)
    out.data[1]["extra"] = ViewerData_extra(res)
    return out
end


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
