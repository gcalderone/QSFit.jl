using GModelFitViewer
import GModelFitViewer: ViewerData

function ViewerData(res::Results; kws...)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    i = findall(isnothing.(match.(r"SpecLine", ctypes))  .&   isnothing.(match.(r"GaussConv", ctypes)))
    keep = string.(keys(res.bestfit))[i]
    meta = GModelFitViewer.Meta(; title=res.spec.label,
                                xlabel="Wavelength",
                                xunit=string(unit(res.spec.unit_x)),
                                xscale=ustrip(res.spec.unit_x),
                                ylabel="Lum. density",
                                yunit=string(unit(res.spec.unit_y)),
                                yscale=ustrip(res.spec.unit_y),
                                keep=keep, kws...)

    d = GModelFitViewer.ViewerData(res.bestfit, res.fsumm, res.data, meta=meta)

    # Add further content
    extra = Vector{OrderedDict{Symbol, Any}}()

    if length(res.post) > 0
        m = OrderedDict{Symbol, Any}()
        m[:post] = OrderedDict{Symbol, Any}()
        m[:post][:label] = "Post-analysis"
        m[:post][:fields] = OrderedDict{Symbol, Any}()
        m[:post][:fields][:Label] = OrderedDict{Symbol, Any}()
        m[:post][:fields][:Label][:meta] = OrderedDict{Symbol, Any}()
        m[:post][:fields][:Label][:meta][:desc] = "Key"
        m[:post][:fields][:Label][:data] = collect(keys(res.post))
        m[:post][:fields][:Value] = OrderedDict{Symbol, Any}()
        m[:post][:fields][:Value][:meta] = OrderedDict{Symbol, Any}()
        m[:post][:fields][:Value][:meta][:desc] = "Value"
        m[:post][:fields][:Value][:data] = collect(values(res.post))
        push!(extra, m)
    end

    d.data[1]["extra"] = extra
    return d
end



function ViewerData(res::MultiResults; kws...)
    meta = Vector{GModelFitViewer.Meta}()
    for ith in 1:length(res.bestfit)
        ctypes = [comptype(res.bestfit[ith], cname) for cname in keys(res.bestfit[ith])]
        i = findall(isnothing.(match.(r"SpecLine", ctypes))  .&   isnothing.(match.(r"GaussConv", ctypes)))
        keep = string.(keys(res.bestfit[ith]))[i]
        push!(meta, GModelFitViewer.Meta(; title=res.spec[ith].label,
                                         xlabel="Wavelength",
                                         xunit=string(unit(res.spec[ith].unit_x)),
                                         xscale=ustrip(res.spec[ith].unit_x),
                                         ylabel="Lum. density",
                                         yunit=string(unit(res.spec[ith].unit_y)),
                                         yscale=ustrip(res.spec[ith].unit_y),
                                         keep=keep, kws...))
    end

    d = GModelFitViewer.ViewerData(res.bestfit, res.fsumm, res.data, meta=meta)
    return d
end
