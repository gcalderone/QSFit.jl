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
    if haskey(res.post, :EW)
        m = OrderedDict{Symbol, Any}()
        m[:EW] = OrderedDict{Symbol, Any}()
        m[:EW][:label] = "Em. lines EW"
        m[:EW][:fields] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Label] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Label][:meta] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Label][:meta][:desc] = "Em. line"
        m[:EW][:fields][:Label][:data] = collect(keys(res.post[:EW]))
        m[:EW][:fields][:Value] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Value][:meta] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Value][:meta][:desc] = "EW [A]"
        m[:EW][:fields][:Value][:data] = collect(values(res.post[:EW]))
        d.data[1]["extra"] = [m]
    end

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
