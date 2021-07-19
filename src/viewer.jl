using GFitViewer
import GFitViewer: ViewerData, viewer

function ViewerData(res::QSFitResults{T}; kw...) where T
    vd = ViewerData(res.model, res.pspec.data, res.fitres; kw...)
    vd.dict[:meta][:banner] = "QSFit (v0.1)<br />Date: " * string(trunc(res.fitres.timestamp, Second))

    id = 1
    m = vd.dict[:models][id][:meta]
    m[:label] = res.source.name * ", z=" * string(res.source.z) * ", E(B-V)=" * string(res.source.mw_ebv)
    m[:label_x] = "Rest frame wavelength"
    m[:unit_x]  = string(QSFit.unit_λ())
    m[:log10scale_x] = 0
    m[:label_y] = "Lum. density"

    if "UNITFUL_FANCY_EXPONENTS" in keys(ENV)
        orig = ENV["UNITFUL_FANCY_EXPONENTS"]
        ENV["UNITFUL_FANCY_EXPONENTS"] = true
        m[:unit_y]  = string(QSFit.unit_lum_density())
        ENV["UNITFUL_FANCY_EXPONENTS"] = orig
    else
        ENV["UNITFUL_FANCY_EXPONENTS"] = true
        m[:unit_y]  = string(QSFit.unit_lum_density())
        delete!(ENV, "UNITFUL_FANCY_EXPONENTS")
    end
    m[:log10scale_y] = QSFit.log10_scale_lum()
    vd.dict[:data][id][:meta][:label] = res.source.specs[id].label

    m = vd.dict[:extra][id]
    m[:EW] = GFitViewer.MDict()
    m[:EW][:label] = "Em. lines EW"
    m[:EW][:fields] = GFitViewer.MDict()
    m[:EW][:fields][:Label] = GFitViewer.MDict()
    m[:EW][:fields][:Label][:meta] = GFitViewer.MDict()
    m[:EW][:fields][:Label][:meta][:desc] = "Em. line"
    m[:EW][:fields][:Label][:data] = collect(keys(res.EW))
    m[:EW][:fields][:Value] = GFitViewer.MDict()
    m[:EW][:fields][:Value][:meta] = GFitViewer.MDict()
    m[:EW][:fields][:Label][:meta][:desc] = "EW [A]"
    m[:EW][:fields][:Value][:data] = collect(values(res.EW))

    m[:extratab_2] = GFitViewer.MDict()
    m[:extratab_2][:label] = "Second extra table"
    m[:extratab_2][:fields] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_1] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_1][:meta] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_1][:meta][:desc] = "Optional. In case we want to add metedata."
    m[:extratab_2][:fields][:fname_1][:data] = [102, 203, 304]
    m[:extratab_2][:fields][:fname_2] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_2][:meta] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_2][:data] = [10.2, 20.3, 30.4]
    m[:extratab_2][:fields][:fname_3] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_3][:meta] = GFitViewer.MDict()
    m[:extratab_2][:fields][:fname_3][:data] = ["String_1", "String_2", "String_3"]
    return vd
end

viewer(res::QSFitResults{T}; filename=nothing, offline=false, kw...) where T =
    viewer(ViewerData(res; kw...); filename=filename, offline=offline)


function ViewerData(res::QSFitMultiResults{T}; kw...) where T
    vd = ViewerData(res.multi, res.source.data, res.fitres; kw...)
    vd.dict[:meta][:banner] = "QSFit (v0.1)<br />Date: " * string(trunc(res.fitres.timestamp, Second))

    for id in 1:length(res.multi)
        m = vd.dict[:models][id][:meta]
        m[:label] = res.source.name * ", z=" * string(res.source.z) * ", E(B-V)=" * string(res.source.mw_ebv)
        m[:label_x] = "Rest frame wavelength"
        m[:unit_x]  = string(QSFit.unit_λ())
        m[:log10scale_x] = 0
        m[:label_y] = "Lum. density"

        if "UNITFUL_FANCY_EXPONENTS" in keys(ENV)
            orig = ENV["UNITFUL_FANCY_EXPONENTS"]
            ENV["UNITFUL_FANCY_EXPONENTS"] = true
            m[:unit_y]  = string(QSFit.unit_lum_density())
            ENV["UNITFUL_FANCY_EXPONENTS"] = orig
        else
            ENV["UNITFUL_FANCY_EXPONENTS"] = true
            m[:unit_y]  = string(QSFit.unit_lum_density())
            delete!(ENV, "UNITFUL_FANCY_EXPONENTS")
        end
        m[:log10scale_y] = QSFit.log10_scale_lum()
        vd.dict[:data][id][:meta][:label] = res.source.spectra[id].label

        m = vd.dict[:extra][id]
        m[:EW] = GFitViewer.MDict()
        m[:EW][:label] = "Em. lines EW"
        m[:EW][:fields] = GFitViewer.MDict()
        m[:EW][:fields][:Label] = GFitViewer.MDict()
        m[:EW][:fields][:Label][:meta] = GFitViewer.MDict()
        m[:EW][:fields][:Label][:meta][:desc] = "Em. line"
        m[:EW][:fields][:Label][:data] = collect(keys(res.EW[id]))
        m[:EW][:fields][:Value] = GFitViewer.MDict()
        m[:EW][:fields][:Value][:meta] = GFitViewer.MDict()
        m[:EW][:fields][:Label][:meta][:desc] = "EW [A]"
        m[:EW][:fields][:Value][:data] = collect(values(res.EW[id]))

        m[:extratab_2] = GFitViewer.MDict()
        m[:extratab_2][:label] = "Second extra table"
        m[:extratab_2][:fields] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_1] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_1][:meta] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_1][:meta][:desc] = "Optional. In case we want to add metedata."
        m[:extratab_2][:fields][:fname_1][:data] = [102, 203, 304]
        m[:extratab_2][:fields][:fname_2] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_2][:meta] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_2][:data] = [10.2, 20.3, 30.4]
        m[:extratab_2][:fields][:fname_3] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_3][:meta] = GFitViewer.MDict()
        m[:extratab_2][:fields][:fname_3][:data] = ["String_1", "String_2", "String_3"]
    end
    return vd
end

viewer(res::QSFitMultiResults{T}; filename=nothing, offline=false, kw...) where T =
    viewer(ViewerData(res; kw...); filename=filename, offline=offline)
