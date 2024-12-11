import GModelFit: _serialize, _deserialize

function _serialize(vv::Quantity)
    @warn "Units can not be serialized: deserialized ones may be inconsistent!"
    string(vv)
end

_serialize(vv::Results) = GModelFit._serialize_struct(vv)
_serialize(vv::Spectrum) = GModelFit._serialize_struct(vv)

_deserialize(::Val{Symbol("QSFit.Results")},
             dd::AbstractDict) =
                 Results(_deserialize(dd["timestamp"]),
                         _deserialize(dd["elapsed"]),
                         _deserialize(dd["spec"]),
                         _deserialize(dd["data"]),
                         _deserialize(dd["bestfit"]),
                         _deserialize(dd["fitstats"]),
                         _deserialize(dd["reduced"]))

function _deserialize(::Val{Symbol("QSFit.Spectrum")}, dd::AbstractDict)
    out = Spectrum(_deserialize(dd["x"]),
                   _deserialize(dd["y"]),
                   _deserialize(dd["err"]),
                   label=_deserialize(dd["label"]),
                   good=convert(Vector{Bool}, _deserialize(dd["good"])),
                   resolution=_deserialize(dd["resolution"]))
    out.isrestframe = _deserialize(dd["isrestframe"])
    out.meta = _deserialize(dd["meta"])
    if  (string(out.unit_x) != dd["unit_y"])  ||
        (string(out.unit_y) != dd["unit_y"])
        @warn "Can't deserialize spectrum units!"
        println("Serialized units were:")
        println("  x: ", dd["unit_x"])
        println("  y: ", dd["unit_y"])
        println("Deserialized units are:")
        println("  x: ", string(out.unit_x))
        println("  y: ", string(out.unit_y))
    end
    return out
end
    
