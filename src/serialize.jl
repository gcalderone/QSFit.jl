import GModelFit: _serialize

_serialize(vv::Quantity) = string(vv)
_serialize(vv::Results) = GModelFit._serialize_struct(vv)
_serialize(vv::Spectrum) = GModelFit._serialize_struct(vv)
