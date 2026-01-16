using TypedJSON, DataStructures
import TypedJSON: lower, reconstruct

TypedJSON.reconstruct(::Val{Symbol("QSFit.Results")}, dict) = Results(values(dict)...)

