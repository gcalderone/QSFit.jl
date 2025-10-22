using JSON

function serialize(filename, res::Results)
    v = GModelFit._serialize(res.bestfit, res.fsumm, res.data)
    push!(v, OrderedDict{String, Any}("elapsed" => res.elapsed, "timestamp" => GModelFit._serialize(res.timestamp)))
    push!(v, GModelFit._serialize(res.spec))
    push!(v, GModelFit._serialize(res.post))

    filename = GModelFit.ensure_file_extension(filename, "json")
    io = open(filename, "w")
    JSON.print(io, v)
    close(io)
    return filename
end


function deserialize(filename::String)
    io = open(filename)
    v = JSON.parse(io, dicttype=OrderedDict)
    close(io)

    return Results(GModelFit._deserialize(v[4]["timestamp"]), v[4]["elapsed"],
                   GModelFit._deserialize(v[5]),
                   GModelFit._deserialize(v[3]),
                   GModelFit._deserialize(v[1]),
                   GModelFit._deserialize(v[2]),
                   GModelFit._deserialize(v[6]))
end
