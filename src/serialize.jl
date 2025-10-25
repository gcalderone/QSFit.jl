using JSON, GZip

function serialize(filename, res::Results; compress=false)
    v = GModelFit._serialize(res.bestfit, res.fsumm, res.data)
    push!(v, OrderedDict{String, Any}("elapsed" => res.elapsed, "timestamp" => GModelFit._serialize(res.timestamp)))
    push!(v, GModelFit._serialize(res.spec))
    push!(v, GModelFit._serialize(res.post))

    filename = GModelFit.ensure_file_extension(filename, "json")

    filename = ensure_file_extension(filename, "json")
    if compress
        filename = ensure_file_extension(filename, "gz")
        io = GZip.open(filename, "w")
    else
        io = open(filename, "w")
    end
    JSON.print(io, v)
    close(io)
    return filename
end


function deserialize(filename::String)
    if filename[end-2:end] == ".gz"
        io = GZip.open(filename)
    else
        io = open(filename)
    end
    v = JSON.parse(io, dicttype=OrderedDict)
    close(io)

    return Results(GModelFit._deserialize(v[4]["timestamp"]), v[4]["elapsed"],
                   GModelFit._deserialize(v[5]),
                   GModelFit._deserialize(v[3]),
                   GModelFit._deserialize(v[1]),
                   GModelFit._deserialize(v[2]),
                   GModelFit._deserialize(v[6]))
end
