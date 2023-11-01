
module InternalPrinting

using DocStringExtensions

# copied from COBREXA v1.5.1

"""
$(TYPEDSIGNATURES)

Nicely prints keys and values.
"""
_pretty_print_keyvals(io, def::String, payload; kwargs...) =
    _pretty_print_keyvals(io, def, isnothing(payload) ? "---" : string(payload); kwargs...)

"""
$(TYPEDSIGNATURES)

Specialization of `_pretty_print_keyvals` for plain strings.
"""
function _pretty_print_keyvals(io, def::String, payload::String)
    print(io, def)
    if isempty(payload)
        println(io, "---")
    else
        println(io, payload)
    end
end

"""
$(TYPEDSIGNATURES)

Specialization of `_pretty_print_keyvals` for dictionaries.
"""
function _pretty_print_keyvals(io, def::String, payload::Dict)

    print(io, def)
    if isempty(payload)
        println(io, "---")
    else
        println(io, "")
        for (k, v) in payload
            if length(v) > 2 && length(v[1]) < 20
                println(io, "\t", k, ": ", v[1], ", ..., ", v[end])
            elseif length(v[1]) > 20 # basically for envipath annotations... or long notes
                println(io, "\t", k, ": ", v[1][1:20], "...")
            else
                println(io, "\t", k, ": ", v)
            end
        end
    end
end

"""
$(TYPEDSIGNATURES)

Print out long reaction stoichiometries nicely.
"""
function _pretty_substances(ss::Vector{String})::String
    if isempty(ss)
        "∅"
    elseif length(ss) > 5
        join([ss[1], ss[2], "...", ss[end-1], ss[end]], " + ")
    else
        join(ss, " + ")
    end
end

end
