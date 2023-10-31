import .CanonicalModel as C
import .InternalPrinting as IP

# copied from COBREXA v1.5.1

function Base.show(io::IO, ::MIME"text/plain", g::C.Gene)
    for fname in fieldnames(C.Gene)
        IP._pretty_print_keyvals(io, "Gene.$(string(fname)): ", getfield(g, fname))
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::C.Metabolite)
    for fname in fieldnames(C.Metabolite)
        if fname == :charge
            c = isnothing(getfield(m, fname)) ? nothing : string(getfield(m, fname))
            IP._pretty_print_keyvals(io, "Metabolite.$(string(fname)): ", c)
        elseif fname == :formula
            c = isnothing(getfield(m, fname)) ? nothing : join([k*string(v) for (k, v) in getfield(m, fname)])
            IP._pretty_print_keyvals(io, "Metabolite.$(string(fname)): ", c)    
        else
            IP._pretty_print_keyvals(io, "Metabolite.$(string(fname)): ", getfield(m, fname))
        end
    end
end

function _pretty_substances(ss::Vector{String})::String
    if isempty(ss)
        "∅"
    elseif length(ss) > 5
        join([ss[1], ss[2], "...", ss[end-1], ss[end]], " + ")
    else
        join(ss, " + ")
    end
end

function Base.show(io::IO, ::MIME"text/plain", r::C.Reaction)
    if r.upper_bound > 0.0 && r.lower_bound < 0.0
        arrow = " ↔  "
    elseif r.upper_bound <= 0.0 && r.lower_bound < 0.0
        arrow = " ←  "
    elseif r.upper_bound > 0.0 && r.lower_bound >= 0.0
        arrow = " →  "
    else
        arrow = " →|←  " # blocked reaction
    end
    substrates =
        ["$(-v) $k" for (k, v) in Iterators.filter(((_, v)::Pair -> v < 0), r.stoichiometry)]
    products =
        ["$v $k" for (k, v) in Iterators.filter(((_, v)::Pair -> v >= 0), r.stoichiometry)]

    for fname in fieldnames(C.Reaction)
        if fname == :stoichiometry
            IP._pretty_print_keyvals(
                io,
                "Reaction.$(string(fname)): ",
                _pretty_substances(substrates) * arrow * _pretty_substances(products),
            )
        elseif fname == :gene_association_dnf
            c = isnothing(getfield(r, fname)) ? nothing : join("(".*[join(x, " and ") for x in getfield(r, fname)].*")", " or ")
            IP._pretty_print_keyvals(
                io,
                "Reaction.$(string(fname)): ",
                c,
            )
        elseif fname in (:lb, :ub, :objective_coefficient)
            IP._pretty_print_keyvals(
                io,
                "Reaction.$(string(fname)): ",
                string(getfield(r, fname)),
            )
        else
            IP._pretty_print_keyvals(io, "Reaction.$(string(fname)): ", getfield(r, fname))
        end
    end
end
