
"""
$(TYPEDSIGNATURES)

Check if the file at `path` has the expected hash.
"""
function check_data_file_hash(path, expected_checksum)
    actual_checksum = bytes2hex(sha256(open(path)))
    if actual_checksum != expected_checksum
        @error "The downloaded data file `$path' seems to be different from the expected one. Tests will likely fail." actual_checksum expected_checksum
    end
end

"""
$(TYPEDSIGNATURES)

Download the file at `url` and save it at `path`, also check if this file is the
expected file by calling [`check_data_file_hash`](@ref). If the file has already
been downloaded, and stored at `path`, then it is not downloaded again.
"""
function download_data_file(url, path, hash)
    if isfile(path)
        check_data_file_hash(path, hash)
        @info "using cached `$path'"
        return path
    end

    Downloads.download(url, path)
    check_data_file_hash(path, hash)
    return path
end

"""
$(TYPEDSIGNATURES)

Parse a formula in format `C2H6O` into a [`MetaboliteFormula`](@ref), which is
basically a dictionary of atom counts in the molecule.
"""
function parse_formula(f::String)::MetaboliteFormula
    res = Dict{String,Int}()
    pattern = @r_str "([A-Z][a-z]*)([1-9][0-9]*)?"

    for m in eachmatch(pattern, f)
        res[m.captures[1]] = isnothing(m.captures[2]) ? 1 : parse(Int, m.captures[2])
    end

    return res
end

"""
$(TYPEDSIGNATURES)

Apply a function to `x` only if it is not `nothing`.
"""
function maybemap(f, x::Maybe)::Maybe
    isnothing(x) ? nothing : f(x)
end


"""
$(TYPEDSIGNATURES)
Parse `SBML.GeneProductAssociation` structure and convert it to a strictly
positive DNF [`GeneAssociationsDNF`](@ref). Negation (`SBML.GPANot`) is not
supported.
"""
function parse_grr(gpa::SBML.GeneProductAssociation)::GeneAssociationsDNF

    function fold_and(dnfs::Vector{Vector{Vector{String}}})::Vector{Vector{String}}
        if isempty(dnfs)
            [String[]]
        else
            unique(unique(String[l; r]) for l in dnfs[1] for r in fold_and(dnfs[2:end]))
        end
    end

    dnf(x::SBML.GPARef) = [[x.gene_product]]
    dnf(x::SBML.GPAOr) = collect(Set(vcat(dnf.(x.terms)...)))
    dnf(x::SBML.GPAAnd) = fold_and(dnf.(x.terms))
    dnf(x) = throw(
        DomainError(
            x,
            "unsupported gene product association contents of type $(typeof(x))",
        ),
    )
    return dnf(gpa)
end

"""
$(TYPEDSIGNATURES)

Evaluate the SBML `GeneProductAssociation` in the same way as
[`Accessors.eval_reaction_gene_association`](@ref).
"""
function eval_grr(
    gpa::SBML.GeneProductAssociation;
    falses::Maybe{AbstractSet{String}} = nothing,
    trues::Maybe{AbstractSet{String}} = nothing,
)::Bool
    @assert !(isnothing(falses) && isnothing(trues)) "at least one of 'trues' and 'falses' must be defined"

    e(x::SBML.GeneProductAssociation) =
        throw(DomainError(x, "evaluating unsupported GeneProductAssociation contents"))
    e(x::SBML.GPARef) =
        !isnothing(falses) ? !(x.gene_product in falses) : (x.gene_product in trues)
    e(x::SBML.GPAAnd) = all(e, x.terms)
    e(x::SBML.GPAOr) = any(e, x.terms)

    e(gpa)
end

"""
$(TYPEDSIGNATURES)

Convert [`GeneAssociationsDNF`](@ref) to the corresponding `SBML.jl` structure.
"""
function unparse_grr(
    ::Type{SBML.GeneProductAssociation},
    x::GeneAssociationsDNF,
)::SBML.GeneProductAssociation
    SBML.GPAOr([SBML.GPAAnd([SBML.GPARef(j) for j in i]) for i in x])
end

"""
$(TYPEDSIGNATURES)

Parse a DNF gene association rule in format `(YIL010W and YLR043C) or (YIL010W
and YGR209C)` to `GeneAssociationsDNF`. Also accepts `OR`, `|`, `||`, `AND`, `&`,
and `&&`.

# Example
```
julia> parse_grr("(YIL010W and YLR043C) or (YIL010W and YGR209C)")
2-element Array{Array{String,1},1}:
 ["YIL010W", "YLR043C"]
 ["YIL010W", "YGR209C"]
```
"""
parse_grr(s::String) = maybemap(parse_grr, parse_grr_to_sbml(s))::Maybe{GeneAssociationsDNF}

"""
PikaParser grammar for stringy GRR expressions.
"""
const grr_grammar = begin
    # characters that typically form the identifiers
    isident(x::Char) =
        isletter(x) ||
        isdigit(x) ||
        x == '_' ||
        x == '-' ||
        x == ':' ||
        x == '.' ||
        x == '\'' ||
        x == '[' ||
        x == ']' ||
        x == '\x03' # a very ugly exception for badly parsed MAT files

    # scanner helpers
    eat(p) = m -> begin
        last = 0
        for i in eachindex(m)
            p(m[i]) || break
            last = i
        end
        last
    end

    # eat one of keywords
    kws(w...) = m -> begin
        last = eat(isident)(m)
        m[begin:last] in w ? last : 0
    end

    PP.make_grammar(
        [:expr],
        PP.flatten(
            Dict(
                :space => PP.first(PP.scan(eat(isspace)), PP.epsilon),
                :id => PP.scan(eat(isident)),
                :orop =>
                    PP.first(PP.tokens("||"), PP.token('|'), PP.scan(kws("OR", "or"))),
                :andop => PP.first(
                    PP.tokens("&&"),
                    PP.token('&'),
                    PP.scan(kws("AND", "and")),
                ),
                :expr => PP.seq(:space, :orexpr, :space, PP.end_of_input),
                :orexpr => PP.first(
                    :or => PP.seq(:andexpr, :space, :orop, :space, :orexpr),
                    :andexpr,
                ),
                :andexpr => PP.first(
                    :and => PP.seq(:baseexpr, :space, :andop, :space, :andexpr),
                    :baseexpr,
                ),
                :baseexpr => PP.first(
                    :id,
                    :parenexpr => PP.seq(
                        PP.token('('),
                        :space,
                        :orexpr,
                        :space,
                        PP.token(')'),
                    ),
                ),
            ),
            Char,
        ),
    )
end

grr_grammar_open(m, _) =
    m.rule == :expr ? Bool[0, 1, 0, 0] :
    m.rule == :parenexpr ? Bool[0, 0, 1, 0, 0] :
    m.rule in [:or, :and] ? Bool[1, 0, 0, 0, 1] :
    m.rule in [:andexpr, :orexpr, :notexpr, :baseexpr] ? Bool[1] :
    (false for _ in m.submatches)

grr_grammar_fold(m, _, subvals) =
    m.rule == :id ? SBML.GPARef(m.view) :
    m.rule == :and ? SBML.GPAAnd([subvals[1], subvals[5]]) :
    m.rule == :or ? SBML.GPAOr([subvals[1], subvals[5]]) :
    m.rule == :parenexpr ? subvals[3] :
    m.rule == :expr ? subvals[2] : isempty(subvals) ? nothing : subvals[1]

"""
$(TYPEDSIGNATURES)

Internal helper for parsing the string GRRs into SBML data structures. More
general than [`parse_grr`](@ref).
"""
function parse_grr_to_sbml(str::String)::Maybe{SBML.GeneProductAssociation}
    all(isspace, str) && return nothing
    tree = PP.parse_lex(grr_grammar, str)
    match = PP.find_match_at!(tree, :expr, 1)
    if match > 0
        return PP.traverse_match(
            tree,
            match,
            open = grr_grammar_open,
            fold = grr_grammar_fold,
        )
    else
        throw(DomainError(str, "cannot parse GRR"))
    end
end

"""
$(TYPEDSIGNATURES)

Converts a nested string gene reaction array back into a gene reaction rule
string.

# Example
```
julia> unparse_grr(String, [["YIL010W", "YLR043C"], ["YIL010W", "YGR209C"]])
"(YIL010W && YLR043C) || (YIL010W && YGR209C)"
```
"""
function unparse_grr(
    ::Type{String},
    grr::GeneAssociationsDNF;
    and = " && ",
    or = " || ",
)::String
    return join(("(" * join(gr, and) * ")" for gr in grr), or)
end

"""
$(TYPEDSIGNATURES)

Unfortunately, many model types that contain dictionares do not have
standardized field names, so we need to try a few possibilities and guess the
best one.
"""
function guesskey(avail, possibilities)
    x = intersect(possibilities, avail)

    if isempty(x)
        @debug "could not find any of keys: $possibilities"
        return nothing
    end

    if length(x) > 1
        @debug "Possible ambiguity between keys: $x"
    end
    return x[1]
end

"""
$(TYPEDSIGNATURES)

Return `fail` if key in `keys` is not in `collection`, otherwise
return `collection[key]`. Useful if may different keys need to be
tried due to non-standardized model formats.
"""
function gets(collection, fail, keys)
    for key in keys
        haskey(collection, key) && return collection[key]
    end
    return fail
end

