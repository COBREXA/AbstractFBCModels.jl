
"""
$(TYPEDSIGNATURES)

Rebuild the given [`Model`](@ref) with all identifiers changed accordingly to
the given individual identifier-mapping functions in arguments. These take a
single `String` identifier as an argument, and return a new `String`
identifier. By default, the identifier maps are `identity`, i.e., no
identifiers are changed.

The identifier maps *must be injective and pure within `identifier_map`*, i.e.,
they must not create identifier name conflicts, and the results returned for a
given argument must be the same during the whole call of `identifier_map`.
Errors stemming from use of non-injective and impure identifier maps are not
handled.

Internal data structures are copied only by reference wherever possible.
"""
function identifier_map(
    m::Model;
    reaction_map = identity,
    metabolite_map = identity,
    gene_map = identity,
    coupling_map = identity,
    compartment_map = identity,
)
    return Model(;
        reactions = Dict(
            reaction_map(k) => Reaction(;
                name = v.name,
                lower_bound = v.lower_bound,
                upper_bound = v.upper_bound,
                stoichiometry = Dict(
                    metabolite_map(sk) => sv for (sk, sv) in v.stoichiometry
                ),
                objective_coefficient = v.objective_coefficient,
                gene_association_dnf = isnothing(v.gene_association_dnf) ? nothing :
                                       [
                    String[gene_map(g) for g in gs] for gs in v.gene_association_dnf
                ],
                annotations = v.annotations,
                notes = v.notes,
            ) for (k, v) in m.reactions
        ),
        metabolites = Dict(
            metabolite_map(k) => Metabolite(;
                name = v.name,
                compartment = isnothing(v.compartment) ? nothing :
                              compartment_map(v.compartment),
                formula = v.formula,
                charge = v.charge,
                balance = v.balance,
                annotations = v.annotations,
                notes = v.notes,
            ) for (k, v) in m.metabolites
        ),
        genes = Dict(gene_map(k) => v for (k, v) in m.genes),
        couplings = Dict(
            coupling_map(k) => Coupling(;
                name = v.name,
                reaction_weights = Dict(
                    reaction_map(rk) => rv for (rk, rv) in v.reaction_weights
                ),
                lower_bound = v.lower_bound,
                upper_bound = v.upper_bound,
                annotations = v.annotations,
                notes = v.notes,
            ) for (k, v) in m.couplings
        ),
    )
end
