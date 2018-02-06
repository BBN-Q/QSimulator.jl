export CompositeQSystem,
       hamiltonian

# tensor products of quantum systems
mutable struct CompositeQSystem
    subsystems::Vector{Any}
    fixed_Hs::Vector{Any} # tuple of Matrix and exansion indices
    parametric_Hs::Vector{Any} # tuple of Functions and expansion indices
    dim::Int
end

CompositeQSystem(qs) = CompositeQSystem(qs, [], [], prod(dim(q) for q in qs))

""" Calculate the drift or natural Hamiltonian of a CompositeQSystem """
function hamiltonian(cqs::CompositeQSystem)
    ham = zeros(Complex128, (dim(cqs), dim(cqs)))
    for h = cqs.fixed_Hs
        expand_add!(ham, new_ham, expand_idxs)
    end
end


""" In place addition of an operator embedded into a larger Hilbert space given a set of expansion indices"""
function expand_add!(op, new_op, expand_idxs)
    for ct = 1:length(expand_idxs)
        for idx = expand_idxs[ct]
            op[idx] += new_op[ct]
        end
    end
end
