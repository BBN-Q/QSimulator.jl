#=
Helper function for creating and manipulating qubit paulis.

Some code from https://github.com/blakejohnson/Qlab.jl 

=#

export ## Methods
		pauli_mats,
		pauli_strs,
		pauli_decomp

const pX = complex128([[0.0 1.0;1.0 0.0]]) 
const pY = complex128([[0.0 -1.0im;1.0im 0.0]]) 
const pZ = complex128([[1.0 0.0;0.0 -1.0]]) 
const pI = eye(Complex128, 2) 

function base_paulis(dim::Int)
	#Pauli projectors onto the bottom 2 levels of a dim-level system
	map({pI, pX, pY, pZ}) do p
		fullPauli = zeros(Complex128, dim, dim)
		fullPauli[1:2, 1:2] = p
		fullPauli
	end
end

function pauli_mats(n::Int, dim=2)
	#Multi-qubit Pauli operators for n qubits
	@assert n > 0 "You need at least 1 qubit!"
	n == 1 && return base_paulis(dim)
	map(x->kron(x...), product(base_paulis(dim), pauli_mats(n-1, dim)))
end

function pauli_mats(dims::Vector{Int})
	#Multi-qubit pauli operators for n multi-level systems
	length(dims) == 1 && return base_paulis(dims[1])
	map(x->kron(x...), product(base_paulis(dims[1]), pauli_mats(dims[2:end])))
end

function pauli_strs(n::Int; sortByWeight=false)
	#Return the strings of pauli operators, e.g.: II, IX, IY, IZ, ....
	if n == 0
		return [""]
	end
	ops = ["I" "X" "Y" "Z"]
	paulis = vec(ops .* pauli_strs(n-1))

	if (sortByWeight)
		sort!(paulis, by=hamming)
	end

	paulis
end

function hamming(pauliStr)
	# compute the Hamming weight of a Pauli operator (number of non-identity components)
	weight = 0
	for pauli in pauliStr
		weight += (pauli == 'I') ? 0 : 1
	end
	weight
end

function pauli_decomp(rho::Array{Complex128,2}, dims::Array{Int, 1}, cutoff=1e-3)
	#Provide a decomposition of a density matrix into qubit pauli operators
	dimScale = 2^length(dims)
	for (str, mat) in zip(pauli_strs(length(dims)), pauli_mats(dims))
		val = real(trace(rho*mat)) / dimScale
		if val > cutoff
			println(str, ": ", val)
		end
	end

end