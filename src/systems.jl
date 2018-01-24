export Resonator,
       FixedTransmon

# atomic  quantum systems such as resonators or transmons
abstract type QSystem end

label(q::QSystem) = q.label
dim(q::QSystem) = q.dim

# some basic operators
raising(q::QSystem) = diagm(sqrt.(1:(dim(q)-1)), -1)
const create = raising
lowering(q::QSystem) = diagm(sqrt.(1:(dim(q)-1)), 1)
const destroy = lowering
number(q::QSystem) = raising(q) * lowering(q)
X(q::QSystem) = raising(q) + lowering(q)
Y(q::QSystem) = 1im*(raising(q) - lowering(q))

# linear resonator
mutable struct Resonator <: QSystem
    label::AbstractString
    frequency::Float64
    dim::Int
end
hamiltonian(r::Resonator) = r.frequency * number(r)

mutable struct FixedTransmon <: QSystem
    label::AbstractString
    frequency::Float64
    anharmonicity::Float64
    dim::Int
end
hamiltonian(t::FixedTransmon) = (t.frequency - 0.5*t.anharmonicity)*number(t) + 0.5*t.anharmonicity * number(t)^2


# tensor products of quantum systems
mutable struct CompositeQSystem
    subsystems::Vector{QSystem}
    hamiltonians::Vector{Any}
end
