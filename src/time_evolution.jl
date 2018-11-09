using LinearAlgebra: I, rmul!, mul!
import Base.Iterators

# We the more specific OrdinaryDiffEq for a significantly lighter dependency tree and build time.
# However, this means we loose the automated solver picking. For now we use a recommended Tsit5
# solver but it is not clear whether it is the right choice for the stiffness of the
# Schrodinger/Lindblad equations.
# * `Rodas4/5` does not appear to support complex equations.
# * `DP5` used by QuantumOptics.jl also gave strange errors
# * the `reltol` and `abstol` were chosen somewhat arbirtarily to be "good enough"
using OrdinaryDiffEq: ODEProblem, solve, Tsit5

export unitary_propagator, unitary_state, me_propagator, me_state

"""
    unitary_propagator(cqs::CompositeQSystem, ts::AbstractVector{<:Real})

Compute the unitary propagator of a CompositeQSystem at given times by solving the differential
equation `dU/dt = -iHU`.

## args
* `cqs`: a CompositeQSystem.
* `ts`: a sorted array of times at which to compute the unitary propagator.

## returns
An array of unitary propagators at the specified times.
"""
function unitary_propagator(cqs::CompositeQSystem, ts::AbstractVector{<:Real})
    @assert issorted(ts)
    function ode(du, u, p, t)
        ham = p[3] # preallocated workspace array
        ham .= p[2] # start from fixed_ham
        add_parametric_hamiltonians!(ham, p[1], t)
        rmul!(ham, -2π * 1im)
        mul!(du, ham, u)
    end
    fixed_ham = hamiltonian(cqs)
    d = dim(cqs)
    u0 = Matrix{ComplexF64}(I, d, d) # start with identity
    work_ham = similar(fixed_ham) # scratch space
    prob = ODEProblem(ode, u0, (float(ts[1]), float(ts[end])), (cqs, fixed_ham, work_ham)spe
    sol = solve(prob, Tsit5(); saveat=ts, save_start=true, reltol=1e-6, abstol=1e-8)
    return sol.u
end

unitary_propagator(cqs::CompositeQSystem, t::Real) = unitary_propagator(cqs, [0.0, t])[end]

"""
    unitary_state(cqs::CompositeQSystem, ts::AbstractVector{<:Real}, ψ0::Vector{<:Number})

Compute the unitary state evolution of a `CompositeQSystem` from a given initial state and at given
times by solving the differential equation `dψ/dt = -iHψ`.

## args
* `cqs`: a CompositeQSystem.
* `ts`: a sorted array of times at which to compute the state.
* `ψ0`: a vector indicating the initial state.

## returns
An array of state vectors for the system at the specified times.
"""
function unitary_state(cqs::CompositeQSystem, ts::AbstractVector{<:Real}, ψ0::Vector{<:Number})
    @assert issorted(ts)
    function ode(dψ, ψ, p, t)
        ham = p[3] # preallocated workspace array
        ham .= p[2] # start from fixed_ham
        add_parametric_hamiltonians!(ham, p[1], t)
        rmul!(ham, -2π * 1im)
        mul!(dψ, ham, ψ)
    end
    fixed_ham = hamiltonian(cqs)
    work_ham = similar(fixed_ham)
    prob = ODEProblem(ode, ψ0, (float(ts[1]), float(ts[end])), (cqs, fixed_ham, work_ham))
    sol = solve(prob, Tsit5(); saveat=ts, save_start=true, reltol=1e-6, abstol=1e-8)
    return sol.u
end

unitary_state(cqs::CompositeQSystem, t::Real, ψ0::Vector{<:Number}) = unitary_state(cqs, [0.0, t], ψ0)[end]

"""
    me_propagator(cqs::CompositeQSystem, ts::AbstractVector{<:Real})

Compute the master equation propagator evolution of a `CompositeQSystem` at given times by solving
the vectorized master equation `du/dt = (-i(I ⊗ H - transpose(H) ⊗ I) + Σ conj(L) ⊗ L - I ⊗ L^†L/2 -
transpose(L^†L/2) ⊗ I)u` which can be derived using the identity `vec(AXB) = (transpose(B) ⊗
A)vec(X)`.

## args
* `cqs`: a CompositeQSystem.
* `ts`: a sorted array of times at which to compute the density matrix.

## returns
An array of propagators for the system at the specified times.
"""
function me_propagator(cqs::CompositeQSystem, ts::AbstractVector{<:Real})
    @assert issorted(ts)
    function ode(du, u, p, t)
        ham = p[3] # preallocated workspace array
        ham .= p[2] # start from fixed_ham
        add_parametric_hamiltonians!(ham, p[1], t)
        d = size(ham, 1)
        I_mat = Matrix{ComplexF64}(I, d, d)
        mul!(du, -2π * 1im * (I_mat ⊗ ham - transpose(ham) ⊗ I_mat), u)
        lind_mat = p[5] # preallocated workspace array
        for (lind_op, idxs) in Iterators.flatten((p[1].fixed_Ls, ((l(t), idxs) for (l, idxs) in p[1].parametric_Ls)))
            lind_mat .= p[4] # start with empty array
            embed_add!(lind_mat, lind_op, idxs)
            du .+= 2π * (conj(lind_mat) ⊗ lind_mat .- 0.5 .* I_mat ⊗ (lind_mat'*lind_mat) .- 0.5 .* transpose(lind_mat'*lind_mat) ⊗ I_mat) * u
        end
    end
    fixed_ham = hamiltonian(cqs)
    work_ham = similar(fixed_ham)
    bare_lind = zeros(ComplexF64, size(fixed_ham))
    work_lind = similar(fixed_ham)
    d = dim(cqs)^2
    u0 = Matrix{ComplexF64}(I, d, d) # start with identity
    prob = ODEProblem(ode, u0, (float(ts[1]), float(ts[end])), (cqs, fixed_ham, work_ham, bare_lind, work_lind))
    sol = solve(prob, Tsit5(); saveat=ts, save_start=true, reltol=1e-6, abstol=1e-8)
    return sol.u
end

me_propagator(cqs::CompositeQSystem, t::Real) = me_propagator(cqs, [0.0, t])[end]

"""
    me_state(cqs::CompositeQSystem, ts::AbstractVector{<:Real}, ρ0::Matrix{<:Number})

Compute the master equation state evolution of a `CompositeQSystem` from a given initial density
matrix and at given times by solving the differential equation `dρ/dt = -i[H, ρ] + Σ LρL^† - L^†Lρ/2
- ρL^†L/2`.

## args
* `cqs`: a CompositeQSystem.
* `ts`: a sorted array of times at which to compute the density matrix.
* `ρ0`: a matrix indicating the initial density matrix.

## returns
An array of density matrices for the system at the specified times.
"""
function me_state(cqs::CompositeQSystem, ts::AbstractVector{<:Real}, ρ0::Matrix{<:Number})
    @assert issorted(ts)
    function ode(dρ, ρ, p, t)
        ham = p[3] # preallocated workspace array
        ham .= p[2] # start from fixed_ham
        add_parametric_hamiltonians!(ham, p[1], t)
        dρ .= -2π * 1im * (ham*ρ - ρ*ham)
        lind_mat = p[5] # preallocated workspace array
        for (index, (lind_op, idxs)) in enumerate([p[1].fixed_Ls; p[1].parametric_Ls])
            lind_mat .= p[4] # start with empty array
            l = index <= length(p[1].fixed_Ls) ? lind_op : lind_op(t)
            embed_add!(lind_mat, l, idxs)
            dρ .+= 2π * (lind_mat*ρ*lind_mat' .- .5 .* lind_mat'*lind_mat*ρ .- .5 .* ρ*lind_mat'*lind_mat)
        end
    end
    fixed_ham = hamiltonian(cqs)
    work_ham = similar(fixed_ham)
    bare_lind = zeros(ComplexF64, size(fixed_ham))
    work_lind = similar(fixed_ham)
    prob = ODEProblem(ode, ρ0, (float(ts[1]), float(ts[end])), (cqs, fixed_ham, work_ham, bare_lind, work_lind))
    sol = solve(prob, Tsit5(); saveat=ts, save_start=true, reltol=1e-6, abstol=1e-8)
    return sol.u
end

me_state(cqs::CompositeQSystem, t::Real, ρ0::Matrix{<:Number}) = me_state(cqs, [0.0, t], ρ0)[end]
