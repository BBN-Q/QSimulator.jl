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
export floquet_propagator, floquet_rise_fall_propagator, choose_times_floquet, decompose_times_floquet

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
    prob = ODEProblem(ode, u0, (float(ts[1]), float(ts[end])), (cqs, fixed_ham, work_ham))
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

######################################################
# Define a `propagator_function` to be a function that takes a CompositeQSystem
# and an array of times and returns an array of propagators at those times
# starting at the first given time.
######################################################

# This value is used in periodic problems to produce a small time
# `dt = TIME_TOL_FRACTION * t_period` where t_period is the periodicity.
# This value is used as a small time-scale on which errors in the time
# do not matter. The default value here says we do not need precision on
# the time better then the periodicity over ten billion.
TIME_TOL_FRACTION = 1e-10

"""
    floquet_propagator(propagator_func::Function, t_period::Real; time_tol_fraction::Real=TIME_TOL_FRACTION)

Given a propagator function and a time period, create a new propagator
function that applies correctly to CompositeQSystems that are periodic
with the given time period. It makes use of the identity
`U(nτ + dt) = U(dt)U(τ)^n` where `U` is the propagator for a system with period `τ`.

## args
* `propagator_func`: a propagator function, e.g. `unitary_propagator`.
* `t_period`: the periodicity of the system.
* `time_tol_fraction`: a tolerance to use in unique_tol for the times mod the period
    expressed as a fraction of `t_period`.

## returns
A new propagator function.
"""
function floquet_propagator(propagator_func::Function, t_period::Real; time_tol_fraction::Real=TIME_TOL_FRACTION)
    function p(cqs::CompositeQSystem, ts::Vector{<:Real})
        @assert issorted(ts)
        quotients, unique_remainders, unique_inds = decompose_times_floquet(ts, t_period, time_tol_fraction=time_tol_fraction)
        # perform time evolution for unique remainders as well as one full period
        us_remainders = propagator_func(cqs, [unique_remainders; t_period + ts[1]])
        u_period = us_remainders[end]
        # create desired unitaries
        us = Matrix{ComplexF64}[]
        quotient = quotients[1]
        u_floquet = u_period^quotient
        for i in 1:length(ts)
            if quotients[i] > quotient
                u_floquet *= u_period^(quotients[i] - quotient)
                quotient = quotients[i]
            end
            push!(us, us_remainders[unique_inds[i]] * u_floquet)
        end
        return us
    end
    return p
end

"""
    choose_times_floquet(center::Real, width::Real, t_period::Real, dt::Real)

Choose times for integration that maximize remainder overlap for Floquet integration.
`center` is guaranteed to be one of the times.

## args
* `center`: the center of the array of times.
* `width`: the desired approximate width of the array of times.
* `t_period`: the periodicity of the system.
* `dt`: the desired approximate time increment.

## returns
An array of times.
"""
function choose_times_floquet(center::Real, width::Real, t_period::Real, dt::Real)
    dt = t_period / ceil(Int, t_period / dt) # reset dt so it divides time_period, making it smaller
    num_times = floor(Int, width / dt)
    num_times += iseven(num_times)
    width = dt * num_times
    times = collect(range(center - width/2, stop=center + width/2, length=num_times))
    times[ceil(Int, num_times/2)] = center
    return times
end

"""
    decompose_times_floquet(ts::Vector{<:Real}, t_period::Real; time_tol_fraction::Real=TIME_TOL_FRACTION)

Decompose the given times as `ts = quotients * t_period + unique_remainders[unique_inds]`
where `unique_remainders` is as small of an array as possible and has `t[1]` as its first element.

## args
* `ts`: an array of times.
* `t_period`: the periodicity of the system.
* `time_tol_fraction`: a tolerance to use in unique_tol for the times mod the period
    expressed as a fraction of `t_period`.

## returns
`quotients`, `unique_remainders`, and `unique_inds`.
"""
function decompose_times_floquet(ts::Vector{<:Real}, t_period::Real; time_tol_fraction::Real=TIME_TOL_FRACTION)
    # find unique times mod a period
    quotients = fld.(ts .- ts[1], t_period)
    remainders = mod.(ts .- ts[1], t_period) .+ ts[1]
    # remove collisions in the times mod a period up to time_tol
    unique_remainders, unique_inds = unique_tol(remainders, time_tol_fraction * t_period)
    # sort the unique remainders
    sort_inds = sortperm(unique_remainders)
    unique_remainders = unique_remainders[sort_inds]
    unique_inds = sortperm(sort_inds)[unique_inds]
    return quotients, unique_remainders, unique_inds
end

"""
    floquet_rise_fall_propagator(propagator_func::Function, t_period::Real, rise_time::Real, fall_time::Real; time_tol_fraction::Real=TIME_TOL_FRACTION)

Create a propagator function that uses the floquet evolution but also
allows a rise time and a fall time during which the pulse is not periodic.
The rise time is assumed to be at the beginning of the vector `ts` and the
fall time is assumed to be at the end.

## args
* `propagator_func`: a propagator function, e.g. `unitary_propagator`.
* `t_period`: the periodicity of the system.
* `rise_time`: the rise time of the pulse (the time at the beginning that is non-periodic).
* `fall_time`: the fall time of the pulse (the time at the end that is non-periodic).
* `time_tol_fraction`: a tolerance to use in unique_tol for the times mod the period
    expressed as a fraction of `t_period`.

## returns
A propagator function. The times passed into this
"""
function floquet_rise_fall_propagator(propagator_func::Function, t_period::Real, rise_time::Real, fall_time::Real; time_tol_fraction::Real=TIME_TOL_FRACTION)
    @assert all([t_period, rise_time, fall_time] .>= 0.0)
    floquet_prop = floquet_propagator(propagator_func, t_period, time_tol_fraction=time_tol_fraction)
    function p(cqs::CompositeQSystem, ts::Vector{<:Real})
        t0, t1 = ts[1] + rise_time, ts[end] - fall_time
        @assert t0 <= t1
        us_risetime = propagator_func(cqs, [ts[ts .< t0]; t0])
        u0 = us_risetime[end]
        us = us_risetime[1:end-1]
        us_floquet = floquet_prop(cqs, [t0; ts[t0 .<= ts .<= t1]; t1])
        u1 = us_floquet[end] * u0
        append!(us, [u * u0 for u in us_floquet[2:end-1]])
        us_falltime = propagator_func(cqs, [t1; ts[t1 .< ts]])
        append!(us, [u * u1 for u in us_falltime[2:end]])
        return us
    end
    return p
end

"""
    unique_tol(ts::Vector{<:Real}, dt::Real)

Compute an array `a` of values where each is within `dt/2` of some value in `ts`
and no two differ by less than `dt`. Further return a vector of indices `inds`
such that `a[inds]` is an approximation of `ts` up to `dt`.

## args
* `ts`: an array of numbers.
* `dt`: a small number within which errors are unimportant.

## returns
An array of the unique values and an array of indices.
"""
function unique_tol(ts::Vector{<:Real}, dt::Real)
    vals = round.(ts ./ dt) .* dt
    unique_vals = unique(vals)
    unique_inds = [findfirst(unique_vals .== v) for v in vals]
    return unique_vals, unique_inds
end
