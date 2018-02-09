# solvers for time evolution of quantum systems

using DifferentialEquations

import QSimulator.add_parametric_hamiltonians!

export unitary_propagator,
       unitary_state

"""
    schrodinger(cqs::CompositeQSystem, ts::Float64; u0::Matrix=Matrix{Complex128}(0,0))

Compute the unitary propagator evolution of a CompositeQSystem evaluted at ts.
"""
function unitary_propagator(cqs::CompositeQSystem, ts::Vector; u0=Matrix{Complex128}(0,0))
    # schrodinger differential equation for unitary with in place update
    # dU/dt = -iHU
    function ode(du, u, p, t)
        ham = p[1] # fixed_ham
        add_parametric_hamiltonians!(ham, p[2], t)
        du[:] = vec(-1im * ham * u)
    end
    # scale Hamiltonian from Hz to rad.
    fixed_ham = 2pi * hamiltonian(cqs)
    # if initial condition not passed start with identity
    if isempty(u0)
        u0 = eye(Complex128, dim(cqs))
    end
    prob = ODEProblem(ode, u0, (0, ts[end]), (fixed_ham, cqs))
    sol = solve(prob; saveat=ts)
    sol.u
end


"""
    unitary_state(cqs::CompositeQSystem, ts::Float64, ψ0::Vector)

Compute the unitary state evolution of a CompositeQSystem from initial state ψ0 evaluted at ts.
"""
function unitary_state(cqs::CompositeQSystem, ts::Vector, ψ0::Vector)
    # schrodinger differential equation for state vector with in place update
    # dψ/dt = -iHψ
    function ode(dψ, ψ, p, t)
        ham = copy(p[1])
        add_parametric_hamiltonians!(ham, p[2], t)
        dψ[:] = vec(-1im * ham * ψ)
    end
    # scale Hamiltonian from Hz to rad.
    fixed_ham = 2pi * hamiltonian(cqs)
    prob = ODEProblem(ode, ψ0, (0, ts[end]), (fixed_ham, cqs))
    sol = solve(prob; saveat=ts, reltol=1e-6)
    sol.u
end


"""
    unitary_state(cqs::CompositeQSystem, ts::Float64, ρ0::Matrix)

Compute the unitary state evolution of a CompositeQSystem from initial density matrix ρ0 evaluted at ts.
"""
function unitary_state(cqs::CompositeQSystem, ts::Vector, ρ0::Matrix)
    # schrodinger differential equation for density matrix with in place update
    # dρ/dt = -i[H, ρ]
    function ode(dρ, ρ, p, t)
        ham = p[1]
        add_parametric_hamiltonians!(ham, p[2], t)
        dρ[:] = vec(-1im * (ham*ρ - ρ*ham))
    end
    # scale Hamiltonian from Hz to rad.
    fixed_ham = 2pi * hamiltonian(cqs)
    prob = ODEProblem(ode, ρ0, (0, ts[end]), (fixed_ham, cqs))
    sol = solve(prob; saveat=ts)
    sol.u
end
