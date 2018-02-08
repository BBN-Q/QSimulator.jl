# solvers for time evolution of quantum systems

using DifferentialEquations

export unitary_propagator,
       unitary_state

"""
    schrodinger(cqs::CompositeQSystem, ts::Float64; u0::Matrix=Matrix{Complex128}(0,0))

Compute the unitary propagator evolution of a CompositeQSystem evaluted at ts.
"""
function unitary_propagator(cqs::CompositeQSystem, ts::Vector; u0=Matrix{Complex128}(0,0))
    # schrodinger differential equation for unitary with in place update
    # dU/dt = -iHU
    function ode(du, u, ham, t)
        du[:] = vec(-1im * ham * u)
    end
    # scale Hamiltonian from Hz to rad.
    ham = 2pi * hamiltonian(cqs)
    # if initial condition not passed start with identity
    if isempty(u0)
        u0 = eye(Complex128, dim(cqs))
    end
    prob = ODEProblem(ode, u0, (0, ts[end]), ham)
    sol = solve(prob; saveat=ts)
    sol.u
end


"""
    unitary_state(cqs::CompositeQSystem, ts::Float64, ρ0::Matrix)

Compute the unitary state evolution of a CompositeQSystem from initial state ρ_init evaluted at ts.
"""
function unitary_state(cqs::CompositeQSystem, ts::Vector, ρ0::Matrix)
    # schrodinger differential equation for density matrix with in place update
    # dρ/dt = -i[H, ρ]
    function ode(dρ, ρ, ham, t)
        dρ[:] = vec(-1im * (ham*ρ - ρ*ham))
    end
    # scale Hamiltonian from Hz to rad.
    ham = 2pi * hamiltonian(cqs)
    prob = ODEProblem(ode, ρ0, (0, ts[end]), ham)
    sol = solve(prob; saveat=ts)
    sol.u
end
