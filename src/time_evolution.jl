# solvers for time evolution of quantum systems

using DifferentialEquations

export schrodinger

"""
    schrodinger(cqs::CompositeQSystem, ts::Float64)

Compute the unitary evolution of a CompositeQSystem evaluted at ts.
"""
function schrodinger(cqs::CompositeQSystem, ts::Vector)
    # dU/dt = -iHU
    ham = 2pi * hamiltonian(cqs)
    function schrodinger_eqn(du, u, ham, t)
        du[:] = vec(-1im * ham * u)
    end
    # initial condition of identity
    u0 = eye(Complex128, dim(cqs))
    prob = ODEProblem(schrodinger_eqn, u0, (0, ts[end]), ham)
    sol = solve(prob; saveat=ts)
    return sol
end


"""
    schrodinger(cqs::CompositeQSystem, ts::Float64, ρ_init::Matrix{Complex128})

Compute the unitary evolution of a CompositeQSystem from initial state ρ_init evaluted at ts.
"""
function schrodinger(cqs::CompositeQSystem, ts::Vector, ρ_init::Matrix{Complex128})
    # dρ/dt = -i[H, ρ]
    ham = 2pi * hamiltonian(cqs)
    function schrodinger_eqn(dρ, ρ, ham, t)
        dρ[:] = vec(-1im * (ham*ρ - ρ*ham))
    end
    prob = ODEProblem(schrodinger_eqn, ρ_init, (0, ts[end]), ham)
    sol = solve(prob; saveat=ts)
    return sol
end
