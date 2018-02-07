# solvers for time evolution of quantum systems

using DifferentialEquations

export schrodinger

"""
    schrodinger(cqs::CompositeQSystem, ts::Float64)

Compute the unitary evolution of a CompositeQSystem evaluted at ts.
"""
function schrodinger(cqs::CompositeQSystem, ts::Vector)
    # dU = -iHU
    ham = hamiltonian(cqs)
    function schrodinger_eqn(du, u, ham, t)
        du[:] = vec(-1im * ham * u)
    end
    # initial condition of identity
    u0 = eye(Complex128, dim(cqs))
    prob = ODEProblem(schrodinger_eqn, u0, (0, ts[end]), ham)
    sol = solve(prob)
    return sol
end
