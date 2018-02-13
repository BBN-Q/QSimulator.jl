# ramsey evolution of a single transmon should produce sinusoidal oscillations at the qubit frequency
qubit_freq = 5.0
q0 = FixedDuffingTransmon("q0", qubit_freq, -0.2, 3)
cqs = CompositeQSystem([q0]);
add_hamiltonian!(cqs, hamiltonian(q0), "q0")
times = collect(linspace(0,1,201))
ψ_init = (1/sqrt(2)) * Complex128[1; 1; 0]
ψs = unitary_state(cqs, times, ψ_init);
signal = Float64[real(ψ_init'*ψ) for ψ in ψs]
expected = 0.5 + 0.5*cos.(2π*qubit_freq * (linspace(0,1,201)))
@test isapprox(signal, expected; rtol=1e-4, atol=1e-4)


# rabi flops should produce sinusoidal oscillations between ground and excited state
function qubit_drive(q::QSimulator.QSystem, drive::Function)
    function add_drive_ham!(ham, idxs, t)
        pulse = 2π * drive(t)
        drive_ham = real(pulse) * X(q) + imag(pulse) * Y(q)
        QSimulator.expand_add!(ham, drive_ham, idxs)
    end
end

qubit_freq = 5.0
nutation_freq = 0.02
q0 = FixedDuffingTransmon("q0", qubit_freq, -0.2, 3)
cqs = CompositeQSystem([q0]);
add_hamiltonian!(cqs, hamiltonian(q0), "q0")
add_hamiltonian!(cqs, qubit_drive(q0, t -> nutation_freq*cos(2π*qubit_freq * t)), q0);
ψ_init = Complex128[1; 0; 0]
times = collect(linspace(0,100,101))
ψs = unitary_state(cqs, times, ψ_init);

g_sim = [abs2(s[1]) for s in ψs]
e_sim = [abs2(s[2]) for s in ψs]

# we would normally expect a factor of 0.5 from lab -> rotating frame but by defining the qubit
# drive in terms of X rather than 0.5 X undoes that
g_expected = 0.5 + 0.5*cos.(2π*nutation_freq * times)
@test isapprox(g_sim, g_expected; rtol=1e-2, atol=1e-2)
@test isapprox(e_sim, 1-g_expected; rtol=1e-2, atol=1e-2)
