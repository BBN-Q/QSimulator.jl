using Test, QSimulator
using QSimulator: index
using LsqFit: curve_fit

############################ Ramsey ##########################

# ramsey evolution of a single transmon should produce sinusoidal oscillations at the qubit frequency
@testset "Ramsey" begin
    dims = 3
    qubit_freq = 5.0
    anharm = -0.2
    q0 = DuffingTransmon("q0", dims, DuffingSpec(qubit_freq, anharm))
    cqs = CompositeQSystem([q0])
    add_hamiltonian!(cqs, q0)
    times = collect(range(0,stop=1,length=201))
    ψ_init = (1/sqrt(2)) * ComplexF64[1; 1; 0]
    ψs = unitary_state(cqs, times, ψ_init)
    signal = Float64[real(ψ_init'*ψ) for ψ in ψs]
    expected = 0.5 .+ 0.5*cos.(2π*qubit_freq * times)
    @test isapprox(signal, expected, rtol=1e-3)

    # check that the propagator gives the same result
    us = unitary_propagator(cqs, times)
    ψs = [u*ψ_init for u in us]
    signal = Float64[real(ψ_init'*ψ) for ψ in ψs]
    @test isapprox(signal, expected, rtol=1e-4)
end

########################### Rabi ##########################

# rabi flops should produce sinusoidal oscillations between ground and excited state
@testset "Rabi" begin
    dims = 3
    qubit_freq = 5.0
    anharm = -0.2
    nutation_freq = 0.02
    q0 = DuffingTransmon("q0", dims, DuffingSpec(qubit_freq, anharm))
    cqs = CompositeQSystem([q0])
    add_hamiltonian!(cqs, q0)
    add_hamiltonian!(cqs, dipole_drive(q0, t -> nutation_freq*cos(2π*qubit_freq * t)), q0)
    ψ_init = ComplexF64[1; 0; 0]
    times = collect(range(0,stop=100,length=101))
    ψs = unitary_state(cqs, times, ψ_init)

    g_sim = [abs2(s[1]) for s in ψs]
    e_sim = [abs2(s[2]) for s in ψs]

    # we would normally expect a factor of 0.5 from lab -> rotating frame but by defining the qubit
    # drive in terms of X rather than 0.5 X undoes that
    g_expected = 0.5 .+ 0.5*cos.(2π*nutation_freq * times)
    @test isapprox(g_sim, g_expected, rtol=1e-2)
    @test isapprox(e_sim, 1 .- g_expected, rtol=1e-2)

    # check that the propagator gives the same result
    # use floquet_propagator to speed it up since the Hamiltonian is periodic
    floquet_prop = floquet_propagator(unitary_propagator, 1/qubit_freq)
    us = floquet_prop(cqs, times)
    ψs = [u*ψ_init for u in us]
    g_sim = [abs2(s[1]) for s in ψs]
    e_sim = [abs2(s[2]) for s in ψs]
    @test isapprox(g_sim, g_expected; rtol=1e-2)
    @test isapprox(e_sim, 1 .- g_expected; rtol=1e-2)
end

########################### Parmetric Flux Drive ##########################

# parametric flux drive should produce flops between resonant states
# parameters from Blue Launch paper
@testset "parametric flux drive" begin
    # create Hamiltonian
    dims = [3,3]
    duffing = DuffingSpec(3.94015, -0.1807)
    transmon = TransmonSpec(0.172, 12.71, 3.69)
    q0 = DuffingTransmon("q0", dims[1], duffing)
    q1 = PerturbativeTransmon("q1", dims[2], transmon)
    # go to rotating frame at duffing.frequency for speed
    # a doubly rotating frame is also possible, but then you can't use floquet
    q0 = RotatingFrameSystem(q0, duffing.frequency)
    q1 = RotatingFrameSystem(q1, duffing.frequency)
    cqs = CompositeQSystem([q0, q1])
    add_hamiltonian!(cqs, q0)
    g = 0.006
    add_hamiltonian!(cqs, t -> g/2 * X_Y([q0, q1], [duffing.frequency, duffing.frequency] * t), [q0,q1])

    # compute iSWAP parameters
    mod_park = 0.0
    mod_amp = 0.323
    ϕ(t, freq) = mod_park + mod_amp * sin(2π * freq * t)
    f(t, freq) = perturbative_transmon_freq(transmon.EC, transmon.EJ1, transmon.EJ2, ϕ(t, freq))
    dummy_freq = 1.0 # use this for computing average frequency
    avg_freq = real(fourier_coefficient(t -> f(t, dummy_freq), dummy_freq, 0))
    detuning = avg_freq - duffing.frequency
    harmonic = 2
    freq = -detuning/harmonic
    fs = FourierSeries(t -> 2π * f(t, freq), freq, collect(-50:50))
    g_eff = abs(g * rotating_frame_series(fs, [harmonic]).terms[harmonic])
    t_gate = 1/(4 * g_eff)

    # create basis states
    basis = TensorProductBasis(dims)
    ψ₀₁ = TensorProductBasisState(basis, (0,1))
    ψ₁₀ = TensorProductBasisState(basis, (1,0))

    # perform time evolution
    add_hamiltonian!(cqs, parametric_drive(q1, t -> ϕ(t, freq)), q1)
    times = collect(range(0, stop=2 * t_gate, length=1000))
    prop = floquet_propagator(unitary_propagator, 1/abs(freq))
    ψs = [u * ψ0 for u in prop(cqs, times)]
    pop_10 = [abs2(ψ[index(ψ₁₀)]) for ψ in ψs]
    pop_01 = [abs2(ψ[index(ψ₀₁)]) for ψ in ψs]
    # population should oscillate between 01 and 10
    expected_10 = cos.(2π * g_eff * times).^ 2
    # there is additional lab frame jaggedness so relax tolerance
    @test isapprox(pop_10, expected_10, rtol=2e-2)
    @test isapprox(pop_01, 1 .- expected_10, rtol=2e-2)
end

########################## ME Solver #################################
# This tests the decay of a driven two level system. Added to the hamiltonian is
# is a drive term of the form (J σ exp(iωt) + h.c.). When ω is chosen to
# coincide with the splitting of the two level system, J expressed in Hz is the
# rate of oscillation. For a system that only has T1 decay, the exponential rate
# of decay is Td = ⁴/₃ T1.

@testset "decaying driven two levels" begin
    dims = 3
    qubit_freq = 5.0
    anharm = -0.2
    q0 = DuffingTransmon("q0", dims, DuffingSpec(qubit_freq, anharm))
    cqs = CompositeQSystem([q0])
    add_hamiltonian!(cqs, q0)
    add_hamiltonian!(cqs, dipole_drive(q0, t -> 0.02*cos(2π*qubit_freq * t)), q0)

    T1 = 50. # in ns
    γ = 1/(2π * T1) # in GHz
    add_lindblad!(cqs, decay(q0, γ), [q0])

    ψ0 = ComplexF64[1; 0; 0]
    ρ0 = ψ0 * ψ0'
    times = collect(range(0, stop=100, length=101))
    # use floquet and me_propagator for speed
    prop = floquet_propagator(me_propagator, 1/qubit_freq)
    ρs = [reshape(u * vec(ρ0), dims, dims) for u in prop(cqs, times)]

    ρ00 = [real(ρ[1, 1]) for ρ in ρs]

    model_t2(x, p) = (exp.(-x ./ p[1]) .* cos.(2π .* p[2] .* x .- p[3]) .+ 1.0) / 2.0

    fit = curve_fit(model_t2, times, ρ00, [T1, .02, 0])
    @test abs(fit.param[1] - 4. * T1 / 3.) < 1.
end
