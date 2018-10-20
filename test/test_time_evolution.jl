using Test, QSimulator
<<<<<<< HEAD
using QSimulator: index
using SpecialFunctions: besselj
=======
using SpecialFunctions: besselj, erf
>>>>>>> af31f09... added floquet_propagator and updated test_time_evolution.jl
import PyPlot
const plt = PyPlot
plt.ioff()

@testset "unitary propagator time independent" begin
    # For a constant Hamiltonian compare exponentiation to full unitary_propagator
    # first, create a block diagonal Hamiltonian matrix
    m1 = [[0,1] [1,0]]
    m2 = [[0,1im,0] [-1im,.1,2+1im] [0,2-1im,0]]
    m3 = [3]
    mat = cat(m1, m2, m3, dims=[1,2])

    # Make a CompositeQSystem from the Hamiltonian
    lh = LiteralHermitian("label", HermitianSpec(mat))
    cqs = CompositeQSystem([lh])
    add_hamiltonian!(cqs, lh)
    times = range(0, stop=10, length=100)

    # compare the unitary propagator to matrix exponentiation
    us_prop = unitary_propagator(cqs, times)
    us_expm = [exp(-2π *1im * mat * t) for t in times]
    @test all(isapprox(u_expm, u_prop; atol=1e-4) for (u_expm, u_prop) in zip(us_expm, us_prop))
end

# parameters for an iSWAP interaction with linear frequency modulation
freq_fixed = 4.0
anharm = -0.2
drive_freq = 0.3

dims = (3,3)
basis = TensorProductBasis(dims)
ψ₀₁ = TensorProductBasisState(basis, (0,1))
ψ₁₀ = TensorProductBasisState(basis, (1,0))

# start with the tunable transmon at zero frquency and add frequency term with parametric_drive
q0 = DuffingTransmon("q0", dims[1], DuffingSpec(0.0, anharm))
q1 = DuffingTransmon("q1", dims[2], DuffingSpec(freq_fixed, anharm))
cqs = CompositeQSystem([q0, q1])
add_hamiltonian!(cqs, q1)

# park the tunable `drive_freq` above the tunable and linearly modulate its frequency with an
# amplitude equal to the modulation frequency
a = freq_fixed + drive_freq
b = drive_freq
add_hamiltonian!(cqs, parametric_drive(q0, t -> a + b * cos(2π * drive_freq * t)), q0)
g = 0.01
add_hamiltonian!(cqs, .5 * g * X_Y([q0, q1]), [q0, q1])
g_eff = g * besselj(1, b/drive_freq)
t_final = 1/(4 * g_eff)

@testset "unitary state time dependent" begin
    # compare state evolution to closed form RWA calculation
    make_plot = "plot" in ARGS
    num_times = make_plot ? 200 : 25
    times = range(0, stop=t_final, length=num_times)
    ψs = unitary_state(cqs, times, vec(ψ₁₀))
    pop₁₀ = [abs2(ψ[index(ψ₁₀)]) for ψ in ψs]
    pop₀₁ = [abs2(ψ[index(ψ₀₁)]) for ψ in ψs]
    pop₁₀_RWA = cos.(2π*g_eff*times).^2
    pop₀₁_RWA = 1 .- pop₁₀_RWA
    if make_plot
        plt.plot(times, pop₁₀, label="10")
        plt.plot(times, pop₀₁, label="01")
        plt.plot(times, pop₁₀_RWA, label="10 RWA")
        plt.plot(times, pop₀₁_RWA, label="01 RWA")
        # mark the period of the drive with vertical dashed lines
        for period in (1:floor(times[end] * drive_freq))./drive_freq
            plt.axvline(period, linestyle="--", color="grey")
        end
        plt.legend()
        plt.show()
    end
    @test isapprox(pop₁₀, pop₁₀_RWA, rtol=.03)
    @test isapprox(pop₀₁, pop₀₁_RWA, rtol=.03)
end


@testset "unitary propagator time dependent" begin
    # compare consistency of state evolution and propagator evolution
    times = range(0, stop=t_final, length=3)
    us_full = unitary_propagator(cqs, times)
    ψ₀ = TensorProductBasisState(basis, (1,0))
    ψs = unitary_state(cqs, times, vec(ψ₀))
    pop_10_state = [abs2(ψ[index(ψ₀)]) for ψ in ψs]
    pop_10_prop = [abs2(vec(ψ₀)' * u * vec(ψ₀)) for u in us_full]
    @test isapprox(pop_10_state, pop_10_prop, rtol=1e-4)
    # test helper function for saving at a point
    us_end = unitary_propagator(cqs, times[end])
    @test us_end == us_full[end]
    ψ_end = unitary_state(cqs, times[end], vec(ψ₀))
    @test ψ_end == ψs[end]
end

@testset "choose_times_floquet and decompose_times_floquet" begin
    t_period = 1/drive_freq
    dt = t_final/10
    times = choose_times_floquet(t_final/2, t_final, t_period, dt)
    @test times[1] > 0
    @test times[end] <= t_final
    @test t_final/2 in times
    prepend!(times, 0.0)
    @test times[1] == 0
    quotients, unique_remainders, unique_inds = decompose_times_floquet(times, t_period)
    @test length(quotients) == length(times)
    @test length(unique_remainders) < length(times)
    @test length(unique_inds) == length(times)
    @test isapprox(t_period * quotients + unique_remainders[unique_inds], times)
end

@testset "floquet with pulseshape" begin
    cqs = CompositeQSystem([q0, q1])
    drive_freq = 0.3
    a = freq_1 + drive_freq
    b = drive_freq
    g = 0.01
    g_eff = g * besselj(1, b/drive_freq)
    t_final = 1/(4 * g_eff)
    add_hamiltonian!(cqs, q1)
    add_hamiltonian!(cqs, .5 * g * X_Y([q0, q1]), [q0, q1])
    # create pulseshape
    rise_time = .1 * t_final
    fwhm = 0.25 * rise_time
    t1 = fwhm
    t2 = t_final - fwhm
    σ = 0.5 * fwhm / sqrt(2*log(2))
    function p(t)
        envelope = 0.5 * (erf((t - t1)/σ) - erf((t - t2)/σ)) # erf_squared
        drive = a + b * cos(2π * drive_freq * t)
        return envelope * drive + (1 - envelope) * (a+b)
    end
    times = collect(linspace(0.0, t_final, 500))
    if "plot" in ARGS
        plt.plot(times, p.(times))
        plt.axvline(rise_time, color="red")
        plt.axvline(t_final-rise_time, color="red")
        plt.show()
    end
    add_hamiltonian!(cqs, flux_drive(q0, p), q0)
    t_period = 1/drive_freq
    prop = floquet_propagator(unitary_propagator, t_period, rise_time, rise_time)
    us_floquet_pulseshape = prop(cqs, times)
    us_unitary = unitary_propagator(cqs, times)
    if "plot" in ARGS
        ind_10 = photons_to_index([1,0], dims)
        ind_01 = photons_to_index([0,1], dims)
        plt.plot(times, [abs2(u[ind_10, ind_10]) for u in us_unitary], color="lightcoral", label="10 unitary")
        plt.plot(times, [abs2(u[ind_01, ind_10]) for u in us_unitary], color="lightblue", label="01 unitary")
        plt.plot(times, [abs2(u[ind_10, ind_10]) for u in us_floquet_pulseshape], color="red", label="10 floquet")
        plt.plot(times, [abs2(u[ind_01, ind_10]) for u in us_floquet_pulseshape], color="blue", label="01 floquet")
        plt.legend(loc="best")
        plt.show()
    end
    for i in 1:length(times)
        @test isapprox(us_unitary[i], us_floquet_pulseshape[i], rtol=1e-5)
    end
end

@testset "me_propagator" begin
    # compare master equation state and propagator evolution for consistency with on-resonance CW
    # microwave
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

    times = range(0,stop=100,length=101)
    ψ0 = ComplexF64[1; 0; 0]
    ρ0 = ψ0 * ψ0'
    us_prop = me_propagator(cqs, times)
    ρs_prop = [reshape(u * vec(ρ0), dims, dims) for u in us_prop]
    ρs = me_state(cqs, times, ρ0)
    prop = floquet_propagator(me_propagator, 1/qubit_freq)
    us_floquet = prop(cqs, times)
    ρs_floquet = [reshape(u * vec(ρ0), dims, dims) for u in us_floquet]

<<<<<<< HEAD
    @test all(isapprox(ρ_prop, ρ_state, rtol=1e-4) for (ρ_prop, ρ_state) in zip(ρs_prop, ρs))

=======
    for i in 1:length(times)
        @test isapprox(ρs_prop[i], ρs[i], rtol=1e-5)
        @test isapprox(ρs_floquet[i], ρs[i], rtol=1e-5)
    end
>>>>>>> af31f09... added floquet_propagator and updated test_time_evolution.jl
    if "plot" in ARGS
        plt.plot(times, [real(ρ[1, 1]) for ρ in ρs], label="state")
        plt.plot(times, [real(ρ[1, 1]) for ρ in ρs_prop], label="propagator")
        plt.plot(times, [real(ρ[1, 1]) for ρ in ρs_floquet], label="floquet")
        plt.legend(loc="best")
        plt.show()
    end
    # test helper function for saving at a point
    us_end = me_propagator(cqs, times[end])
    @test us_end == us_prop[end]
    ρ_end = me_state(cqs, times[end], ρ0)
    @test ρ_end == ρs[end]
end

@testset "decay and dephasing" begin
    # compare master equation free evolution decay under T1 and T2 to exponential predictions
    dims = 3
    qubit_freq = 5.0
    anharm = -0.2
    q0 = DuffingTransmon("q0", dims, DuffingSpec(qubit_freq, anharm))
    cqs = CompositeQSystem([q0])
    add_hamiltonian!(cqs, q0)
    T1 = 50.0 # in ns
    γ1 = 1/(2π * T1) # in GHz
    Tϕ = 200.0 # ns
    γϕ = 1/(2π * Tϕ) # in GHz
    T2 = 1/(.5/T1 + 1/Tϕ)
    add_lindblad!(cqs, decay(q0, γ1), [q0])
    add_lindblad!(cqs, dephasing(q0, γϕ), [q0])

    times = range(0, stop=100, length=101)
    ψ0 = ComplexF64[1; 1; 0]/sqrt(2)
    ρ0 = ψ0 * ψ0'
    ρs = me_state(cqs, times, ρ0)

    analytical_T1 = exp.(-times/T1)/2
    analytical_T2 = exp.(-times/T2)/2

    # check diagonal and off-diagonal against exponential decay with T1/T2 rates
    @test all(isapprox(ρ[2,2], check_T1, rtol=1e-10) for (ρ,check_T1) in zip(ρs, analytical_T1))
    @test all(isapprox(ρ[1,2], check_T2, rtol=1e-2) for (ρ,check_T2) in zip(ρs, analytical_T2))

    if "plot" in ARGS
        plt.plot(times, [real(ρ[2, 2]) for ρ in ρs], label="excited population")
        plt.plot(times, analytical_T1, "--", label="exponential T₁")
        plt.plot(times, [real(ρ[1, 2]) for ρ in ρs], label="coherence")
        plt.plot(times, analytical_T2, "--", label="exponential T₂")
        plt.legend(loc="best")
        plt.show()
    end
end
