using Test, QSimulator
using SpecialFunctions: besselj
import PyPlot
const plt = PyPlot
plt.ioff()

@testset "unitary propagator time independent" begin
    # compare exponentiation, full unitary_propagator
    # and unitary_propagator with projection
    mat_type = Array{ComplexF64,2}
    m1 = convert(mat_type, [[0,1] [1,0]])
    m2 = convert(mat_type, [[0,1im,0] [-1im,.1,2+1im] [0,2-1im,0]])
    m3 = convert(mat_type, [3][:,:]) # no simple syntax for creating a 1x1 array
    mat = cat(m1, m2, m3, dims=[1,2])
    inds1 = [1,2]
    inds2 = [3,4,5]
    inds3 = [6]
    @test mat[inds1, inds1] == m1
    @test mat[inds2, inds2] == m2
    @test mat[inds3, inds3] == m3
    lh = LiteralHermitian("label", HermitianSpec(mat))
    cqs = CompositeQSystem([lh])
    add_hamiltonian!(cqs, lh)
    times = collect(range(0, stop=10, length=100))
    us_full = unitary_propagator(cqs, times)
    us_expm = [exp(-2π *1im * mat * t) for t in times]
    @test isapprox(us_expm[end], us_full[end]; atol=1e-5)
end

# iSWAP interaction
freq_1 = 4.0
anharm = -0.2
dims = [3,3]
q0 = DuffingTransmon("q0", dims[1], DuffingSpec(0.0, anharm)) # tunable qubit
q1 = DuffingTransmon("q1", dims[2], DuffingSpec(freq_1, anharm))
cqs = CompositeQSystem([q0, q1])
drive_freq = 0.3
a = freq_1 + drive_freq
b = drive_freq
add_hamiltonian!(cqs, parametric_drive(q0, t -> a + b * cos(2π * drive_freq * t)), q0)
add_hamiltonian!(cqs, q1)
g = 0.01
add_hamiltonian!(cqs, .5 * g * X_Y([q0, q1]), [q0, q1])
g_eff = g * besselj(1, b/drive_freq)
t_final = 1/(4 * g_eff)

@testset "unitary state time dependent" begin
    # compare state evolution to closed form RWA calculation
    make_plot = "plot" in ARGS
    num_times = make_plot ? 200 : 25
    times = collect(range(0, stop=t_final, length=num_times))
    ψ0 = photons_to_state([1,0], dims)
    ψs = unitary_state(cqs, times, ψ0)
    pop_10 = [abs2(ψ[photons_to_index([1,0], dims)]) for ψ in ψs]
    pop_01 = [abs2(ψ[photons_to_index([0,1], dims)]) for ψ in ψs]
    pop_10_RWA = cos.(2π*g_eff*times).^2
    pop_01_RWA = 1 .-pop_10_RWA
    if make_plot
        periods = collect(1:floor(times[end] * drive_freq))./drive_freq
        plt.plot(times, pop_10, label="10")
        plt.plot(times, pop_01, label="01")
        plt.plot(times, pop_10_RWA, label="10 RWA")
        plt.plot(times, pop_01_RWA, label="01 RWA")
        for period in periods
            plt.axvline(period, linestyle="--", color="grey")
        end
        plt.show()
    end
    @test isapprox(pop_10, pop_10_RWA, rtol=.03)
end


@testset "unitary propagator time dependent" begin
    times = collect(range(0, stop=t_final, length=3))
    us_full = unitary_propagator(cqs, times)
    # compare state evolution to propagator evolution
    ψ0 = photons_to_state([1,0], dims)
    ψs = unitary_state(cqs, times, ψ0)
    pop_10_state = [abs2(ψ[photons_to_index([1,0], dims)]) for ψ in ψs]

    pop_10_prop = [abs2(ψ0' * u * ψ0) for u in us_full]
    @test isapprox(pop_10_state, pop_10_prop, rtol=1e-5)
end

@testset "me_propagator" begin
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

    times = collect(range(0,stop=100,length=101))
    ψ0 = ComplexF64[1; 0; 0]
    ρ0 = ψ0 * ψ0'
    us_prop = me_propagator(cqs, times)
    ρs_prop = [reshape(u * vec(ρ0), dims, dims) for u in us_prop]
    ρs = me_state(cqs, times, ρ0)

    for i in 1:length(times)
        @test isapprox(ρs_prop[i], ρs[i], rtol=1e-5)
    end
    if "plot" in ARGS
        plt.plot(times, [real(ρ[1, 1]) for ρ in ρs], label="state")
        plt.plot(times, [real(ρ[1, 1]) for ρ in ρs_prop], label="propagator")
        plt.legend(loc="best")
        plt.show()
    end
end

@testset "decay and dephasing" begin
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

    times = collect(range(0,stop=100,length=101))
    ψ0 = ComplexF64[1; 1; 0]/sqrt(2)
    ρ0 = ψ0 * ψ0'
    ρs = me_state(cqs, times, ρ0)

    analytical_T1 = exp.(-times/T1)/2
    analytical_T2 = exp.(-times/T2)/2
    for i in 1:length(times)
        @test isapprox(ρs[i][2,2], analytical_T1[i], rtol=1e-10)
        @test isapprox(ρs[i][1,2], analytical_T2[i], rtol=1e-4)
    end

    if "plot" in ARGS
        plt.plot(times, [real(ρ[2, 2]) for ρ in ρs], label="excited population")
        plt.plot(times, analytical_T1, "--", label="T1")
        plt.plot(times, [real(ρ[1, 2]) for ρ in ρs], label="coherence")
        plt.plot(times, analytical_T2, "--", label="T2")
        plt.legend(loc="best")
        plt.show()
    end
end
