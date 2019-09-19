using Test, QSimulator
using LinearAlgebra: diag, diagm, I
const pnt = QSimulator.PERTURBATIVE_NUM_TERMS
const cnt = QSimulator.CHARGE_NUM_TERMS

@testset "DuffingTransmon" begin
	freq = 4.0
	anharm = -0.2
	dims = [3]
	q = DuffingTransmon("q", dims[1], DuffingSpec(freq, anharm))
	h = hamiltonian(q)
	@test isapprox(h, diagm(0 => [0.0, freq, 2 * freq + anharm]))
end

@testset "PerturbativeTransmon" begin
	fmax = 4.0
	anharm_max = -0.2
	fmin = 3.0
	t = fit_transmon(fmax, fmin, anharm_max, PerturbativeTransmon, pnt)
	q = PerturbativeTransmon("q", 3, t)
	h = hamiltonian(q, 0.0)
	@test isapprox(h, diagm(0 => [0.0, fmax, 2 * fmax + anharm_max]))

	ξ = QSimulator.xi_effective(2,160,0,0)
	U = perturbative_transmon_Ueigen(ξ)
	@test size(U) == (5,25)
	Neigen = U*transmon_N(25, ξ)*U'

	σy = 1im*(raising(2) - lowering(2))
	λ = perturbative_transmon_λ(2,160,0,0; num_terms=6)
	Λ = perturbative_transmon_Λ(2,160,0,0; num_terms=6)

	# confirm eq. (24) in arXiv:1706.06566v2
	@test isapprox(Neigen[1:2,1:2], λ*σy/(2*√ξ); rtol=1e-6)
	@test isapprox(Neigen[2:3,2:3], Λ*σy/√(2ξ); rtol=1e-5)
end

@testset "PerturbativeTransmon vs DiagonalChargeBasisTransmon" begin
	t = TransmonSpec(.172, 12.71, 3.69)
	pt = PerturbativeTransmon("q", 3, t, num_terms=pnt)
	ct = DiagonalChargeBasisTransmon("q", 3, t, num_terms=cnt)
	for ϕ in range(0, stop=1.0, length=10)
		ham_pt = hamiltonian(pt, ϕ)
		ham_ct = hamiltonian(ct, ϕ)
	    @test isapprox(ham_pt + ham_ct[1,1] * I, ham_ct, rtol=1e-5)
	end
end

@testset "fit_transmon fixed" begin
	freq = 5.0
	anharm = -0.2
	ham = QSimulator.duffing_hamiltonian(freq, anharm, 3)

	function test_spec(s, model)
		@test s.EJ2 == 0.0
		ham_1 = hamiltonian(model("Test", 3, s))
		@test isapprox(ham + ham_1[1,1] * I, ham_1, rtol=1e-7)
	end

	fit_spec = fit_transmon(freq, freq, anharm, PerturbativeTransmon, pnt)
	test_spec(fit_spec, PerturbativeTransmon)
	test_spec(fit_spec, DiagonalChargeBasisTransmon)
	fit_spec = fit_transmon(freq, freq, anharm, DiagonalChargeBasisTransmon, cnt)
	test_spec(fit_spec, PerturbativeTransmon)
	test_spec(fit_spec, DiagonalChargeBasisTransmon)
end

@testset "fit_transmon tunable" begin
	# tunable transmon
	freq_max = 5.0
	freq_min = 4.0
	anharm_max = -0.2
	ham_max = QSimulator.duffing_hamiltonian(freq_max, anharm_max, 3)
	ham_min = QSimulator.duffing_hamiltonian(freq_max, 0.0, 2)


	function test_spec(s, model)
		@test s.EJ2 != 0.0
		ham_max_1 = hamiltonian(PerturbativeTransmon("Test", 3, s))
		@test isapprox(ham_max + ham_max_1[1,1] * I, ham_max_1, rtol=1e-7)

		ham_min_1 = hamiltonian(PerturbativeTransmon("Test", 2, s))
		@test isapprox(ham_min + ham_min_1[1,1] * I, ham_min_1, rtol=1e-8)
	end

	fit_spec = fit_transmon(freq_max, freq_min, anharm_max, PerturbativeTransmon, pnt)
	test_spec(fit_spec, PerturbativeTransmon)
	test_spec(fit_spec, DiagonalChargeBasisTransmon)
	fit_spec = fit_transmon(freq_max, freq_min, anharm_max, DiagonalChargeBasisTransmon, cnt)
	test_spec(fit_spec, PerturbativeTransmon)
	test_spec(fit_spec, DiagonalChargeBasisTransmon)
end

@testset "DuffingSpec <-> TransmonSpec" begin
	freq = 4.0
	anharm = -0.2
	d = DuffingSpec(freq, anharm)
	t = TransmonSpec(d)
	d1 = DuffingSpec(t)
	@test isapprox(d1.frequency, freq, rtol=1e-8)
	@test isapprox(d1.anharmonicity, anharm, rtol=1e-6)
end
