using LinearAlgebra: eigvals

########################### Fixed Transmon Fit ##########################

f_01 = 5 # f_01 transition frequency
α = -0.2 # anharmonicity
d = 101 # dimensionality of system
# Fit EC and EJ
fit_EC, fit_EJ = QSimulator.fit_fixed_transmon(f_01, α, d)

# Verify that fit EC and EJ are consisent with f_01 and α
q0 = FixedTransmon("Test", fit_EC, fit_EJ, d)
fit_levels = eigvals(hamiltonian(q0))
fit_f_01 = fit_levels[2]-fit_levels[1]
fit_f_12 = fit_levels[3]-fit_levels[2]
fit_α = fit_f_12 - fit_f_01
@test isapprox(f_01, fit_f_01; rtol=1e-3, atol=1e-3)
@test isapprox(α, fit_α; rtol=1e-3, atol=1e-3)

# Test against analytical formula for  eigenenergies (see Eq. 2.11 in 10.1103/PhysRevA.76.042319)
Em(m, EC, EJ) = -EJ + sqrt(8 * EC * EJ) * (m + 0.5) - (EC / 12) * (6 * m^2 + 6 * m + 3)

analytic_levels = [Em(m, fit_EC, fit_EJ) for m=0:2]
f01_analytic = analytic_levels[2] - analytic_levels[1]
α_analytic = analytic_levels[3] - 2 * analytic_levels[2] +analytic_levels[1]
@test isapprox(f_01, f01_analytic; rtol=1e-2)
@test isapprox(α, α_analytic; rtol=1e-3, atol=2e-2)

############################ Tunable Transmon Fit ##########################

f_01_max = 5 # f_01 transition frequency at 0 flux
f_01_min = 4 # f_01 transition frequency at 1/2 flux
α_max = -0.2 # anharmonicity at 0 flux
# dims = [101, 3]
# fit_model = [TunableTransmon, TunableDuffingTransmon]

ds = [3]
fit_model = [PerturbativeTransmon]

for (model, d) = zip(fit_model, ds)
    # Fit EC and EJ
    fit_EC, fit_EJ, fit_d = QSimulator.fit_tunable_transmon(f_01_max, f_01_min, α_max, d, model)

    fit_EJ1, fit_EJ2 = QSimulator.asymmetry_to_EJs(fit_EJ, fit_d)

    # Verify that fit EC and EJ are consisent with f_01 and α
    q1 = model("Test", d, TransmonSpec(fit_EC, fit_EJ1, fit_EJ2))
    fit_levels_max = eigvals(hamiltonian(q1, 0))
    fit_levels_min = eigvals(hamiltonian(q1, 0.5))
    fit_f_01_max = fit_levels_max[2]-fit_levels_max[1]
    fit_f_12_max = fit_levels_max[3]-fit_levels_max[2]
    fit_f_01_min = fit_levels_min[2]-fit_levels_min[1]
    fit_α_max = fit_f_12_max - fit_f_01_max
    @test isapprox(f_01_max, fit_f_01_max; rtol=1e-3, atol=1e-3)
    @test isapprox(f_01_min, fit_f_01_min; rtol=1e-3, atol=1e-3)
    @test isapprox(α_max, fit_α_max; rtol=1e-3, atol=1e-3)

    fit_EJ_extremas = [QSimulator.scale_EJ(fit_EJ, phi, fit_d) for phi = [0, 0.5]]
    fit_levels_extremas = [fit_levels_max[1:3], fit_levels_min[1:3]]

    # Test against analytical formula at min and max (see Eq. 2.11 in 10.1103/PhysRevA.76.042319)
    for (fit_EJ_extrema, fit_levels_extrema)= zip(fit_EJ_extremas, fit_levels_extremas)
        analytic_levels = [Em(m, fit_EC, fit_EJ_extrema) for m=0:2][1:3]
        @test isapprox(analytic_levels[2]-analytic_levels[1],
                       fit_levels_extrema[2]-fit_levels_extrema[1];
                       rtol=1e-3, atol=3e-2)
        @test isapprox(analytic_levels[3]-analytic_levels[2],
                       fit_levels_extrema[3]-fit_levels_extrema[2];
                       rtol=1e-2)
    end
end
