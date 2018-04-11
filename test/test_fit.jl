########################### Fixed Transmon Fit ##########################

f_01 = 5 # f_01 transition frequency
α = -0.2 # anharmonicity
dim = 3 # dimension

# Fit EC and EJ
EC_fit, EJ_fit = QSimulator.fit_fixed_transmon(f_01, α, dim)

# Verify that fit EC and EJ are consisent with f_01 and α
q0 = FixedTransmon("Test", EC_fit, EJ_fit, dim)
fit_levels = eigvals(hamiltonian(q0))
fit_f_01 = fit_levels[2]-fit_levels[1]
fit_f_12 = fit_levels[3]-fit_levels[2]
fit_α = fit_f_12 - fit_f_01
@test isapprox(f_01, fit_f_01; rtol=1e-3, atol=1e-3)
@test isapprox(α, fit_α; rtol=1e-3, atol=1e-3)

############################ Tunable Transmon Fit ##########################

f_01_max = 5 # f_01 transition frequency at 0 flux
f_01_min = 4 # f_01 transition frequency at 1/2 flux
α_max = -0.2 # anharmonicity at 0 flux
dim = 3 # dimension
model_fit = [TunableTransmon, TunableDuffingTransmon]

for model = model_fit
    # Fit EC and EJ
    EC_fit, EJ_fit, d_fit = QSimulator.fit_tunable_transmon(f_01_max, f_01_min, α_max, dim, model)

    # # Verify that fit EC and EJ are consisent with f_01 and α
    q1 = model("Test", EC_fit, EJ_fit, d_fit, dim)
    fit_levels_max = eigvals(hamiltonian(q1, 0))
    fit_levels_min = eigvals(hamiltonian(q1, 0.5))
    fit_f_01_max = fit_levels_max[2]-fit_levels_max[1]
    fit_f_12_max = fit_levels_max[3]-fit_levels_max[2]
    fit_f_01_min = fit_levels_min[2]-fit_levels_min[1]
    fit_α_max = fit_f_12_max - fit_f_01_max
    @test isapprox(f_01_max, fit_f_01_max; rtol=1e-3, atol=1e-3)
    @test isapprox(f_01_min, fit_f_01_min; rtol=1e-3, atol=1e-3)
    @test isapprox(α_max, fit_α_max; rtol=1e-3, atol=1e-3)
end
