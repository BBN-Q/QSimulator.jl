# For comparing benchmarks between two branches
# This file is meant to be run from compare_benchmarks.sh from the top-level directory of QSimulator.
# This file benchmarks the "new" branch, for instance the branch that has fixed a performance regression.
# It runs the benchmark on the new branch, then loads the results from the old branch and compares.

push!(LOAD_PATH, pwd())  # This allows using a copy of the repository.
include("benchmark.jl")
tune!(suite)
loadparams!(suite, BenchmarkTools.load("benchmark_params.json")[1], :evals, :samples);
results = run(suite; verbose=true)
old_results = BenchmarkTools.load("old_benchmark.json")[1]
using Statistics: median
show(judge(median(results), median(old_results)))
