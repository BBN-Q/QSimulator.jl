# For comparing benchmarks between two branches
# This file should be run from the top-level directory of QSimulator via compare_benchmarks.sh.
# This file benchmarks the "old" branch, for instance the branch with a performance regression.
# It saves the benchmarking parameters and the results so they can be used to compare
# with the "new" branch.

push!(LOAD_PATH, pwd())  # This allows using a copy of the repository.
include("benchmark.jl")
tune!(suite)
BenchmarkTools.save("benchmark_params.json", params(suite))
results = run(suite; verbose=true)
BenchmarkTools.save("old_benchmark.json", results)
show(results)
