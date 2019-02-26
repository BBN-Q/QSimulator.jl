# QSimulator.jl

Package for simulating time dynamics of quantum systems with a focus on superconducting qubits.
Based off of [BBN QSimulator package](https://github.com/BBN-Q/QSimulator.jl)

## Installation

```
Pkg.clone("ssh://git@bitbucket.lab.rigetti.com:7999/qos/qsimulator.jl.git", "QSimulator")
```

## Unit tests

```julia
Pkg.test("QSimulator")
```

## Benchmarks

We can track the code performance using the benchmarking suite in `test/benchmark.jl`.

```julia
julia> import QSimulator

julia> include(joinpath(dirname("QSimulator"), "test", "benchmark.jl"))

julia> tune!(suite)
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "unitary" => 2-element BenchmarkTools.BenchmarkGroup:
	  tags: []
	  "propagator" => 4-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "lab frame parametric 2Q gate" => 2-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "3 transmons" => Benchmark(evals=1, seconds=5.0, samples=10000)
			  "2 transmons" => Benchmark(evals=1, seconds=5.0, samples=10000)
			  ⋮

julia> results = run(suite; verbose=true)
(1/1) benchmarking "unitary"...
  (1/1) benchmarking "propagator"...
    (1/1) benchmarking "free evolution"...
      (1/3) benchmarking "2 level transmon"...
      done (took 2.286840473 seconds)
      ⋮

julia> show(results)
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "unitary" => 1-element BenchmarkTools.BenchmarkGroup:
	  tags: []
	  "propagator" => 1-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "free evolution" => 3-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "2 level transmon" => Trial(138.139 μs)
			  "3 level transmon" => Trial(208.912 μs)
			  "4 level transmon" => Trial(391.713 μs)
        ⋮
```

### Comparing benchmarks

To track performance regressions or improvements run the benchmark suite as above on the two
different branches. Save the results from one session and then load them into another to compare.

Save tuning results and benchmark results
```julia
julia> include(joinpath(dirname("QSimulator"), "test", "benchmark.jl"))

julia> tune!(suite)
⋮
julia> BenchmarkTools.save("test/benchmark_params.json", params(suite))

julia> results = run(suite; verbose=true)
⋮
julia> BenchmarkTools.save("test/master_20190226.json", results)
```

Then change branches and compare.
```julia
julia> include(joinpath(dirname("QSimulator"), "test", "benchmark.jl"))

julia> loadparams!(suite, BenchmarkTools.load("params.json")[1], :evals, :samples);

julia> results = run(suite; verbose=true)
⋮
julia> master_results = BenchmarkTools.load("test/master_benchmark.json")[1]

julia> using Statistics: median

julia> show(judge(median(results), median(master_results)))
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "unitary" => 2-element BenchmarkTools.BenchmarkGroup:
	  tags: []
	  "propagator" => 4-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "lab frame parametric 2Q gate" => 2-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "3 transmons" => TrialJudgement(+30.75% => regression)
			  "2 transmons" => TrialJudgement(+141.74% => regression)
		  "free evolution" => 6-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "single transmon (2 levels)" => TrialJudgement(-49.67% => improvement)
			  "dipole chain (3 transmons)" => TrialJudgement(-57.51% => improvement)
			  "dipole chain (2 transmons)" => TrialJudgement(-61.38% => improvement)
			  "single transmon (4 levels)" => TrialJudgement(-44.73% => improvement)
			  "dipole chain (4 transmons)" => TrialJudgement(-85.14% => improvement)
			  "single transmon (3 levels)" => TrialJudgement(-55.76% => improvement)
		  "rabi flops" => 3-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "2 level transmon" => TrialJudgement(+358.12% => regression)
			  "3 level transmon" => TrialJudgement(+348.04% => regression)
			  "4 level transmon" => TrialJudgement(+276.86% => regression)
		  "rotating frame parametric 2Q gate" => 3-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "3 transmons" => TrialJudgement(+45.21% => regression)
			  "4 transmons" => TrialJudgement(-48.81% => improvement)
			  "2 transmons" => TrialJudgement(+313.31% => regression)
	  "pure state" => 4-element BenchmarkTools.BenchmarkGroup:
		  tags: []
		  "lab frame parametric 2Q gate" => 2-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "3 transmons" => TrialJudgement(+112.82% => regression)
			  "2 transmons" => TrialJudgement(+322.79% => regression)
		  "free evolution" => 6-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "single transmon (2 levels)" => TrialJudgement(-38.65% => improvement)
			  "dipole chain (3 transmons)" => TrialJudgement(-39.28% => improvement)
			  "dipole chain (2 transmons)" => TrialJudgement(-44.91% => improvement)
			  "single transmon (4 levels)" => TrialJudgement(-38.68% => improvement)
			  "dipole chain (4 transmons)" => TrialJudgement(+12.23% => regression)
			  "single transmon (3 levels)" => TrialJudgement(-49.57% => improvement)
		  "rabi flops" => 3-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "2 level transmon" => TrialJudgement(+388.11% => regression)
			  "3 level transmon" => TrialJudgement(+398.11% => regression)
			  "4 level transmon" => TrialJudgement(+427.58% => regression)
		  "rotating frame parametric 2Q gate" => 3-element BenchmarkTools.BenchmarkGroup:
			  tags: []
			  "3 transmons" => TrialJudgement(+267.74% => regression)
			  "4 transmons" => TrialJudgement(+40.88% => regression)
			  "2 transmons" => TrialJudgement(+679.59% => regression)
```
