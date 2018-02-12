# QSimulator.jl

Package for simulating time dynamics of quantum systems with a focus on superconducting qubits.
Based off of [BBN QSimulator package](https://github.com/BBN-Q/QSimulator.jl)


## Benchmarks

We can track the code performance using the benchmarking suite in `test/benchmark.jl`.

```julia
julia> include(joinpath(Pkg.dir("QSimulator"), "test", "benchmark.jl"))

julia> results = run(suite; verbose=true)
(1/1) benchmarking "unitary"...
  (1/1) benchmarking "propagator"...
    (1/1) benchmarking "free evolution"...
      (1/3) benchmarking "2 level transmon"...
      done (took 2.286840473 seconds)
      (2/3) benchmarking "3 level transmon"...
      done (took 3.195876166 seconds)
      (3/3) benchmarking "4 level transmon"...
      done (took 5.418358854 seconds)
    done (took 11.267108865 seconds)
  done (took 11.637083596 seconds)
done (took 12.000488034 seconds)
1-element BenchmarkTools.BenchmarkGroup:
  tags: []
  "unitary" => 1-element BenchmarkGroup([])

julia> showall(results)
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
```
