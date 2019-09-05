# QSimulator.jl

[![Build Status](https://travis-ci.org/BBN-Q/QSimulator.jl.svg?branch=master)](https://travis-ci.org/BBN-Q/QSimulator.jl)

Package for simulating time dynamics of quantum systems with a focus on superconducting qubits.

## Installation

```
(v1.2) pkg> add https://github.com/BBN-Q/QSimulator.jl
```

## Unit tests

```julia
Pkg.test("QSimulator")
```

## Benchmarks

We can track the code performance between commits by running the benchmarking suite in
`benchmark/benchmarks.jl` using [PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl).

```julia
julia> import QSimulator

julia> using PkgBenchmark

julia> results = benchmarkpkg("QSimulator")
PkgBenchmark: Running benchmarks...
PkgBenchmark: using benchmark tuning data in /home/cryan/repos/QSimulator.jl/benchmark/tune.json
Benchmarking:  54%|███████████████████████████████████████████████████████████████▎                                                      |  ETA: 0:01:00
    [2/2]:        "unitary"
      [1/2]:      "propagator"
        [1/4]:    "lab frame parametric 2Q gate"
          [2/2]:  "2 transmons"
          ⋮

julia> export_markdown(stdout, results)
# Benchmark Report for *.*

## Job Properties
* Time of benchmark: 3 Sep 2019 - 22:11
* Package commit: dirty
* Julia commit: c6da87
* Julia command flags: None
* Environment variables: None

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                                | time            | GC time    | memory          | allocations |
|---------------------------------------------------|----------------:|-----------:|----------------:|------------:|
| `["operators", "X(q) (2 levels)"]`                |  97.304 ns (5%) |            |  272 bytes (1%) |           4 |
| `["operators", "X(q) (3 levels)"]`                | 102.005 ns (5%) |            |  320 bytes (1%) |           4 |
| `["operators", "X(q) (4 levels)"]`                | 108.733 ns (5%) |            |  384 bytes (1%) |           4 |
| `["operators", "X(q) (5 levels)"]`                | 113.097 ns (5%) |            |  464 bytes (1%) |           4 |
| `["operators", "X(q, 0.123) (2 levels)"]`         |  95.520 ns (5%) |            |  272 bytes (1%) |           4 |
| `["operators", "X(q, 0.123) (3 levels)"]`         | 101.485 ns (5%) |            |  320 bytes (1%) |           4 |
| `["operators", "X(q, 0.123) (4 levels)"]`         | 105.386 ns (5%) |            |  384 bytes (1%) |           4 |
| `["operators", "X(q, 0.123) (5 levels)"]`         | 116.289 ns (5%) |            |  464 bytes (1%) |           4 |
| `["operators", "lowering(q) (2 levels)"]`         |  77.935 ns (5%) |            |  240 bytes (1%) |           3 |
| `["operators", "lowering(q) (3 levels)"]`         |  80.843 ns (5%) |            |  288 bytes (1%) |           3 |
⋮
```

### Comparing benchmarks

To track performance regressions or improvements run the benchmark suite as above on the two
different branches.

```julia
julia> using Statistics: median

julia> judgement = judge("QSimulator", "faster-branch", "master"; f=median)

julia> export_markdown(stdout, judgement)
⋮
```
