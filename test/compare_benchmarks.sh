#!/bin/env sh

# Run benchmarks on two branches and compare the results

# This file is meant to be run from the top-level directory of QSimulator.
# Benchmarking may take several minutes. So, you may want to work from
# a copy of the repository.
# LOAD_PATH is altered to look first in the current directory for modules

# Edit these branch names manually
old_branch=more-benchmarks
new_branch=parametric-types

git checkout ${old_branch}
julia test/benchmark_old_branch.jl

git checkout ${new_branch}
julia test/benchmark_new_branch.jl

