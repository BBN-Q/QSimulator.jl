#!/bin/bash
set -ex

# Switch to same directory that this script lives in
cd "$(dirname "$0")"

# Switch to project base directory
pushd ../

# Create a Julia depot on the worker, if none exists
# This acts as a global cache for Julia dependencies
# The declarative nature of Pkg in Julia 1.0 should mean
# that this is a safe optimization.
mkdir -p ~/.julia

# Build base with system deps
docker build -f Dockerfile -t julia-base .

# Build deps. Cache on host
docker run -v ~/.julia:/opt/.julia -v $(pwd):/project-mount julia-base julia --project setup.jl

# Run tests
docker run -v ~/.julia:/opt/.julia -v $(pwd):/project-mount julia-base /project-mount/scripts/runtests.sh
