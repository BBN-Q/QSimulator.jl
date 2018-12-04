#!/bin/bash
set -ex

# Switch to same directory that this script lives in
cd "$(dirname "$0")"

# Switch to project base directory
pushd ../

julia --project -e "using Pkg; Pkg.test()"
