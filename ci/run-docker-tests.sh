#!/bin/bash
set -ex

# Switch to same directory that this script lives in
cd "$(dirname "$0")"

# Switch to project base directory
pushd ../

docker build -f Dockerfile -t qsimulator-tests .

docker run qsimulator-tests
