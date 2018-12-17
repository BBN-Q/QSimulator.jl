# QSimulator
Travis CI:

master: [![Build Status](https://travis-ci.org/BBN-Q/QSimulator.jl.svg?branch=master)](https://travis-ci.org/BBN-Q/QSimulator.jl)

v0.2: [![Build Status](https://travis-ci.org/BBN-Q/QSimulator.jl.svg?branch=v0.2)](https://travis-ci.org/BBN-Q/QSimulator.jl)

Unitary and Lindbladian evolution of quantum states in Julia.

In order to use it, run

```julia
Pkg.clone("git@github.com:BBN-Q/QSimulator.jl.git")
```

## Dependencies

QSimulator uses Julia 0.4 compatible syntax.  It also depends on one Julia package not yet registered with METADATA.jl:
* [ExpmV.jl](https://github.com/marcusps/ExpmV.jl)
