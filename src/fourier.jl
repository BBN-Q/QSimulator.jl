using QuadGK: quadgk

export FourierSeries, eval_series, rotating_frame_series

const OFFSET = 2.0 # a good default offset for fourier_coefficient
# both for calculating FourierSeries for frequency under modulation
# and for calculating renormalization couplings

"""
    fourier_coefficient(f::Function, frequency::Real, harmonic::Int, offset::Number=OFFSET)

Calculate the fourier coefficient
`1/T ∫_0^T f(τ) e^(-2πi * frequency * harmonic * τ) dτ`
of a periodic function `f`. Optionally add a constant to the integrand so
that the integral is nonzero (this helps with convergence).

# args
* `f`: the function of one variable.
* `frequency`: the frequency of periodicity of `f`.
* `harmonic`: the desired harmonic.
* `offset`: the optional offset of the integral.

# returns
The fourier coefficient.
"""
function fourier_coefficient(f::Function, frequency::Real, harmonic::Int, offset::Number=OFFSET)
    T = 1/frequency
    g(τ) = offset/T + f(τ) * exp(-2π * 1im * frequency * harmonic * τ)
    integral, _ = quadgk(g, 0.0, T)
    return (integral - offset)/T
end

"""
Representation of a Fourier series `f(t) = ∑ A_n e^(iωnt)`
"""
struct FourierSeries
    frequency::Real
    terms::Dict{Int, <:Number} # harmonic => coefficient
end

"""
    FourierSeries(f::Function, frequency::Real, harmonics::Vector{Int}, offsets::Union{Vector{<:Number}, Nothing}=nothing)

Construct a FourierSeries from a periodic function `f` on given harmonics.

# args
* `f`: the periodic function.
* `frequency`: the frequency of periodicity of f.
* `harmonics`: a vector of harmonics to compute.
* `offsets`: a corresponding ovector of offsets to use in `fourier_coefficient`. Defaults to all OFFSET.

# returns
A FourierSeries.
"""
function FourierSeries(f::Function, frequency::Real, harmonics::Vector{Int}, offsets::Union{Vector{<:Number}, Nothing}=nothing)
    if offsets == nothing; offsets = OFFSET * ones(length(harmonics)); end
    terms = Dict(h => fourier_coefficient(f, frequency, h, o) for (h, o) in zip(harmonics, offsets))
    return FourierSeries(frequency, terms)
end

"""
    eval_series(fs::FourierSeries, times::Vector{<:Real})

Evaluate a FourierSeries on given times.

# args
* `fs`: a FourierSeries.
* `times`: a vector of times.

# returns
A vector of values.
"""
function eval_series(fs::FourierSeries, times::Vector{<:Real})
    v = collect(fs.terms)
    harmonic_vector, coefficient_vector = [p[1] for p in v], [p[2] for p in v]
    return exp.(2π * 1im * fs.frequency * times * transpose(harmonic_vector)) * coefficient_vector
end

"""
    rotating_frame_series(fs::FourierSeries, harmonics::Vector{Int})

Compute the fourier series for `exp(i∫_0^t f(τ) dτ)` where `f` is
the function given by the FourierSeries `fs`. This function ignores a
constant term in `fs` since that maps onto an overall exponential factor
which is not part of a FourierSeries.

# args
* `fs`: a FourierSeries for the frequency of a qubit.
* `harmonics`: the desired harmonics.

# returns
The FourierSeries.
"""
function rotating_frame_series(fs::FourierSeries, harmonics::Vector{Int})
    iω = 2π * 1im * fs.frequency
    function exp_of_integral(t::Real)
        int = 0.0
        for (h, c) in fs.terms
            if h != 0
                int += c * (exp(iω * h * t) - 1)/(iω * h)
            end
        end
        return exp(1im * int)
    end
    # note the integrals are between -1 and 1 so adding 2 will suffice
    return FourierSeries(exp_of_integral, fs.frequency, harmonics, OFFSET * ones(length(harmonics)))
end
