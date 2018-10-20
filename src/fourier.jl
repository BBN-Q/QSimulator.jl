using QuadGK: quadgk

export FourierSeries, eval_series, rotating_frame_series

"""
    fourier_coefficient(f::Function, frequency::Real, harmonic::Int)

Calculate the fourier coefficient `1/T ∫_0^T f(τ) e^(-2πi * frequency * harmonic * τ) dτ`
of a periodic function `f`.

# args
* `f`: the function of one variable.
* `frequency`: the frequency of periodicity of `f`.
* `harmonic`: the desired harmonic.

# returns
The fourier coefficient.
"""
function fourier_coefficient(f::Function, frequency::Real, harmonic::Int)
    T = 1/frequency
    g(t) = f(t) * exp(-2π * 1im * frequency * harmonic * t)
    integral, _ = quadgk(g, 0.0, T, atol=sqrt(eps(T)))
    return integral/T
end

"""
Representation of a Fourier series `f(t) = ∑ Aₙ e^(iωnt)`
"""
struct FourierSeries
    frequency::Real
    terms::Dict{Int, <:Number} # harmonic => coefficient
end

"""
    FourierSeries(f::Function, frequency::Real, harmonics::AbstractVector{Int})

Construct a FourierSeries from a periodic function `f` on given harmonics.

# args
* `f`: the periodic function.
* `frequency`: the frequency of periodicity of f.
* `harmonics`: a vector of harmonics to compute.

# returns
A FourierSeries.
"""
function FourierSeries(f::Function, frequency::Real, harmonics::AbstractVector{Int})
    terms = Dict(h => fourier_coefficient(f, frequency, h) for h in harmonics)
    return FourierSeries(frequency, terms)
end

"""
    eval_series(fs::FourierSeries, times::AbstractVector{<:Real})

Evaluate a FourierSeries on given times.

# args
* `fs`: a FourierSeries.
* `times`: a vector of times.

# returns
A vector of values.
"""
function eval_series(fs::FourierSeries, times::AbstractVector{<:Real})
    return reduce(+, coeff * exp.(2π * 1im * fs.frequency * harmonic * times) for (harmonic, coeff) in fs.terms)
end

"""
    rotating_frame_series(fs::FourierSeries, harmonics::AbstractVector{Int})

Compute the fourier series for `exp(i∫_0^t f(τ) dτ)` where `f` is the function given by
the FourierSeries `fs`. This function ignores a constant term in `fs` since that maps onto
an overall exponential factor which is not part of a FourierSeries.

# args
* `fs`: a FourierSeries.
* `harmonics`: the desired harmonics.

# returns
The FourierSeries.
"""
function rotating_frame_series(fs::FourierSeries, harmonics::AbstractVector{Int})
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
    return FourierSeries(exp_of_integral, fs.frequency, harmonics)
end
