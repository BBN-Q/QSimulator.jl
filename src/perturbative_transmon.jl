# Perturbative expansion of the Mathieu functions for the energy levels of a tunable transmon. The
# expansion is performed in the dimensionless paramter ξ = √(2EC/EJ). Each 3 element vector v below
# gives a term  (v[1]/2^v[2]) * ξ^v[3]
# For details see  Didier, N., Sete, E. A., da Silva, M. P., & Rigetti, C. (2017). Analytical
# modeling of parametrically-modulated transmon qubits. http://arxiv.org/abs/1706.06566

export perturbative_transmon_freq, perturbative_transmon_anharm, PERTURBATIVE_NUM_TERMS

# pertubative expansion coefficients for the 0 ↔ 1 transisition
const PT_FREQ = [
    [4., 0, -1],
    [-1., 0, 0],
    [-1., 2, 1],
    [-21., 7, 2],
    [-19., 7, 3],
    [-5319., 15, 4],
    [-6649., 15, 5],
    [-1180581., 22, 6],
    [-446287., 20, 7],
    [-1489138635., 31, 8],
    [-648381403., 29, 9],
    [-614557854099., 38, 10],
    [-75265839129., 34, 11],
    [-637411859250147., 46, 12],
    [-86690561488017., 42, 13],
    [-405768570324517701., 53, 14],
    [-15191635582891041., 47, 15],
    [-2497063196283456607731., 63, 16],
    [-102281923716042917215., 57, 17],
    [-2292687293949773041433127., 70, 18],
    [-25544408245062216574759., 62, 19],
    [-4971071120163260007203175705., 78, 20],
    [-59956026877695226936825271., 70, 21],
    [-6299936888270974385982624367587., 85, 22],
    [-20465345194746565030172477629., 75, 23],
    [-36984324599399309412347250837528543., 94, 24],
    [-128862667153189778842334459173303., 84, 25],
    [-62313306363032484263243187605857135455., 101, 26],
    [-57979140436623262897403437875845329., 89, 27],
    [-239123145585215826671902664445701932931163., 109, 28],
    [-236766982175345008541229542667501202161., 97, 29],
    [-518667194120793070334115565427490753019904133., 116, 30]
]

# pertubative expansion coefficients for the anharmonicity
const PT_ANH = [
    [0., 0, -1],
    [-1., 0, 0],
    [-9., 4, 1],
    [-81., 7, 2],
    [-3645., 12, 3],
    [-46899., 15, 4],
    [-1329129., 19, 5],
    [-20321361., 22, 6],
    [-2648273373., 28, 7],
    [-45579861135., 31, 8],
    [-1647988255539., 35, 9],
    [-31160327412879., 38, 10],
    [-2457206583272505., 43, 11],
    [-50387904068904927., 46, 12],
    [-2145673984043982897., 50, 13],
    [-47368663010124907041., 53, 14],
    [-17329540083222030375645., 60, 15],
    [-410048712835835979799431., 63, 16],
    [-20066784213453521778111375., 67, 17],
    [-507447585299180759749453827., 70, 18],
    [-53019019946496461235728807475., 75, 19],
    [-1429754157181172012054040903645., 78, 20],
    [-79571741391885949104006842758911., 82, 21],
    [-2283773190022904454409743892590327., 85, 22],
    [-540565733415401595950277192471356985., 91, 23],
    [-16479511149218202447739080120870460083., 94, 24],
    [-1034743270413623494225962156243473940687., 98, 25],
    [-33436402163767100825528499521512789823595., 101, 26],
    [-4445866752212713247387096882061024305480817., 106, 27],
    [-151943275517394187333713519962683890787229463., 109, 28],
    [-10671890872277478133986830879693284605566580129., 113, 29],
    [-384886904723357410697832985058012690322303378753., 116, 30]
]

const PT_FREQ_PRE = [[par[3], par[1] / (2 ^ par[2])] for par in PT_FREQ]
const PT_ANH_PRE = [[par[3], par[1] / (2 ^ par[2])] for par in PT_ANH]
const MAX_NUM_TERMS = length(PT_FREQ)
const PERTURBATIVE_NUM_TERMS = 10 # default number of terms to use
const MAXDER = 20 # expansion for fourier coefficients

"""
    xi_effective(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real)

Calculate the dimensionless transmon paramter ξ = √(2EC/EJ) where EJ is
the effective Josephson energy at some flux.

## args
* `EC`: the charging enery.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
"""
function xi_effective(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real)
    return sqrt((2 * EC) / sqrt((EJ1 ^ 2 + EJ2 ^ 2) + (2 * EJ1 * EJ2) * cos(2π*ϕ)))
end

function perturbative_sum_ξ(ξ::Real, PT::Array{<:Array{<:Real, 1}, 1})
    series_term = par -> par[2] * (ξ ^ par[1])
    return mapreduce(series_term, +, PT)
end

"""
    perturbative_sum(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real,
                          PT::Array{<:Array{<:Real, 1}, 1})

Compute the series expansion in the small parameter ξ.

## args
* `EC`: the charging enery.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `PT`: an array of two element arrays containing a power and a coefficient.

## returns
The sum of `coefficient * ξ^power` for `[power, coefficient]` in `PT`.
"""
function perturbative_sum(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real,
                          PT::Array{<:Array{<:Real, 1}, 1})
    ξ = xi_effective(EC, EJ1, EJ2, ϕ)
    return EC * perturbative_sum_ξ(ξ, PT)
end

"""
    perturbative_transmon_freq(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=PERTURBATIVE_NUM_TERMS)

Calculate the 01 transition frequency of a tunable or fixed transmon.

## args
* `EC`: the charging enery.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `num_terms`: the number of terms in the perturbative expansion to use.
    Must be at least 1 and at most MAX_NUM_TERMS.

## returns
The frequency of the transmon.
"""
function perturbative_transmon_freq(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=PERTURBATIVE_NUM_TERMS)
    return perturbative_sum(EC, EJ1, EJ2, ϕ, PT_FREQ_PRE[1:num_terms])
end

"""
    perturbative_transmon_anharm(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=PERTURBATIVE_NUM_TERMS)

Calculate the anharmonicity a tunable or fixed transmon. The sign convention
is that the anharmonicity of a transmon is negative.

## args
* `EC`: the charging enery.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `num_terms`: the number of terms in the perturbative expansion to use.
    Must be at least 1 and at most MAX_NUM_TERMS.

## returns
The anharmonicity of the transmon.
"""
function perturbative_transmon_anharm(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=PERTURBATIVE_NUM_TERMS)
    return perturbative_sum(EC, EJ1, EJ2, ϕ, PT_ANH_PRE[1:num_terms])
end

######################################################
# Backwards compatibility
######################################################
export mathieu_f01, mathieu_η

function mathieu_f01(t_params::Tuple{Float64, Float64, Float64}, ϕ::Real)
    warn("Deprecation warning: mathieu_f01.")
    EC, EJ1, EJ2 = t_params
    return perturbative_transmon_freq(EC, EJ1, EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
end

function mathieu_η(t_params::Tuple{Float64, Float64, Float64}, ϕ)
    warn("Deprecation warning: mathieu_η.")
    EC, EJ1, EJ2 = t_params
    return perturbative_transmon_anharm(EC, EJ1, EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
end
