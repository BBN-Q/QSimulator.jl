# Perturbative expansion of the Mathieu functions for the energy levels of a tunable transmon. The
# expansion is performed in the dimensionless paramter ξ = √(2EC/EJ). Each 3 element vector v below
# gives a term  (v[1]/2^v[2]) * ξ^v[3]
# For details see  Didier, N., Sete, E. A., da Silva, M. P., & Rigetti, C. (2017). Analytical
# modeling of parametrically-modulated transmon qubits. http://arxiv.org/abs/1706.06566

export perturbative_transmon_freq, perturbative_transmon_anharm
export perturbative_transmon_λ, perturbative_transmon_Λ
export transmon_N, transmon_φ

include("perturbative_transmon_ueigen.jl") # kind of messy!

"""
    transmon_φ(dim::Integer, ξ, ϕ::Real=0.0)
Transmon phase operator, φ̂. This function returns a matrix in the basis defined
by the Fock states of the linearized transmon Hamiltonian. See eq. (5) in
arXiv:1706.06566.

Tip: type φ with `\\varphi`.
"""
function transmon_φ(dim::Integer, ξ, ϕ::Real=0.0)
    return √(ξ) * (raising(dim, ϕ) + lowering(dim, ϕ))
end

"""
    transmon_N(dim::Integer, ξ, ϕ::Real=0.0)
Transmon charge number operator, N̂. This function returns a matrix in the basis
defined by the Fock states of the linearized transmon Hamiltonian. See eq. (5)
in arXiv:1706.06566.
"""
function transmon_N(dim::Integer, ξ, ϕ::Real=0.0)
    return (im/(2*√ξ)) * (raising(dim, ϕ) - lowering(dim, ϕ))
end

# pertubative expansion coefficients for the 0 ↔ 1 transisition
const PT_FREQ = [
    (4., 0, -1),
    (-1., 0, 0),
    (-1., 2, 1),
    (-21., 7, 2),
    (-19., 7, 3),
    (-5319., 15, 4),
    (-6649., 15, 5),
    (-1180581., 22, 6),
    (-446287., 20, 7),
    (-1489138635., 31, 8),
    (-648381403., 29, 9),
    (-614557854099., 38, 10),
    (-75265839129., 34, 11),
    (-637411859250147., 46, 12),
    (-86690561488017., 42, 13),
    (-405768570324517701., 53, 14),
    (-15191635582891041., 47, 15),
    (-2497063196283456607731., 63, 16),
    (-102281923716042917215., 57, 17),
    (-2292687293949773041433127., 70, 18),
    (-25544408245062216574759., 62, 19),
    (-4971071120163260007203175705., 78, 20),
    (-59956026877695226936825271., 70, 21),
    (-6299936888270974385982624367587., 85, 22),
    (-20465345194746565030172477629., 75, 23),
    (-36984324599399309412347250837528543., 94, 24),
    (-128862667153189778842334459173303., 84, 25),
    (-62313306363032484263243187605857135455., 101, 26),
    (-57979140436623262897403437875845329., 89, 27),
    (-239123145585215826671902664445701932931163., 109, 28),
    (-236766982175345008541229542667501202161., 97, 29),
    (-518667194120793070334115565427490753019904133., 116, 30)
]

# pertubative expansion coefficients for the anharmonicity
const PT_ANH = [
    (0., 0, -1),
    (-1., 0, 0),
    (-9., 4, 1),
    (-81., 7, 2),
    (-3645., 12, 3),
    (-46899., 15, 4),
    (-1329129., 19, 5),
    (-20321361., 22, 6),
    (-2648273373., 28, 7),
    (-45579861135., 31, 8),
    (-1647988255539., 35, 9),
    (-31160327412879., 38, 10),
    (-2457206583272505., 43, 11),
    (-50387904068904927., 46, 12),
    (-2145673984043982897., 50, 13),
    (-47368663010124907041., 53, 14),
    (-17329540083222030375645., 60, 15),
    (-410048712835835979799431., 63, 16),
    (-20066784213453521778111375., 67, 17),
    (-507447585299180759749453827., 70, 18),
    (-53019019946496461235728807475., 75, 19),
    (-1429754157181172012054040903645., 78, 20),
    (-79571741391885949104006842758911., 82, 21),
    (-2283773190022904454409743892590327., 85, 22),
    (-540565733415401595950277192471356985., 91, 23),
    (-16479511149218202447739080120870460083., 94, 24),
    (-1034743270413623494225962156243473940687., 98, 25),
    (-33436402163767100825528499521512789823595., 101, 26),
    (-4445866752212713247387096882061024305480817., 106, 27),
    (-151943275517394187333713519962683890787229463., 109, 28),
    (-10671890872277478133986830879693284605566580129., 113, 29),
    (-384886904723357410697832985058012690322303378753., 116, 30)
]

const PT_λ = [
    (1., 0, 0),
    (-1., 3, 1),
    (-11., 8, 2),
    (-65., 11, 3),
    (-4203., 17, 4),
    (-40721., 20, 5),
    (-1784885., 25, 6),
    (-21465147., 28, 7),
    (-4455462653., 35, 8),
    (-61698199851., 38, 9),
    (-3623317643901., 43, 10),
    (-56143119646191., 46, 11),
    (-7321743985484303., 52, 12),
    (-125280019793719221., 55, 13),
    (-8984438512815167237., 60, 14),
    (-168544684286400995331., 63, 15),
    (-105741913308715347076701., 71, 16),
    (-2164311753394257835891059., 74, 17),
    (-184798694135089048676718297., 79, 18),
    (-4109869091672376619457585371., 82, 19),
    (-761062061371895548979377743237., 88, 20),
    (-18317012159331390907042783219855., 91, 21),
    (-1831630981593132690479908285273395., 96, 22),
    (-47512263370928552970648689915451821., 99, 23),
    (-20440707519371829420653298425077482201., 106, 24),
    (-569157711742925565406447462105395143103., 109, 25)
]

const PT_Λ = [
    (1., 0, 0),
    (-1., 2, 1),
    (-73., 9, 2),
    (-79., 9, 3),
    (-113685., 19, 4),
    (-747533., 21, 5),
    (-175422349., 28, 6),
    (-698471247., 29, 7),
    (-1520876829389., 39, 8),
    (-13668058962903., 41, 9),
    (-4122722770459287., 48, 10),
    (-2534488707574995., 46, 11),
    (-26543348405245135937., 58, 12),
    (-281548290669062665101., 60, 13),
    (-98933257452818263360213., 67, 14),
    (-561603848629069641896937., 68, 15),
    (-3372037991404912212166296765., 79, 16),
    (-40819563311626093062783992331., 81, 17),
    (-16314102788878455728540034311379., 88, 18),
    (-52535388424912627194648863334467., 88, 19),
    (-178610931461508948221684711385383067., 98, 20),
    (-2444937960639526361173164055382471707., 100, 21),
    (-1103567409503040799217165335410059740779., 107, 22),
    (-8017554417550804194373089101907638666069., 108, 23),
    (-30711842188423912661533983529887505235301321., 118, 24),
    (-473069922042437374183190305740304564254754227., 120, 25)
]

_prepare(x) = [(par[3], par[1] / (2.0 ^ par[2])) for par in x]
const PT_FREQ_PRE = _prepare(PT_FREQ)
const PT_ANH_PRE = _prepare(PT_ANH)
const PT_λ_PRE = _prepare(PT_λ)
const PT_Λ_PRE = _prepare(PT_Λ)

const PERTURBATIVE_NUM_TERMS = 10 # default number of terms to use

"""
    xi_effective(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real)

Calculate the dimensionless transmon parameter ξ = √(2EC/EJ) where EJ is
the effective Josephson energy at some flux.

## args
* `EC`: the charging energy.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
"""
function xi_effective(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real)
    return sqrt((2 * EC) / sqrt((EJ1 ^ 2 + EJ2 ^ 2) + (2 * EJ1 * EJ2) * cos(2π*ϕ)))
end

function perturbative_sum_ξ(ξ::Real, PT::Vector{<:Tuple{<:Integer, <:Real}})
    series_term = par -> par[2] * (ξ ^ par[1])
    return mapreduce(series_term, +, PT)
end

"""
    perturbative_sum(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real,
                          PT::Vector{<:Tuple{<:Integer, <:Real}})

Compute the series expansion in the small parameter ξ.

## args
* `EC`: the charging energy.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `PT`: an array of two element tuples containing a power and a coefficient.

## returns
`EC * ∑ (coefficient * ξ^power)` for `(power, coefficient)` in `PT`.
"""
function perturbative_sum(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real,
                          PT::Vector{<:Tuple{<:Integer, <:Real}})
    ξ = xi_effective(EC, EJ1, EJ2, ϕ)
    return EC * perturbative_sum_ξ(ξ, PT)
end

"""
    perturbative_transmon_freq(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=$(PERTURBATIVE_NUM_TERMS))

Calculate the 01 transition frequency of a tunable or fixed transmon.

## args
* `EC`: the charging energy.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `num_terms`: the number of terms in the perturbative expansion to use.
    Must be at least 1 and at most $(length(PT_FREQ_PRE)).

## returns
The frequency of the transmon.
"""
function perturbative_transmon_freq end

"""
    perturbative_transmon_anharm(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=$(PERTURBATIVE_NUM_TERMS))

Calculate the anharmonicity a tunable or fixed transmon. The sign convention
is that the anharmonicity of a transmon is negative.

## args
* `EC`: the charging energy.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `num_terms`: the number of terms in the perturbative expansion to use.
    Must be at least 1 and at most $(length(PT_ANH_PRE)).

## returns
The anharmonicity of the transmon.
"""
function peturbative_transmon_anharm end

"""
    perturbative_transmon_λ(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=$(PERTURBATIVE_NUM_TERMS))

Calculate the weight of the number operator on the g-e σ_y operator in the
three-level transmon eigenbasis, i.e. `N̂ = λ σ_y / (2*√ξ) + ⋯`

## args
* `EC`: the charging energy.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `num_terms`: the number of terms in the perturbative expansion to use.
    Must be at least 1 and at most $(length(PT_λ_PRE)).

## returns
The value of λ for the transmon.
"""
function perturbative_transmon_λ end

"""
    perturbative_transmon_Λ(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=$(PERTURBATIVE_NUM_TERMS))

Calculate the weight of the number operator on the e-f s_y operator in the
three-level transmon eigenbasis, i.e. `N̂ = ⋯ + Λ s_y / √(2ξ) + ⋯`

## args
* `EC`: the charging energy.
* `EJ1`: the Josephson energy of one junction.
* `EJ2`: the Josephson energy of the other junction.
* `ϕ`: the external flux through the SQUID in units of the flux quantum.
* `num_terms`: the number of terms in the perturbative expansion to use.
    Must be at least 1 and at most $(length(PT_Λ_PRE)).

## returns
The value of Λ for the transmon.
"""
function perturbative_transmon_Λ end

for (name, array, scale) in [(:perturbative_transmon_freq,   PT_FREQ_PRE, true),
                             (:perturbative_transmon_anharm, PT_ANH_PRE,  true),
                             (:perturbative_transmon_λ,      PT_λ_PRE,   false),
                             (:perturbative_transmon_Λ,      PT_Λ_PRE,   false)]
    @eval function ($name)(EC::Real, EJ1::Real, EJ2::Real, ϕ::Real; num_terms::Int=PERTURBATIVE_NUM_TERMS)
        max_num_terms = length($array)
        num_terms < 1 && throw(ArgumentError("num_terms < 1"))
        (num_terms > max_num_terms) && @warn "only using $max_num_terms terms."
        t = min(num_terms, max_num_terms)
        ξ = xi_effective(EC, EJ1, EJ2, ϕ)
        return *(perturbative_sum_ξ(ξ, ($array)[1:t]), $(ifelse(scale, :((EC,)), ()))...)
    end
end

######################################################
# Backwards compatibility
######################################################
export mathieu_f01, mathieu_η

function mathieu_f01(t_params::Tuple{Float64, Float64, Float64}, ϕ::Real)
    @warn "Deprecation warning: mathieu_f01."
    EC, EJ1, EJ2 = t_params
    return perturbative_transmon_freq(EC, EJ1, EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
end

function mathieu_η(t_params::Tuple{Float64, Float64, Float64}, ϕ)
    @warn "Deprecation warning: mathieu_η."
    EC, EJ1, EJ2 = t_params
    return perturbative_transmon_anharm(EC, EJ1, EJ2, ϕ, num_terms=PERTURBATIVE_NUM_TERMS)
end
