export ## Types
       Field,
       MicrowaveControl,
       QuadratureControl,
       ## Methods
       amplitude,
       label,
       load_sequence!

#AWG channels are controls

#A double-balanced mixer driven by an RF/microwave souce and a single AWG channel
mutable struct MicrowaveControl <: Control
    label::AbstractString
    freq::Float64
    phase::Float64
    timeStep::Float64
    sequence::AbstractInterpolation
end
MicrowaveControl(label, freq; phase=0., timeStep=1/1.2, sequence=interpolate([0.0, 0.0], BSpline(Interpolations.Constant()))) =
    MicrowaveControl(label, freq, phase, timeStep, sequence)

# Methods to create the interpolated sequence object
function load_sequence!(mc::MicrowaveControl, sequence::Vector; interpolation=BSpline(Interpolations.Constant()))
    mc.sequence = interpolate(sequence, interpolation)
end

function load_sequence!(mc::MicrowaveControl, seqDict::Dict, n::Int; interpolation=BSpline(Interpolations.Constant()))
    #For double balanced mixer need one AWG channel
    #Hack around off the expected BBNAPSx-xx style
    chan = label(mc)[end]
    APSLabel = label(mc)[1:7]
    mc.sequence = interpolate(seqDict[APSLabel]["chan_"*string(chan)][n], interpolation)
    return nothing
end
function amplitude(mc::MicrowaveControl, t::Float64)
    #Scale time by the timestep before interpolating
    #Note that Julia indexing starting at 1
    if t < 1.0
        return mc.sequence[t/mc.timeStep + 1]*cos(2*pi*mc.freq*t + mc.phase)
    else
        return mc.sequence[t/mc.timeStep]*cos(2*pi*mc.freq*t + mc.phase)
    end
end

#A pair of AWG channels driving an IQ mixer with a microwave source at a given frequency
struct QuadratureControl <: Control
    label::AbstractString
    freq::Float64
    phase::Float64
    timeStep::Float64
    sequence_I::AbstractInterpolation
    sequence_Q::AbstractInterpolation
end
QuadratureControl(label, freq; phase=0., timeStep=1/1.2,
        sequence_I=interpolate( [0.0], BSpline(Constant()), OnGrid() ),
        sequence_Q=interpolate( [0.0], BSpline(Constant()), OnGrid() )) =
    QuadratureControl(label, freq, phase, timeStep, sequence_I, sequence_Q)

#Create the interpolated object
function load_sequence!(qc::QuadratureControl, seqDict::Dict, n::Int; interpolation=BSpline(Constant()))
    #For quadrature controls we need two channels
    #Hack around off the expected BBNAPSx-xx style
    chan_I = label(qc)[end-1]
    chan_Q = label(qc)[end]
    APSLabel = label(qc)[1:7]
    qc.sequence_I = interpolate(seqDict[APSLabel]["chan_"*string(chan_I)][n], interpolation, OnGrid())
    qc.sequence_Q = interpolate(seqDict[APSLabel]["chan_"*string(chan_Q)][n], interpolation, OnGrid())
    return nothing
end
function amplitude(qc::QuadratureControl, t::Float64)
    #Create the complex I/Q pair
    #Scale time by the timestep
    phasor = qc.sequence_I[t/qc.timeStep] - 1im*qc.sequence_Q[t/qc.timeStep]
    return abs(phasor)*cos(2*pi*qc.freq*t + qc.phase + angle(phasor))
end

#AWG controls get connected, possibly through a transfer function, to fields
struct Field
    control::Control
end
#Pull in amplitude from controls
#[todo] - allow for transfer function
amplitude(f::Field, t::Float64) = amplitude(f.control, t)
frequency(f::Field) = f.control.freq
