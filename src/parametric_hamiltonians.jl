# helper function to add time varying control Hamiltonian's

export microwave_drive,
       flux_drive,
       rotating_flip_flop

"""
    microwave_drive(q::QSystem, drive::Function)

Given a function `drive` that returns the drive amplitude at time `t` return a function applying a
microwave drive hamiltonian.
"""
function microwave_drive(q::QSystem, drive::Function)
    const x_ham = X(q)
    # TODO: specialize on whether control is complex to avoid branch below
    const y_ham = Y(q)
    function add_drive_ham!(ham, idxs, t)
        pulse = 2π * drive(t)
        drive_ham = real(pulse) * x_ham + imag(pulse) * y_ham
        QSimulator.embed_add!(ham, drive_ham, idxs)
    end
end


"""
    flux_drive(tt::QSystem, drive::Function)

Given a function for the flux (in units of Φ₀) return a function applyling the Hamiltonian for a
tunable qubit.
"""
function flux_drive(tt::QSystem, drive::Function)
    function add_drive_ham!(ham, idxs, t)
        tt_ham = 2π * hamiltonian(tt, drive(t))
        QSimulator.embed_add!(ham, tt_ham, idxs)
    end
end


"""
    rotating_flip_flop(a::QSystem, b::QSystem, strength, freq)

FlipFlop or XY interaction rotating in an interaction frame of the difference of number operators on
two QSystems.
"""
function rotating_flip_flop(a::QSystem, b::QSystem, strength, freq)
    # find non-zero elements of lower and upper half
    lower_half_op = raising(a) ⊗ lowering(b)
    dims = size(lower_half_op)
    nz = findnz(lower_half_op)
    # convert to linear indices
    const lower_triangle_nz = ([sub2ind(dims, r, c) for (r,c) = zip(nz[1], nz[2])], 2π*strength*nz[3])

    upper_half_op = lowering(a) ⊗ raising(b)
    nz = findnz(upper_half_op)
    # convert to linear indices
    const upper_triangle_nz = ([sub2ind(dims, r, c) for (r,c) = zip(nz[1], nz[2])], 2π*strength*nz[3])

    function add_flip_flop_ham!(ham, idxs, t)
        phase = exp(1im*2π*freq*t)
        # TODO: explore @inbounds
        for (idx, val) = zip(lower_triangle_nz...)
            for idxbis = idxs[idx]
                ham[idxbis] += phase*val
            end
        end
        conj_phase = conj(phase)
        for (idx, val) = zip(upper_triangle_nz...)
            for idxbis = idxs[idx]
                ham[idxbis] += conj_phase*val
            end
        end
    end
end
