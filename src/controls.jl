# helper function to add time varyinig control Hamiltonian's

export microwave_drive,
       flux_drive

"""
    microwave_drive(q::QSystem, drive::Function)

Given a function `drive` that returns the drive amplitude at time `t` return a function applying a
microwave drive hamiltonian.
"""
function microwave_drive(q::QSystem, drive::Function)
    function add_drive_ham!(ham, idxs, t)
        pulse = 2π * drive(t)
        drive_ham = real(pulse) * X(q) + imag(pulse) * Y(q)
        QSimulator.expand_add!(ham, drive_ham, idxs)
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
        QSimulator.expand_add!(ham, tt_ham, idxs)
    end
end
