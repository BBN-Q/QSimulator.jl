export hamiltonian

""" The Hamiltonian of a resonator """
hamiltonian(r::Resonator) = r.frequency * number(r)

""" Calculate the drift or natural Hamiltonian of a CompositeQSystem """
function hamiltonian(cqs::CompositeQSystem)
    ham = zeros(Complex128, (dim(cqs), dim(cqs)))
    for (new_ham, idxs) = cqs.fixed_Hs
        embed_add!(ham, new_ham, idxs)
    end
    return ham
end

""" Tunable Transmon Hamiltonian in the charge basis """
function hamiltonian(t::TunableTransmon, flux)
  N = floor(Int, dim(t)/2)
  scaled_EJ = scale_EJ(t.E_J, flux, t.d)
  4 * t.E_C * diagm((-N:N).^2) - scaled_EJ  * 0.5 * (diagm(ones(dim(t)-1),-1) + diagm(ones(dim(t)-1),1))
end

""" Fixed Transmon Hamiltonian in the charge basis """
function hamiltonian(t::FixedTransmon)
  N = floor(Int, dim(t)/2)
  4 * t.E_C * diagm((-N:N).^2) - 0.5*t.E_J*(diagm(ones(dim(t)-1),-1) + diagm(ones(dim(t)-1),1))
end

""" Dispatch on FixedTransmon for the Hamiltonian of a Fixed Transmon """
hamiltonian(t::FixedDuffingTransmon) = (t.frequency - 0.5*t.anharmonicity)*number(t) + 0.5*t.anharmonicity * number(t)^2

function hamiltonian(t::TunableDuffingTransmon, flux)
    scaled_EJ = scale_EJ(t.E_J, flux, t.d)
    ωₚ = sqrt(8*t.E_C*scaled_EJ)
    return diagm([(ωₚ - t.E_C / 2) * ct - t.E_C / 2 * ct^2 for ct in 0:(t.dim-1)])
end

function hamiltonian(t::MathieuTransmon, flux)
    EJ₁ = .5 * (t.d + 1) * t.E_J
    EJ₂ = t.E_J - EJ1
    t_params = (t.E_C, EJ₁, EJ₂)
    f01 = mathieu_f01(t_params, flux)
    η = mathieu_η(t_params, flux)
    return diagm([0., f01, 2. * f01 - η])
end
