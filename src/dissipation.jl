#using QIP

abstract Dissipation

label(d::Dissipation) = d.label
rate(d::Disspation) = d.rate

type Cooling <: Dissipation
  label::String
  rate::Real
  sys::QSystem
end

function generator(d::Cooling, t::Real)
  local a
  a = lowering(d.sys)
  return rate(d)*QIP.dissipator(a)
end

function generator(sys::CompositeQSystem, t::Real)
  local H, Gsop
  H = hamiltonian(sys, t)
  # [todo: blargh]
  Gsop = QIP.hamiltonian(H)

  for (ct, d) in enumerate(sys.dissipators)
    Gsop += generator(d)
  end
end

function 
end