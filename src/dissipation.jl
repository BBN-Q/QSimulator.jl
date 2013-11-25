#using QIP

export Dissipation,
       label,
       rate,
       Cooling,
       generator

abstract Dissipation

label(d::Dissipation) = d.label
rate(d::Disspation) = d.rate

type Cooling <: Dissipation
  label::String
  rate::Real
  sys::QSystem
end

function generator(d::Cooling, t::Real)
  local a, n, i
  a = lowering(d.sys)
  n = a'*a
  # return pairs of left/right linear operators
  return [(sqrt(rate)*a,sqrt(rate)*a'),
          (-1/2*rate(d)*n,id),
          (id,-1/2*rate(d)*n)]
end

function generator(sys::CompositeQSystem, t::Real=0.)
  local H, Gsop, id, d
  H = hamiltonian(sys, t)
  d = size(H,1)
  id = eye(d)

  Gsop = QIP.liou(-1im*H,id)
  add!(Gsop, QIP.liou(id,1im*H))

  for (ct, diss) in enumerate(sys.dissipators)
      for sop in generator(diss,t)
          add!(Gsop, 
               QIP.liou(expand(sop[1],c.subSystemExpansions[ct],d),
                        expand(sop[2],c.subSystemExpansions[ct],d)))
      end
  end
  return Gsop
end