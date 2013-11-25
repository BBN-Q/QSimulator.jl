#using QIP

#import QIP.hamiltonian

export ## Types
       Cooling,
       ## Methods
       label,
       rate,
       generator

label(d::Dissipation) = d.label
rate(d::Dissipation) = d.rate

type Cooling <: Dissipation
  label::String
  rate::Real
  sys::QSystem
end

function generator(d::Cooling, t::Real)
  local a, n, id, r
  a = lowering(d.sys)
  n = a'*a
  id = eye(Float64,size(a,1))
  r = 2pi*rate(d)
  # return pairs of left/right linear operators
  return [(sqrt(r)*a,sqrt(r)*a'),
          (-1/2*r*n,id),
          (id,-1/2*r*n)]
end

function generator(c::CompositeQSystem, t::Real=0.)
  local H, Gsop, id, d
  H = QSimulator.hamiltonian(c, t)
  d = size(H,1)
  id = eye(d)

  #Gsop = QIP.liou(-1im*H,id)
  #add!(Gsop, QIP.liou(id,1im*H))
  Gsop = QIP.hamiltonian(H)

  for (ct, diss) in enumerate(c.dissipators)
      for sop in generator(diss,t)
          add!(Gsop, 
               QIP.liou(expand(sop[1],c.subSystemExpansions[ct],d),
                        expand(sop[2],c.subSystemExpansions[ct],d)))
      end
  end
  return Gsop
end