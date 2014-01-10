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
