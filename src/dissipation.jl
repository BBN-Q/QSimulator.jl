export ## Types
       Cooling,
       ## Methods
       label,
       rate,
       liouvillian_left,
       liouvillian_right,
       liouvillian_bilat

label(d::Dissipation) = d.label
rate(d::Dissipation) = d.rate

# NOTE: Only single body dissipators are currently supported

struct Cooling <: Dissipation
  label::AbstractString
  rate::Real
  system::QSystem
end

function liouvillian_left(d::Cooling, t::Real)
  # Note these return the DUAL to the time propagator generator
  local a, n, r
  a = lowering(d.system)
  n = a'*a
  r = rate(d)
  return -1/2*r*n
end

function liouvillian_right(d::Cooling, t::Real)
  # Note these return the DUAL to the time propagator generator
  local a, n, r
  a = lowering(d.system)
  n = a'*a
  r = rate(d)
  return -1/2*r*n
end

function liouvillian_bilat(d::Cooling, t::Real)
  # Note these return the DUAL to the time propagator generator
  local a, n, r
  a = lowering(d.system)
  r = rate(d)
  return r*kron(transpose(a),a')
end
