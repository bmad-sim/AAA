module Symplectic
using ..BeamTracking: Coords, Beam
using ..GTPSA: @FastGTPSA!

export track!


struct Drift{T}
  L::T  # drift length / m
end

struct Quadrupole{T}
  L::T   # quadrupole length / m
  B1::T  # quadrupole gradient / (T·m^-1)
end


"""
    track!(ele::Symplectic.Drift, beamf::Beam, beami::Beam) -> beamf

## Description
This function performs exact, hence symplectic, tracking of a `Beam` through a drift.

### Arguments
  - `ele`   -- element of type `Symplectic.Drift`
  - `beamf` -- final `Beam`
  - `beami` -- initial `Beam`
"""
function track!(ele::Drift, beamf::Beam, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  @assert !(beamf.species === beami.species) "Input species must be equal to output!"
  L = ele.L
  zi = beami.z
  zf = beamf.z

  tilde_m = mass(beami.species) / sr_pc(mass(beami.species), beami.beta_gamma_0)

  @. ps = sqrt((1.0 + zi.pz)^2 - zi.px^2 - zi.py^2)
  @. et = sqrt((1.0 + zi.pz)^2 + tilde_m^2)
  @. zf.x = zi.x + zi.px * L / ps
  @. zf.px = zi.px
  @. zf.y = zi.y + zi.py * L / ps
  @. zf.py = zi.py
  @. zf.z = zi.z - (1.0 + zi.pz) * (L / ps - L / (beami.beta_0 * et))
  @. zf.pz = zi.pz

  return beamf
end # function track!(::Drift, ::Beam, ::Beam)


"""
    track!(ele::Quadrupole, beamf::Beam, beami::Beam)

## Description
track quadrupole


### Arguments
  - `ele`   -- element of type `Symplectic.Quadrupole`
  - `beamf` -- final `Beam`
  - `beami` -- initial `Beam`

### Implementation
This integrator uses the so-called Matrix-Kick-Matrix method to construct
an integrator accurate though second-order in the integration step-size.
"""
function track!(ele::Symplectic.Quadrupole, beamf::Beam, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  L = ele.L

  k2_num = charge(species) * ele.B1 / ele.P0

  trackQuadMx!(beamf, beami, k2_num, L / 2)
  trackQuadK!( beami, beamf, L)
  trackQuadMx!(beamf, beami, k2_num, L / 2)

  return beamf
end # function track!(::Quadrupole)


"""
    trackQuadMx!(beamf::Beam, beami::Beam, k2_num::Float64, s::Float64)

track "matrix part" of quadrupole
"""
function trackQuadMx!(beamf::Beam, beami::Beam, k2_num::Float64, s::Float64)
  zi = beami.z
  zf = beamf.z:w


  @. p_red = 1 + zi.pz  # reduced momentum, P/P0 = 1 + δ
  @. xp = px / p_red
  @. yp = py / p_red
  @. k2 = k2_num / p_red
  @. ks = sqrt(abs(k2)) * s

  @. cx =  (k2_num >= 0) * cos(ks)    + (k2_num < 0) * cosh(ks)
  @. cy =  (k2_num >= 0) * cosh(ks)   + (k2_num < 0) * cos(ks)
  @. sx =  (k2_num >= 0) * sincu(ks)  + (k2_num < 0) * sinch(ks)
  @. sy =  (k2_num >= 0) * sinch(ks)  + (k2_num < 0) * sincu(ks)
  @. sx2 = (k2_num >= 0) * sincu(2ks) + (k2_num < 0) * sinch(2ks)
  @. sy2 = (k2_num >= 0) * sinch(2ks) + (k2_num < 0) * sincu(2ks)
  @. sxz = (k2_num >= 0) * sin(ks)^2  + (k2_num < 0) * sinh(ks)^2
  @. syz = (k2_num >= 0) * sinh(ks)^2 + (k2_num < 0) * sin(ks)^2

  @. zf.x  = zi.x  * cx + xp * s * sx
  @. zf.px = zi.px * cx - k2 * p_red * zi.x * s * sx
  @. zf.y  = zi.y  * cy + yp * s * sy
  @. zf.py = zi.py * cy + k2 * p_red * zi.y * s * sy
  @. zf.z  = (zi.z - (s / 4) * (xp^2 * (1 + sx2) + yp^2 * (1 + sy2)
                                + k2 * zi.x^2 * (1 - sx2) - k2 * zi.y^2 * (1 - sy2))
                   + (zi.x * xp * sxz - zi.y * yp * syz) / 2.)
  @. zf.pz = zi.pz

  return beamf
end # function trackQuadMx


"""
    trackQuadK!(beamf::Beam, beami::Beam, s::Float64)

track "remaining part" of quadrupole, a position kick
"""
function trackQuadK!(beamf::Beam, beami::Beam, s::Float64)
  zi = beami.z
  zf = beamf.z

  tilde_m = mass(beami.species) / sr_pc(mass(beami.species), beami.beta_gamma_0)
  β0 = sr_beta(mass(beami.species), beami.beta_gamma_0)

  @. ps = sqrt((1.0 + zi.pz)^2 - zi.px^2 - zi.py^2)
  @. et = sqrt((1.0 + zi.pz)^2 + tilde_m^2)
  @. p_red = 1 + zi.pz  # reduced momentum, P/P0 = 1 + δ
  @. xp = px / p_red
  @. yp = py / p_red

  @. zf.x  = zi.x + s * zi.x * (1 / (sqrt(1 - (xp^2 + yp^2)) - 1))
  @. zf.y  = zi.y + s * zi.y * (1 / (sqrt(1 - (xp^2 + yp^2)) - 1))
  @. zf.z  = zi.z - s * zi.y * (1 / (sqrt(1 - (xp^2 + yp^2)) - 1)
                                - (xp^2 + yp^2) / 2
                                - 1 / (β0 * sqrt(1 + (tilde_m / p_red)^2)))
  @. zf.px = zi.px
  @. zf.py = zi.py
  @. zf.pz = zi.pz

  return beamf
end # function trackQ!M::Quadrupole()


end # module Symplectic

