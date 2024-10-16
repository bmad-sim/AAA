#=

Utility functions needed for tracking.

=#


"""
    sr_gamma(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the relativistic Lorentz factor ``\\gamma``.
"""
function sr_gamma(beta_gamma)
  return hypot(1, beta_gamma)
end


"""
    sr_gamma_m1(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the relativistic Lorentz factor minus one, ``\\gamma - 1``.
"""
function sr_gamma_m1(beta_gamma)
  return beta_gamma^2 / (hypot(1, beta_gamma) + 1)
end


"""
    sr_beta(beta_gamma)

For a particle with relativistic parameter ``\\beta\\gamma``,
compute the normalized velocity ``\\beta = v / c``.
"""
function sr_beta(beta_gamma)
  return beta_gamma / hypot(1, beta_gamma)
end


"""
    sr_pc(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the energy-equivalent momentum, ``pc``.
"""
function sr_pc(e_rest, beta_gamma)
  return e_rest * beta_gamma
end


"""
    sr_ekin(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the kinetic energy,
``E_\\text{kin} = mc^2(\\gamma - 1)``.
"""
function sr_ekin(e_rest, beta_gamma)
  return e_rest * sr_gamma_m1(beta_gamma)
end


"""
    sr_etot(e_rest, beta_gamma)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the total energy, ``E_\\tot = mc^2\\gamma``.
"""
function sr_etot(e_rest, beta_gamma)
  return e_rest * hypot(1, beta_gamma)
end


"""
    brho(e_rest, beta_gamma, ne = 1)

For a particle with a given rest energy and relativistic parameter
``\\beta\\gamma``, compute the magnetic rigidity, ``B\\rho = p / q``.

DTA: Need to handle energy units other than ``\\mathrm{eV}``..
"""
function brho(e_rest, beta_gamma, ne = 1)
  return sr_pc(e_rest, beta_gamma) / (ne * clight)
end


## If given ``E_\text{kin}`` instead of ``\beta\gamma``,
## use the following:
#
#function sr_gamma(e_rest, e_kin)
#  return e_kin / e_rest + 1
#end
#
#function sr_gamma_m1(e_rest, e_kin)
#  return e_kin / e_rest
#end
#
#function sr_beta_gamma(e_rest, e_kin)
#  return sqrt(e_kin / e_rest * (e_kin / e_rest + 2))
#end
#
#function sr_beta(e_rest, e_kin)
#  return sr_beta_gamma(e_rest, e_kin) / sr_gamma(e_rest, e_kin)
#end
#
#function sr_pc(e_rest, e_kin)
#  #return sqrt(e_kin * (e_kin + 2e_rest))
#  return e_rest * sr_beta_gamma(e_rest, e_kin)
#end
#
#function sr_etot(e_rest, e_kin)
#  return e_rest + e_kin
#end
#
#function brho(e_rest, e_kin, ne = 1)
#  return sr_pc(e_rest, e_kin) / (ne * clight)
#end

"""
    sincu(z)

Compute the unnormalized sinc function, ``\\operatorname{sincu}(z) = \\operatorname{sin}(z) / z``,
with a correct treatment of the apparent singularity at the origin.
"""
function sincu(z::Number)
  threshold = sqrt(2eps())
  if abs(z) < threshold
    return 1.
  else
    return sin(z) / z
end

"""
    sinch(z)

Compute the hyperbolic sinc function, ``\\operatorname{sinch}(z) = \\operatorname{sinh}(z) / z``,
with a correct treatment of the apparent singularity at the origin.
"""
function sinch(z::Number)
  threshold = sqrt(2eps())
  if abs(z) < threshold
    return 1.
  else
    return sinh(z) / z
end

