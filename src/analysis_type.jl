""" IPBeam(;effct=:N_ALL, bcond=:CC, bcurv=:CONST, gamma=0.2, elnum= bcurv==:CONST ? 1 : 10,
            ro=2237e-9, kn=1.2, nu=0.3, el=1.06e3, eta=0.1, lam=100, phit=2pi/3, rad0=0.2)

Analysis type: In-plane dynamic analysis of beams.

Parameters
----------
::Int64         elnum       : Total number of elements
::Symbol        effct       : Type of analysis. Assumptions.
                                :N_ALL    all effects (nonlocal)
                                :L_ALL    all effects (local)
                                :N_SHEAR  shear deformation considered (nonlocal)
                                :N_AXIAL  axial extension considered (nonlocal)
                                :N_ROT    rotational inertia considered (nonlocal)
                                :N_NONE   no effect is considered
::Symbol       bcond       : Define boundary condition
                                :CC  clamped - clamped
                                :HH  hinged - hinged
                                :HC  hinged - clamped
                                :CF  clamped - free
                                :FC  free-clamped
                                :FF  free - free
::Symbol       bcurv        : Beam curvature
                                :VAR_PAR    Variable Curvature (parabola)
                                :CONST      Constant Curvature
::Float64      gamma        : Nonlocal scale parameter (constant)
::Float64      ro           : Weight of cross_section of the beam
::Int64        elnum        : Total number of elements
::Float64      kn           : Shear correction factor
::Float64      el           : Elastic modulus
::Float64      nu           : Poisson's ratio
::Float64      lam          : Slenderness ratio
::Float64      phit         : Opening angle (radian)
::Float64      rad0         : R0 radius_init
::Float64      eta          : A parameter for definition of cross-section height (0.1 to 0.4).
"""
type IPBeam

    effct::Symbol
    bcond::Symbol
    bcurv::Symbol
    gamma::Float64
    elnum::Int64
    ro::Float64
    kn::Float64
    nu::Float64
    el::Float64
    eta::Float64
    lam::Float64
    phit::Float64
    rad0::Float64

  function IPBeam(;effct=:N_ALL, bcond=:CC, bcurv=:CONST, gamma=0.2, elnum= bcurv==:CONST ? 1 : 10,
    ro=2237e-9, kn=1.2, nu=0.3, el=1.06e3, eta=0.1, lam=100, phit=2pi/3, rad0=0.2)

        new(effct, bcond, bcurv, gamma, elnum, ro, kn, nu, el, eta, lam, phit, rad0)

  end

end

""" OPBeam(;effct=:N_ALL, bcond=:CC, bcurv=:CONST, elnum= bcurv==:CONST ? 1 : 10,
            gamma=0.2, ro=2237e-9, nu=0.3, kb=1.2, el=1.06e3, eta=0.1, lam=100, phit=2pi/3, rad0=0.2)

Analysis type: Out-of-plane dynamic analysis of beams.

Parameters
----------
::Int64         elnum   : Total number of elements
::Symbol        effct   : Type of analysis. Assumptions.
                            :N_ALL       all effects (nonlocal)
                            :L_ALL       all effects (local)
                            :N_NOSHEAR   no shear. (nonlocal)
                            :N_ROT_BEND  rotatory inertia due to bending vibration considered (nonlocal)
                            :N_ROT_TORS  rotatory inertia due to torsional vibration considered (nonlocal)
                            :N_NONE      no effect is considered
::Symbol       bcond   : Define boundary condition
                            :CC  clamped - clamped
                            :HH  hinged - hinged
                            :HC  hinged - clamped
                            :CF  clamped - free
                            :FC  free - clamped
                            :FF  free - free
::Symbol       bcurv    : Beam curvature
                            :VAR    Variable Curvature
                            :CONST  Constant Curvature
::Int64        elnum    : Number of total elements
::Float64      gamma    : Nonlocal scale parameter
::Float64      ro       : Weight of cross_section of the beam
::Float64      nu       : Poisson's ratio
::Float64      kb       : Shear correction factor
::Float64      el       : Elastic modulus
::Float64      eta      : A parameter for definition of cross-section height (0.1 to 0.4).
::Float64      lam      : Slenderness ratio
::Float64      phit     : Opening angle (radian)
::Float64      rad0     : R0 radius_init
"""
type OPBeam

    effct::Symbol
    bcond::Symbol
    bcurv::Symbol
    elnum::Int64
    gamma::Float64
    ro::Float64
    nu::Float64
    kb::Float64
    el::Float64
    eta::Float64
    lam::Float64
    phit::Float64
    rad0::Float64

  function OPBeam(;effct=:N_ALL, bcond=:CC, bcurv=:CONST, elnum= bcurv==:CONST ? 1 : 10,
    gamma=0.2, ro=2237e-9, nu=0.3, kb=1.2, el=1.06e3, eta=0.1, lam=100, phit=2pi/3, rad0=0.2)

    new(effct, bcond, bcurv, elnum, gamma, ro, nu, kb, el, eta, lam, phit, rad0)

  end

end


typealias Beam Union{IPBeam, OPBeam}  # Analysis type = IPBeam or OPBeam type.
