""" element_fundamental(atype::IPBeam, element_no::Int64, om::Number)

Calculate fundamental matrix for each element of the beam.

Parameters
----------
::IPBeam    aytpe       : In-plane dynamic analysis
::Int64     element_no  : Element no.
::Float64   om          : Frequency

Returns
-------
A(6x6) fundamental matrix for a given element.
"""
function element_fundamental(atype::IPBeam, element_no::Int64, om::Number)

    effct = atype.effct
    bcurv = atype.bcurv
    thetat = atype.thetat
    rad0 = atype.rad0
    lam = atype.lam
    gam = atype.gamma
    eta = atype.eta
    kn = atype.kn
    ro = atype.ro
    el = atype.el
    nu = atype.nu


    h0 = sqrt(12) * rad0 * thetat / lam  # Height_init for cross-section (to be used in height(θ) function)
    b = h0 * 2/3  # Width of the cross-section (constant)
    shear = el / (2*(1+nu))  # Shear modulus

    radius(θ) = bcurv == :CONST ? rad0 : rad0 / (cos(θ))^3  # Gives radius for given angle.
    height(θ) = θ <= 0 ? h0 * (1 - eta * θ/thetat) : h0 * (1 + eta * θ/thetat);
    # Gives height of the cross-section for given angle.

    amat = zeros(6, 6);  # Fundamental matrix initialized.

    rad = (radius(coords(atype, element_no)) + radius(coords(atype, element_no + 1))) / 2  # Radius of this element
    h = (height(coords(atype, element_no)) + height(coords(atype, element_no + 1))) / 2  # Height of this element
    area = b * h  # Cross-sectional area of this element.
    mu = ro * area  # Weight of unit length beam for this particular element.
    ib = (b * h^3) / 12  # Torsional area moment of inertia for this element.

    if effct == :N_ALL  # All effects nonlocal.
        amat[1, 2] = (1 + rad / (el * area) * (gam/rad)^2 * rad * mu * om^2) /
                     (1 - rad / (el * area) * (gam/rad)^2 * rad * mu * om^2)
        amat[1, 5] = (rad / (el * area) * (1 + (gam/rad)^2)) /
                     (1 - rad / (el * area) * (gam/rad)^2 * rad * mu * om^2)
        amat[2, 1] = -(1 + rad * kn / (shear * area)*(gam/rad)^2*rad*mu*om^2) /
                     (1 - rad * kn / (shear * area) * (gam/rad)^2*rad*mu*om^2)
        amat[2, 3] = rad / (1 - rad * kn/(shear*area)*(gam/rad)^2*rad*mu*om^2)
        amat[2, 6] = (rad * kn / (shear * area) * (1+(gam/rad) ^ 2)) /
                     (1 - rad * kn / (shear * area) * (gam/rad)^2*rad*mu*om^2)
        amat[3, 2] = -(rad / (el * ib) * (gam/rad)^2 * rad^2*mu * om^2) /
                     (1 - rad/(el * ib) * (gam/rad)^2*rad*mu*ib/area*om^2)
        amat[3, 4] = (rad / (el * ib)) /
                      (1 - (rad)/(el * ib) * (gam/rad)^2*rad*mu*ib/area*om^2)
        amat[4, 3] = -rad * mu * ib * om^2 / area
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om ^ 2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om ^ 2
        amat[6, 5] = -1

    elseif effct == :L_ALL  # All effects local.
        amat[1, 2] = 1
        amat[1, 5] = rad / (el * area)
        amat[2, 1] = -1
        amat[2, 3] = rad
        amat[2, 6] = rad * kn / (shear * area)
        amat[3, 4] = rad / (el * ib)
        amat[4, 3] = -rad * mu * ib * om^2 / area
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om^2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om^2
        amat[6, 5] = -1

    elseif effct == :N_NOSHEAR
        amat[1, 2] = (1 + rad / (el * area) * (gam/rad)^2 * rad * mu * om^2) /
                     (1 - rad / (el * area) * (gam/rad)^2 * rad * mu * om^2)
        amat[1, 5] = (rad / (el * area) * (1 + (gam/rad)^2)) /
                     (1 - rad / (el * area) * (gam/rad)^2 * rad * mu * om^2)
        amat[2, 1] = -1
        amat[2, 3] = rad
        amat[3, 2] = -(rad / (el * ib) * (gam/rad)^2 * rad^2*mu * om^2) /
                     (1 - rad/(el * ib) * (gam/rad)^2*rad*mu*ib/area*om^2)
        amat[3, 4] = (rad / (el * ib)) /
                      (1 - (rad)/(el * ib) * (gam/rad)^2*rad*mu*ib/area*om^2)
        amat[4, 3] = -rad * mu * ib * om^2 / area
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om ^ 2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om ^ 2
        amat[6, 5] = -1

    elseif effct == :N_NOAXIAL
        amat[1, 2] = 1
        amat[2, 1] = -(1 + rad * kn / (shear * area)*(gam/rad)^2*rad*mu*om^2) /
                     (1 - rad * kn / (shear * area) * (gam/rad)^2*rad*mu*om^2)
        amat[2, 3] = rad / (1 - rad * kn/(shear*area)*(gam/rad)^2*rad*mu*om^2)
        amat[2, 6] = (rad * kn / (shear * area) * (1+(gam/rad) ^ 2)) /
                     (1 - rad * kn / (shear * area) * (gam/rad)^2*rad*mu*om^2)
        amat[3, 2] = -(rad / (el * ib) * (gam/rad)^2 * rad^2*mu * om^2) /
                     (1 - rad/(el * ib) * (gam/rad)^2*rad*mu*ib/area*om^2)
        amat[3, 4] = (rad / (el * ib)) /
                      (1 - (rad)/(el * ib) * (gam/rad)^2*rad*mu*ib/area*om^2)
        amat[4, 3] = -rad * mu * ib * om^2 / area
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om ^ 2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om ^ 2

    elseif effct == :N_ROT
        amat[1, 2] = (1 + rad / (el * area) * (gam/rad)^2 * rad * mu * om^2) /
                     (1 - rad / (el * area) * (gam/rad)^2 * rad * mu * om^2)
        amat[1, 5] = (rad / (el * area) * (1 + (gam/rad)^2)) /
                     (1 - rad / (el * area) * (gam/rad)^2 * rad * mu * om^2)
        amat[2, 1] = -(1 + rad * kn / (shear * area)*(gam/rad)^2*rad*mu*om^2) /
                     (1 - rad * kn / (shear * area) * (gam/rad)^2*rad*mu*om^2)
        amat[2, 3] = rad / (1 - rad * kn/(shear*area)*(gam/rad)^2*rad*mu*om^2)
        amat[2, 6] = (rad * kn / (shear * area) * (1+(gam/rad) ^ 2)) /
                     (1 - rad * kn / (shear * area) * (gam/rad)^2*rad*mu*om^2)
        amat[3, 2] = -(rad / (el * ib) * (gam/rad)^2 * rad^2*mu * om^2)
        amat[3, 4] = rad / (el * ib)
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om ^ 2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om ^ 2
        amat[6, 5] = -1

    elseif effct == :N_NONE
        amat[1, 2] = 1
        amat[2, 1] = -1
        amat[2, 3] = rad
        amat[3, 2] = -(rad / (el * ib) * (gam/rad)^2 * rad^2*mu * om^2)
        amat[3, 4] = rad / (el * ib)
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om ^ 2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om ^ 2
        amat[6, 5] = -1

    elseif effct == :L_NONE
        amat[1, 2] = 1
        amat[2, 1] = -1
        amat[2, 3] = rad
        amat[3, 4] = rad / (el * ib)
        amat[4, 6] = -rad
        amat[5, 1] = -rad * mu * om ^ 2
        amat[5, 6] = 1
        amat[6, 2] = -rad * mu * om ^ 2
        amat[6, 5] = -1

    end
    return amat

end  # end of function element_fundamental(IPBeam).

""" element_fundamental(atype::OPBeam, element_no::Int64, om::Number)

Calculate fundamental matrix for each element of the beam (for out-of-plane analysis).

Parameters
----------
::IPBeam    bytpe       : Out-of-plane dynamic analysis
::Int64     element_no  : Element no.
::Float64   om          : Frequency

Returns
-------
A(6x6) fundamental matrix for a given element.
"""
function element_fundamental(atype::OPBeam, element_no::Int64, om::Number)

    effct = atype.effct
    bcurv = atype.bcurv
    thetat = atype.thetat
    rad0 = atype.rad0
    lam = atype.lam
    gam = atype.gamma
    eta = atype.eta
    kb = atype.kb
    ro = atype.ro
    el = atype.el
    nu = atype.nu

    h0 = sqrt(12) * rad0 * thetat / lam  # Height of cross-section
    b = h0 * 2/3  # Width of the cross_section
    shear = el / (2*(1+nu))  # Shear modulus

    radius(θ) = bcurv == :CONST ? rad0 : rad0 / (cos(θ))^3  # Gives radius for given angle.
    height(θ) = θ <= 0 ? h0 * (1 - eta * θ / thetat) : h0 * (1 + eta * θ / thetat);  # Gives height of cross-section
                                                                                         # for given angle.

    amat = zeros(6, 6);  # Fundamental matrix initialized.

    rad = (radius(coords(atype, element_no)) + radius(coords(atype, element_no + 1))) / 2  # Radius of this element
    h = (height(coords(atype, element_no)) + height(coords(atype, element_no + 1))) / 2  # Height of this element
    area = b * h  # Cross-sectional area of this element.
    mu = ro * area  # Weight of unit length beam for this particular element.
    inn = (h * b^3) / 12  # Moment of inertia (about binormal axis)
    ip = (h * b^3) / 12 + (b * h^3) / 12  # Polar moment of inertia.
    j = (h * b^3)/3 * (1 - 0.63 * b/h * (1 - b^4/(12 * h^4))) #  torsional constant. Important for
                                                              #  rectangular cross-sections.

    if effct == :N_ALL  # All effects nonlocal.
        amat[1, 2] = -rad/(1-(kb*rad)/(shear*area) * (gam/rad)^2 * rad*mu*om^2)
        amat[1, 6] = (kb * rad / (shear * area) * (1 + (gam/rad)^2)) /
                     (1-kb*rad/(shear*area) * (gam/rad)^2 * rad * mu * om^2)
        amat[2, 1] = (rad/(el * inn) * (gam/rad)^2 * rad^2 * mu * om^2) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[2, 3] = -(1 + rad/(el*inn) * (gam/rad)^2 * rad*mu/area*ip*om^2) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[2, 4] = (rad/(el*inn) * (1 + (gam/rad)^2)) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[3, 2] = (1+rad/(shear*j)*0.5*(gam/rad)^2 * rad*mu/area * inn * om^2) /
                     (1-rad/(shear*j)*(gam/rad)^2 * rad*mu/area * ip * om^2)
        amat[3, 5] = (rad/(shear*j) * (1 + 3/2*(gam/rad)^2)) /
                     (1-rad/(shear*j)*(gam/rad)^2 * rad*mu/area * ip * om^2)
        amat[4, 2] = -rad * mu / area * inn * om^2
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 3] = -rad * mu / area * ip * om^2
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

    elseif effct == :L_ALL  # All effects local.
        amat[1, 2] = -rad
        amat[1, 6] = kb * rad / (shear * area)
        amat[2, 3] = -1
        amat[2, 4] = rad/(el*inn)
        amat[3, 2] = 1
        amat[3, 5] = rad/(shear*j)
        amat[4, 2] = -rad * mu / area * inn * om^2
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 3] = -rad * mu/ area * ip * om^2
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

    elseif effct == :N_NOSHEAR  # No shear effect.
        amat[1, 2] = -rad
        amat[2, 1] = (rad/(el * inn) * (gam/rad)^2 * rad^2 * mu * om^2) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[2, 3] = -(1 + rad/(el*inn) * (gam/rad)^2 * rad*mu/area*ip*om^2) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[2, 4] = (rad/(el*inn) * (1 + (gam/rad)^2)) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[3, 2] = (1+rad/(shear*j)*0.5*(gam/rad)^2 * rad*mu/area * inn * om^2) /
                     (1-rad/(shear*j)*(gam/rad)^2 * rad*mu/area * ip * om^2)
        amat[3, 5] = (rad/(shear*j) * (1 + 3/2*(gam/rad)^2)) /
                     (1-rad/(shear*j)*(gam/rad)^2 * rad*mu/area * ip * om^2)
        amat[4, 2] = -rad * mu / area * inn * om^2
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 3] = -rad * mu / area * ip * om^2
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

        elseif effct == :N_NOROT_BEND  # No bending rotatory inertia.
        amat[1, 2] = -rad/(1-(kb*rad)/(shear*area) * (gam/rad)^2 * rad*mu*om^2)
        amat[1, 6] = (kb * rad / (shear * area) * (1 + (gam/rad)^2)) /
                     (1-kb*rad/(shear*area) * (gam/rad)^2 * rad * mu * om^2)
        amat[2, 1] = rad/(el * inn) * (gam/rad)^2 * rad^2 * mu * om^2
        amat[2, 3] = -(1 + rad/(el*inn) * (gam/rad)^2 * rad*mu/area*ip*om^2)
        amat[2, 4] = rad/(el*inn) * (1 + (gam/rad)^2)
        amat[3, 2] = 1 / (1-rad/(shear*j)*(gam/rad)^2 * rad*mu/area * ip * om^2)
        amat[3, 5] = (rad/(shear*j) * (1 + 3/2*(gam/rad)^2)) /
                     (1-rad/(shear*j)*(gam/rad)^2 * rad*mu/area * ip * om^2)
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 3] = -rad * mu / area * ip * om^2
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

        elseif effct == :N_NOROT_TORS  # No torsional rotatory inertia.
        amat[1, 2] = -rad/(1-(kb*rad)/(shear*area) * (gam/rad)^2 * rad*mu*om^2)
        amat[1, 6] = (kb * rad / (shear * area) * (1 + (gam/rad)^2)) /
                     (1-kb*rad/(shear*area) * (gam/rad)^2 * rad * mu * om^2)
        amat[2, 1] = (rad/(el * inn) * (gam/rad)^2 * rad^2 * mu * om^2) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[2, 3] = -1 / (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[2, 4] = (rad/(el*inn) * (1 + (gam/rad)^2)) /
                     (1-rad/(el*inn) * (gam/rad)^2 * rad * mu/area * inn * om^2)
        amat[3, 2] = (1+rad/(shear*j)*0.5*(gam/rad)^2 * rad*mu/area * inn * om^2)
        amat[3, 5] = (rad/(shear*j) * (1 + 3/2*(gam/rad)^2))
        amat[4, 2] = -rad * mu / area * inn * om^2
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

        elseif effct == :N_NONE  # No effects.
        amat[1, 2] = -rad
        amat[2, 1] = rad/(el * inn) * (gam/rad)^2 * rad^2 * mu * om^2
        amat[2, 3] = -1
        amat[2, 4] = rad/(el*inn) * (1 + (gam/rad)^2)
        amat[3, 2] = 1
        amat[3, 5] = rad/(shear*j) * (1 + 3/2*(gam/rad)^2)
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

    elseif effct == :L_NONE  # No effects.
        amat[1, 2] = -rad
        amat[2, 3] = -1
        amat[2, 4] = rad/(el*inn)
        amat[3, 2] = 1
        amat[3, 5] = rad/(shear*j) * (1 + 3/2*(gam/rad)^2)
        amat[4, 5] = -1
        amat[4, 6] = rad
        amat[5, 4] = 1
        amat[6, 1] = -rad * mu * om^2

    end
    return amat
end  # end of function fundamental(OPBeam).


function beam_fundamental(atype::Beam, om::Number)
    elnum = atype.elnum
    fundamatrix = []  # This vector consists of the fundamental matrices of each elements.
    for i = 1:elnum
        push!(fundamatrix, element_fundamental(atype, i, om))
    end
    return fundamatrix
end


function y1(atype::Beam, om::Number)
    no_of_elements = atype.elnum
    mat = []
    for i = 1:no_of_elements
        push!(mat, expm(element_fundamental(atype, i, om) * coords(atype, i)))
    end
    return mat
end  # end of function y1.


function y2(atype::Beam, om::Number)
    no_of_elements = atype.elnum
    mat = []
    for i = 1:no_of_elements
        push!(mat, expm(element_fundamental(atype, i, om) * coords(atype, i+1)))
    end
    return mat
end  # end of function y2.
