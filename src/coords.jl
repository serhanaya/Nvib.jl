""" coords(atype::Beam, nodenum)

Determine the coordinate of the selected node in radians.

Parameters
----------
::Beam      atype   : Analysis type (In-plane or Out-of-plane analysis)
::Int64     nodenum : Node number (e.g. coords(atype, 0) gives the first node's coordinate)

Returns
-------
::Float64   phi     : Coordinate in radians.
"""
function coords(atype::Beam, nodenum)

    elnum = atype.elnum
    phit = atype.phit

    phia = -phit/2
    phib = phit/2
    phiel = (phib - phia) / elnum
    phix = zeros(elnum + 1)  # vector initialization for coordinate of element x
                             # size = [total number of nodes]

    phix[1] = phia  # Coordinate of nodenum=1 is phia rad.
    for i in 2:length(phix)
        phix[i] = phia + (i-1) * phiel
    end
    return phix[nodenum]
end  # end of function coords.
