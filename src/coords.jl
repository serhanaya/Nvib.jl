""" coords(atype::Beam, nodenum)

Determine the coordinate of the selected node in radians.

Parameters
----------
::Beam      atype   : Analysis type (In-plane or Out-of-plane analysis)
::Int64     nodenum : Node number (e.g. coords(atype, 0) gives the first node's coordinate)

Returns
-------
::Float64   theta   : Coordinate in radians.
"""
function coords(atype::Beam, nodenum)

    elnum = atype.elnum
    thetat = atype.thetat

    thetaa = -thetat/2
    thetab = thetat/2
    thetael = (thetab - thetaa) / elnum
    thetax = zeros(elnum + 1)  # vector initialization for coordinate of element x
                             # size = [total number of nodes]

    thetax[1] = thetaa  # Coordinate of nodenum=1 is thetaa rad.
    for i in 2:length(thetax)
        thetax[i] = thetaa + (i-1) * thetael
    end
    return thetax[nodenum]
end  # end of function coords.
