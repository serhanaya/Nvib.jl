""" detrmt(atype::Beam, om::Number)

Calculate determinant of A of size [6 * number of elements, 6 * number of elements]
"""
function detrmt(atype::Beam, om::Number)

    elnum = atype.elnum
    bcond = atype.bcond
    rad0 = atype.rad0
    lam = atype.lam
    phit = atype.phit

    rowdim = 6 * elnum
    coldim = 6 * elnum

    ytop = zeros(3, coldim)
    ybottom = zeros(3, coldim)
    ymid = zeros(rowdim - 6, coldim)

    mat = beam_fundamental(atype, om)
    yy1 = y1(atype, om)
    yy2 = y2(atype, om)

    if bcond == :CC  # clamped-clamped

        ytop[:, 1:6] = yy1[1][1:3, :]
        ybottom[:, (coldim-5):end] = yy2[elnum][1:3, :]

    elseif bcond == :HH  # hinged-hinged

        ytop[1:2, 1:6] = yy1[1][1:2, :]
        ytop[3, 1:6] = yy1[1][4, :]
        ybottom[1:2, (coldim-5):end] = yy2[elnum][1:2, :]
        ybottom[3, (coldim-5):end] = yy2[elnum][4, :]

    elseif bcond == :HC  # hinged-clamped

        ytop[1:2, 1:6] = yy1[1][1:2, :]
        ytop[3, 1:6] = yy1[1][4, :]
        ybottom[:, (coldim - 5):end] = yy2[elnum][1:3, :]

    elseif bcond == :CF  # clamped-free

        ytop[:, 1:6] = yy1[1][1:3, :]
        ybottom[:, (coldim - 5):end] = yy2[elnum][4:end, :]

    elseif bcond == :FC  # free-clamped

        ytop[:, 1:6] = yy1[1][4:end, :]
        ybottom[:, (coldim - 5):end] = yy2[elnum][1:3, :]

    elseif bcond == :FF  # free-free

        ytop[:, 1:6] = yy1[1][4:end, :]
        ybottom[:, (coldim - 5):end] = yy2[elnum][4:end, :]

    end

    for i in 1:(elnum-1)

        ymid[ 6(i-1)+1 : 6i, 6(i-1)+1 : 6i ] = expm(mat[i] * coords(atype, i))
        ymid[ 6(i-1)+1 : 6i, 6i + 1 : 6(i+1)  ] = - expm(mat[i+1] * coords(atype, i))

    end

    ymat = vcat(ytop, ymid, ybottom)

    return det(ymat)

end  # end of function detrmt.
