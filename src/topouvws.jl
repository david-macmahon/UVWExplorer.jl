module TopoUVW

import Rotations: RotY, RotZ
import ..ERFAHelpers: radec2obs

function hadec2uvws(hob, dob, obslla, antxyz, bls)
    # Rotation matrix to transform XYZ frame to UVW frame in direction of
    # (hob, dob) at obslla.lon.
    xyz2uvw = [0 1 0
               0 0 1
               1 0 0] * RotY(dob) * RotZ(hob-deg2rad(obslla.lon))

    # Transform antxyz to antuvw
    antuvw = xyz2uvw * antxyz

    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

function azel2uvws(aob, eob, obslla, antxyz, bls)
    hob, dob = ae2hd(aob, eob, deg2rad(obslla.lat))
    hadec2uvws(hob, dob, obslla, antxyz, bls)
end

function azza2uvws(aob, zob, obslla, antxyz, bls)
    axel2uvws(aob, π/2-zob, obslla, antxyz, bls)
end

function radec2uvws(ra, dec, jd, obslla, antxyz, bls;
                    dut1=nothing, xp=nothing, yp=nothing)
    # Convert nothings to default values
    dut1 = something(dut1, 0.0)
    xp = something(xp, 0.0)
    yp = something(yp, 0.0)

    obs = radec2obs(ra, dec, jd, obslla; dut1, xp, yp)
    hadec2uvws(obs.hob, obs.dob, obslla, antxyz, bls)
end

end # module TopoUVW