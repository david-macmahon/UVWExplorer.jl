module TopoUVW

import Rotations: RotY, RotZ
import ..ERFAHelpers: radec2obs

function hadec2uvws(hob, dob, obslla)
    # Rotation matrix to transform XYZ frame to UVW frame in direction of
    # (hob, dob) at obslla.lon.
    [0 1 0
     0 0 1
     1 0 0] * RotY(dob) * RotZ(hob-deg2rad(obslla.lon))
end

function hadec2uvws(hob, dob, obslla, xyz)
    # Transform antxyz to antuvw
    hadec2uvws(hob, dob, obslla) * xyz
end

function hadec2uvws(hob, dob, obslla, antxyz, bls)
    antuvw = hadec2uvws(hob, dob, obslla, antxyz)
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

function azel2uvws(aob, eob, obslla)
    hob, dob = ae2hd(aob, eob, deg2rad(obslla.lat))
    hadec2uvws(hob, dob, obslla)
end

function azel2uvws(aob, eob, obslla, xyz)
    hob, dob = ae2hd(aob, eob, deg2rad(obslla.lat))
    hadec2uvws(hob, dob, obslla, xyz)
end

function azel2uvws(aob, eob, obslla, antxyz, bls)
    hob, dob = ae2hd(aob, eob, deg2rad(obslla.lat))
    hadec2uvws(hob, dob, obslla, antxyz, bls)
end

function azza2uvws(aob, zob, obslla)
    azel2uvws(aob, π/2-zob, obslla)
end

function azza2uvws(aob, zob, obslla, antxyz)
    azel2uvws(aob, π/2-zob, obslla, antxyz)
end

function azza2uvws(aob, zob, obslla, antxyz, bls)
    azel2uvws(aob, π/2-zob, obslla, antxyz, bls)
end

function radec2uvws(ra, dec, jd, obslla;
                    dut1=nothing, xp=nothing, yp=nothing)
    # Convert nothings to default values
    dut1 = something(dut1, 0.0)
    xp = something(xp, 0.0)
    yp = something(yp, 0.0)

    obs = radec2obs(ra, dec, jd, obslla; dut1, xp, yp)
    hadec2uvws(obs.hob, obs.dob, obslla)
end

function radec2uvws(ra, dec, jd, obslla, xyz;
                    dut1=nothing, xp=nothing, yp=nothing)
    radec2uvws(ra, dec, jd, obslla; dut1, xp, yp) * xyz
end

function radec2uvws(ra, dec, jd, obslla, antxyz, bls;
                    dut1=nothing, xp=nothing, yp=nothing)
    antuvw = radec2uvws(ra, dec, jd, obslla, antxyz; dut1, xp, yp)
    # Compute baselines
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

end # module TopoUVW