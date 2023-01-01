module GCRSUVW

import Rotations: RotY, RotZ
import ERFA: DAS2R, c2t06a, apco13, ldsun, ab, s2c, c2s, utctai, taitt

function radec2uvws(ra, dec, jd, obslla, antxyz, bls;
                    dut1=nothing, xp=nothing, yp=nothing)
    # Convert nothings to default values
    dut1 = something(dut1, 0.0)
    xp = something(xp, 0.0)
    yp = something(yp, 0.0)

    # Get TT
    tai1, tai2 = utctai(jd, 0)
    tt1, tt2 = taitt(tai1, tai2)
    
    # Convert lat/lon to radians
    lat = deg2rad(obslla.lat)
    lon = deg2rad(obslla.lon)
    alt = obslla.alt

    # Convert ra,dec from ICRS to GCRS.  For completeness, this involves space
    # motion, parallax, light deflection, and abberation, but here we just do
    # light deflection due to the sun and abberation.  We use `apco13` rather
    # than `apcg13` because the former includes observatory position/velocity in
    # the astrometry parameters used to calculate light deflection and
    # abberation whereas the latter does not.
    astrom = apco13(jd, 0, dut1, lon, lat, alt, xp, yp, 0, 0, 0, 0)[1]
    srcvec = s2c(ra, dec)
    ldvec = ldsun(srcvec, astrom.eh, astrom.em)
    gcvec = ab(ldvec, astrom.v, astrom.em, astrom.bm1)
    ra_gcrs, dec_gcrs = c2s(gcvec)

    # Rotation matrix to transform GCRS frame to UVW frame in direction of
    # (ra_gcrs, dec_gcrs).
    gcrs2uvw = [0 1 0
                0 0 1
                1 0 0] * RotY(dec_gcrs) * RotZ(-ra_gcrs)

    # Rotation matrix to transform ITRF (or XYZ) to GCRS (oriented) is the
    # transpose of ERFA's celestial-to-terrestrial ("c2t") matrix.
    xyz2gcrs = c2t06a(tt1, tt2, jd, dut1/86400, xp*DAS2R, yp*DAS2R)'

    # Transform antxyz to antuvw
    antuvw = gcrs2uvw * xyz2gcrs * antxyz

    # Compute baselines
    mapreduce(hcat, bls) do (a1i, a2i)
        antuvw[:, a2i] - antuvw[:, a1i]
    end
end

end # module GCRSUVW