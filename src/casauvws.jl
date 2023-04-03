# Based on casatools.measures.measures.touvw.  See:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatools.measures.html?highlight=measure#casatools.measures.measures.touvw

module CASAUVW

import ERFA: gd2gc, WGS84, DJM0
import PyCASA: me, qa

"""
    radec2uvws(ra, dec, jd, obslla, xyz)

Convert `xyz` to `uvw` in UVW frame.  Uses CASA (via PyCASA).
"""
function radec2uvws(ra, dec, jd, obslla, xyz;
                    dut1=nothing, xp=nothing, yp=nothing)
    # Warn once if any EOP values are given
    if any((dut1, xp, yp) .!== nothing)
        @warn "caller-supplied EOP parameters ignored by $(@__MODULE__)" maxlog=1
    end

    # Setup the `measures` frame
    qq = qa.quantity
    obsitrf = gd2gc(WGS84, deg2rad(obslla.lon), deg2rad(obslla.lat), obslla.alt)
    me.doframe(me.position("itrf", qq.(obsitrf, "m")...))
    me.doframe(me.direction("J2000", qq(ra, "rad"), qq(dec, "rad")))
    me.doframe(me.epoch("UTC", qq(jd-DJM0, "d")))

    # Get xyz in measures form
    me_ants = me.position("itrf", qq(xyz[1,:], "m"),
                                  qq(xyz[2,:], "m"),
                                  qq(xyz[3,:], "m"))

    # Transform me_antxyz to antuvw Array
    bl = me.asbaseline(me_ants)
    reshape(me.touvw(bl)[2]["value"], 3, :)
end

"""
    radec2uvws(ra, dec, jd, obslla)

Return 3x3 rotation matrix that can be used to transform from XYZ frame to UVW
frame.  Uses CASA (via PyCASA).
"""
function radec2uvws(ra, dec, jd, obslla;
                    dut1=nothing, xp=nothing, yp=nothing)
    radec2uvws(ra, dec, jd, obslla, [1.0 0.0 0.0
                                     0.0 1.0 0.0
                                     0.0 0.0 1.0]; dut1, xp, yp)
end

"""
    radec2uvws(ra, dec, jd, obslla, antxyz, bls)

Return a 3xN UVW Matrix with one column per `bls` entry.  Uses CASA (via
PyCASA).
"""
function radec2uvws(ra, dec, jd, obslla, antxyz, bls;
                    dut1=nothing, xp=nothing, yp=nothing)
    antuvw = radec2uvws(ra, dec, jd, obslla, antxyz; dut1, xp, yp)
    # Compute baselines
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

function vla_location()
    me.addxvalue(me.observatory("VLA"))["value"]
end

end # module CASAUVW