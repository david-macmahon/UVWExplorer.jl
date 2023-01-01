# Based on casatools.measures.measures.touvw.  See:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatools.measures.html?highlight=measure#casatools.measures.measures.touvw

module CASAUVW

import ERFA: gd2gc, WGS84, DJM0
import PyCASA: me, qa

"""
    radec2uvws(ra, dec, jd, obslla, antxyz, bls)

Return a 3xN UVW Matrix with one column per `bls` entry.  Uses CASA (via
PyCASA).
"""
function radec2uvws(ra, dec, jd, obslla, antxyz, bls;
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

    # Get antxyz in measures form
    me_ants = me.position("itrf", qq(antxyz[1,:], "m"),
                                  qq(antxyz[2,:], "m"),
                                  qq(antxyz[3,:], "m"))

    # Transform me_antxyz to antuvw Array
    bl = me.asbaseline(me_ants)
    antuvw = reshape(me.touvw(bl)[2]["value"], 3, :)

    # Compute UVWs for each baseline
    mapreduce(hcat, bls) do (a1i, a2i)
        antuvw[:, a2i] - antuvw[:, a1i]
    end
end

function vla_location()
    me.addxvalue(me.observatory("VLA"))["value"]
end

end # module CASAUVW