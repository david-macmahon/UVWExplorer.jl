# Based on https://gist.github.com/demorest/8c8bca4ac5860796593ca07006cc3df6
# by Paul Demerost

module PyUVW

import ERFA: gd2gc, WGS84
using PyCall

function __init__()
    py"""
    import numpy as np
    import astropy.coordinates as coord
    import astropy.time
    import astropy.units as u

    # Can do this to get updated IERS B values into astropy
    from astropy.utils import iers
    from astropy.utils import data
    iers_a = iers.IERS_A.open(data.download_file(iers.IERS_A_URL, cache=True))
    iers_b = iers.IERS_B.open(data.download_file(iers.IERS_B_URL, cache=True))
    iers_auto = iers.IERS_Auto.open()

    def vla_location():
        return coord.EarthLocation.of_site('vla')

    def ut1_utc(jd, aorb="B"):
        if aorb == "A":
            return iers_a.ut1_utc(jd)
        else:
            return iers_b.ut1_utc(jd)

    def getxp(jd, aorb="B"):
        if aorb == "A":
            return iers_a.pm_xy(jd)[0]
        else:
            return iers_b.pm_xy(jd)[0]

    def getyp(jd, aorb="B"):
        if aorb == "A":
            return iers_a.pm_xy(jd)[1]
        else:
            return iers_b.pm_xy(jd)[1]

    def compute_uvws(ra, dec, jd, obsxyz, antxyz):
        # Time of observation:
        t = astropy.time.Time(jd-2400000.5, format='mjd', scale='utc')

        # Format antenna positions and VLA center as EarthLocation.
        antpos_ap = coord.EarthLocation(antxyz[:,0]*u.m, antxyz[:,1]*u.m, antxyz[:,2]*u.m)
        vla = coord.EarthLocation(obsxyz[0]*u.m, obsxyz[1]*u.m, obsxyz[2]*u.m)

        # Convert antenna pos terrestrial to celestial.  For astropy use 
        # get_gcrs_posvel(t)[0] rather than get_gcrs(t) because if a velocity 
        # is attached to the coordinate astropy will not allow us to do additional 
        # transformations with it (https://github.com/astropy/astropy/issues/6280)
        vla_p, vla_v = vla.get_gcrs_posvel(t)
        antpos_c_ap = coord.GCRS(antpos_ap.get_gcrs_posvel(t)[0], 
                obstime=t, obsgeoloc=vla_p, obsgeovel=vla_v)

        # Define the UVW frame relative to a certain point on the sky.  There are
        # two versions, depending on whether the sky offset is done in ICRS 
        # or GCRS:
        pnt = coord.SkyCoord(ra*u.rad, dec*u.rad, frame='icrs')
        #frame_uvw = pnt.skyoffset_frame() # ICRS
        frame_uvw = pnt.transform_to(antpos_c_ap).skyoffset_frame() # GCRS

        # Rotate antenna positions into UVW frame.
        antpos_uvw_ap = antpos_c_ap.transform_to(frame_uvw)

        # SkyOffsetFrame coords seem to come out as WUV, so shuffle into UVW order:
        return coord.CartesianRepresentation(
                x = antpos_uvw_ap.cartesian.y,
                y = antpos_uvw_ap.cartesian.z,
                z = antpos_uvw_ap.cartesian.x)
    """
end # function __init__

function radec2uvws(ra, dec, jd, obslla, xyz;
                    dut1=nothing, xp=nothing, yp=nothing)
    # Warn once if any EOP values are given
    if any((dut1, xp, yp) .!== nothing)
        @warn "caller-supplied EOP parameters ignored by $(@__MODULE__)" maxlog=1
    end

    obsitrf = gd2gc(WGS84, deg2rad(obslla.lon), deg2rad(obslla.lat), obslla.alt)
    py"compute_uvws"(ra, dec, jd, obsitrf, xyz').xyz
end

function radec2uvws(ra, dec, jd, obslla;
                    dut1=nothing, xp=nothing, yp=nothing)
    radec2uvws(ra, dec, jd, obslla, [1.0 0.0 0.0
                                     0.0 1.0 0.0
                                     0.0 0.0 1.0]; dut1, xp, yp)
end

function radec2uvws(ra, dec, jd, obslla, antxyz, bls;
                    dut1=nothing, xp=nothing, yp=nothing)
    antuvw = radec2uvws(ra, dec, jd, obslla, antxyz; dut1, xp, yp)
    # Compute baselines
    mapreduce(hcat, bls) do (a1, a2)
        antuvw[:, a2] - antuvw[:, a1]
    end
end

function vla_location()
    [py"vla_location"().value.tolist()...]
end

function getΔUT1(jd, aorb="B")
    py"ut1_utc"(jd, aorb)[]
end

function getxp(jd, aorb="B")
    py"getxp"(jd, aorb)[]
end

function getyp(jd, aorb="B")
    py"getyp"(jd, aorb)[]
end

end # module PyUVW