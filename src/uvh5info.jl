import HDF5

"""
    uvh5info(uvh5_name::AbstractString
    uvh5info(uvh5_file::HDF5.File)
    uvh5info(uvh5_header::HDF5.Group)

Extract info from UVH5 file or `Header` for first time sample.  Extracted items
are returned in a `NamedTuple` with the following fields:

* `ra` - right ascension of phase center in radians
* `dec` - declination of phase center in radians
* `srcpos` - `(; ra, dec)` of phase center in radians
* `jd` - Julian date of `uvws`
* `obslla` - observatory location as an LLA object
* `antxyz` - antenna positions as XYZ with origin at observatory location
* `bls` - baseline tuples `(a1idx, a2idx)` where a1idx/a2idx are indices into
  `antxyz` (NOT antenna numbers).
* `uvws` - UVW array from the UVH5 file

This function can also be passed an HDF5.File object or name of a UVH5 file as
well.
"""
function uvh5info(uvh5_header::HDF5.Group)
    h5h = uvh5_header
    a1s = h5h["ant_1_array"][]
    a2s = h5h["ant_2_array"][]
    jds = h5h["time_array"][]
    sel = [a1s[i]!=a2s[i] && jds[i]==jds[1] for i in eachindex(jds)]
    antnums = h5h["antenna_numbers"][]
    # Map antenna numbers to antenna indices
    antmap = Dict(n=>i for (i,n) in enumerate(antnums))
    a1is = map(a->antmap[a], a1s[sel])
    a2is = map(a->antmap[a], a2s[sel])
    (
        ra=h5h["phase_center_ra"][],
        dec=h5h["phase_center_dec"][],
        jd = jds[1],
        obslla = (lat=h5h["latitude"][], lon=h5h["longitude"][], alt=h5h["altitude"][]),
        antxyz = h5h["antenna_positions"][],
        bls = collect(zip(a1is, a2is)),
        uvws = h5h["uvw_array"][][:,sel]
    )
end

function uvh5info(uvh5_file::HDF5.File)
    uvh5info(uvh5_file["Header"])
end

function uvh5info(uvh5_name::AbstractString)
    uvh5info(HDF5.h5open(uvh5_name))
end
