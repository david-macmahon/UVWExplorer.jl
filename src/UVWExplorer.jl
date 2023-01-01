module UVWExplorer

export radec2uvws

include("erfa_helpers.jl")
import .ERFAHelpers: radec2obs

include("uvh5info.jl")
export uvh5info

include("topouvws.jl")
export TopoUVW

include("gcrsuvws.jl")
export GCRSUVW

include("pyuvws.jl")
export PyUVW

include("casauvws.jl")
export CASAUVW

function radec2uvws(mod, ra, dec, jd, obslla, antxyz, bls; dut1=0, xp=0, yp=0)
    @assert mod in (TopoUVW, GCRSUVW, PyUVW, CASAUVW) "unsupported module: $mod"
    mod.radec2uvws(ra, dec, jd, obslla, antxyz, bls; dut1, xp, yp)
end

end # module UVWExplorer