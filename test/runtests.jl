using UVWExplorer
using Test

ra  = 3.2686162621088966
dec = 0.03582094364245917
jd  = 2.4598932634780095e6

obslla = (
    lat =   34.07861111111111,
    lon = -107.61777777777777,
    alt = 2124.0
)

antxyz = [
     75.3561   -1407.4597
    489.4072     -77.4865
    721.5445    -735.1909
]

topouvw_expected = [-617.4059528467405, -1505.1348066165679, 1412.6458483084875]
gcrsuvw_expected = [-617.7848812531224, -1504.9793145403532, 1412.645848308487]

topo = TopoUVW.radec2uvws(ra, dec, jd, obslla)
gcrs = GCRSUVW.radec2uvws(ra, dec, jd, obslla)
py = PyUVW.radec2uvws(ra, dec, jd, obslla)
casa = CASAUVW.radec2uvws(ra, dec, jd, obslla)

topoantuvw = topo * antxyz
gcrsantuvw = gcrs * antxyz
pyantuvw   = py   * antxyz
casaantuvw = casa * antxyz

topouvw = topoantuvw[:,2] - topoantuvw[:,1]
gcrsuvw = gcrsantuvw[:,2] - gcrsantuvw[:,1]
pyuvw   = pyantuvw[  :,2] - pyantuvw[  :,1]
casauvw = casaantuvw[:,2] - casaantuvw[:,1]

@testset "Topo UVW" begin
    @test topouvw ≈ topouvw_expected
    @test topouvw ≈ gcrsuvw atol=0.5
    @test topouvw ≈ pyuvw atol=0.5
    @test topouvw ≈ casauvw atol=0.5
end

@testset "GCRS UVW" begin
    @test gcrsuvw ≈ gcrsuvw_expected
    @test gcrsuvw ≈ pyuvw atol=1.2e-3
    @test gcrsuvw ≈ casauvw atol=3.5e-3
    @test pyuvw ≈ casauvw atol=3e-3
end