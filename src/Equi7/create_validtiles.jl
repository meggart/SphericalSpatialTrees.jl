import Shapefile
const ZONES = ("AF","AN","AS","EU","NA","OC","SA")
function get_tilecoord(s) 
    r = match(r"E(\d\d\d)N(\d\d\d)T1",s)
    parse(Int,r.captures[1]),parse(Int,r.captures[2])
end

for zone in ZONES
    tab2 = Shapefile.Table("../Equi7Grid-0.2.6/src/equi7grid/grids/$zone/GEOG/EQUI7_V14_$(zone)_GEOG_TILE_T1")
    df2 = DataFrame(tab2)
    allcoords = get_tilecoord.(df2.TILE)
    open("./dev/SphericalSpatialTrees/src/Equi7/tiles/$zone",create=true, write=true) do io
        print(io,allcoords)
    end
end