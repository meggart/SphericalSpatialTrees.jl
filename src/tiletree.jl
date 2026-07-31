#Simple implementation of a node that can collect and group tiles and their spherical bounding boxes
#It can have children, which are new TileNodes, or leaves, which are usually RegularGridNodes.
struct TileNode{T}
    children::Vector{TileNode{T}}
    leaves::Vector{T}
    extent::SphericalCap
end
Base.show(io::IO,node::TileNode) = print(io,"Node with $(length(node.children)) children and $(length(node.leaves)) leaves")
function getchild(t::TileNode, i)
    i > length(t.children) ? t.leaves[i-length(t.children)] : t.children[i]
end
nchild(t::TileNode) = length(t.children) + length(t.leaves)
getchild(t::TileNode) = (getchild(t, i) for i in 1:nchild(t))
isleaf(::TileNode) = false
node_extent(t::TileNode) = t.extent