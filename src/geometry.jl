module geometry
# import base operators that will be overloaded. This is just
# to simplify arithmetic on Points (i.e. scaling and 
# translation)
import Base: *, +, -

# We 

# Export types and methods
export Point, *, +, -, Triangulation, Vertex, Line, Line3, Triangle, 
		Triangle6, Quad, Quad9, Mesh, LoadMesh


"""
	Point{dim}
	Point(x::Float64)
	Point(x::Float64, y::Float64)
	Point(x::Float64, y::Float64, z::Float64)
	Point(x::NTuple{N, Float64}) where N
A parametrized type representing a point in `dim` 
dimensional space.
"""
struct Point{dim}
	x::NTuple{dim, Float64}
	function Point(x::Float64)
		new{1}((x,))
	end
	function Point(x::Float64, y::Float64)
		new{2}((x, y))
	end
	function Point(x::Float64, y::Float64, z::Float64)
		new{3}((x, y, z))
	end
	function Point(x::NTuple{N, Float64}) where N
		@assert (N > 0 && N <= 3) "Points must be in Euclidean space"
		new{N}(x)
	end
end


"""
	*(a::Float64, p::Point{dim}) where dim
	*(p::Point{dim}, a::Float64) where dim	
Multiply a point `p` by a scalar `a` and return a new 
`Point` object. Commutative.
"""
function *(a::Float64, p::Point{dim}) where dim
	return Point(a.*p.x)
end


function *(p::Point{dim}, a::Float64) where dim
	return Point(a.*p.x)
end

"""
	+(p1::Point{dim}, p2::Point{dim}) where dim
Add the components of two points `p1` and `p2`. 
Return a new `Point` object.
"""
function +(p1::Point{dim}, p2::Point{dim}) where dim
	return Point(p1.x .+ p2.x)
end

"""
	-(p1::Point{dim}, p2::Point{dim}) where dim
Subtract the coordinates of point `p2` from point `p1`
component-wise. Return a new `Point` object.
"""
function -(p1::Point{dim}, p2::Point{dim}) where dim
	return Point(p1.x .- p2.x)
end


"""
	Triangulation{dim, spacedim}
Abstract supertype for all `dim` dimensional geometric 
triangulations in `spacedim` dimensional space.
# Example
- `Line` is a 1D geometric entity that can be embedded in 
1D, 2D, or 3D spacedim dimensions
"""
abstract type Triangulation{dim, spacedim} end

"""
	Vertex <: Triangulation{0, spacedim}
0D point element
# Attributes
- `node::NTuple{1, Tuple{Int64, Point{spacedim}}}`
1-Tuple of `(global_node_number, coordinates)`
"""
struct Vertex{spacedim} <: Triangulation{0, spacedim}
	node::NTuple{1, Tuple{Int64, Point{spacedim}}}
	function Vertex(node_ids::NTuple{1, Int64},
		coordinates::NTuple{1, Point{spacedim}}) where spacedim
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end

"""
	Line <: Triangulation{1, spacedim}
1D linear (2 node) line element
# Attributes
- `node::NTuple{2, Tuple{Int64, Point{spacedim}}}`
2-Tuple of `(global_node_number, coordinates)`
"""
struct Line{spacedim} <: Triangulation{1, spacedim}
	node::NTuple{2, Tuple{Int64, Point{spacedim}}}
	function Line(node_ids::NTuple{2, Int64},
		coordinates::NTuple{2, Point{spacedim}}) where spacedim
		@assert spacedim >= 1
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end

"""
	Line3 <: Triangulation{1, spacedim}
1D quadratic (3 node) line element
# Attributes
- `node::NTuple{3, Tuple{Int64, Point{spacedim}}}`
3-tuple of `(global_node_number, coordinates)`
"""
struct Line3{spacedim} <: Triangulation{1, spacedim}
	node::NTuple{3, Tuple{Int64, Point{spacedim}}}
	function Line3(node_ids::NTuple{3, Int64},
		coordinates::NTuple{3, Point{spacedim}}) where spacedim
		@assert spacedim >= 1
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end


"""
	Triangle <: Triangulation{2, spacedim}
2D linear (3 node) triangular element
# Attributes
- `node::NTuple{3, Tuple{Int64, Point{spacedim}}}`
3-tuple of `(global_node_number, coordinates)`
"""
struct Triangle{spacedim} <: Triangulation{2, spacedim}
	node::NTuple{3, Tuple{Int64, Point{spacedim}}}
	function Triangle(node_ids::NTuple{3, Int64},
		coordinates::NTuple{3, Point{spacedim}}) where spacedim
		@assert spacedim >= 2
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end

"""
	Triangle6{spacedim} <: Triangulation{2, spacedim}
2D quadratic (6 node) triangular element
# Attributes
- `node::NTuple{6, Tuple{Int64, Point{spacedim}}}`
6-tuple of `(global_node_number, coordinates)`
"""
struct Triangle6{spacedim} <: Triangulation{2, spacedim}
	node::NTuple{6, Tuple{Int64, Point{spacedim}}}
	function Triangle6(node_ids::NTuple{6, Int64},
		coordinates::NTuple{6, Point{spacedim}}) where spacedim
		@assert spacedim >= 2
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end

"""
	Quad{spacedim} <: Triangulation{2, spacedim}
2D linear (4 node) quadrilateral element
# Attributes
- `node::NTuple{4, Tuple{Int64, Point{spacedim}}}`
4-tuple of `(global_node_number, coordinates)`
"""
struct Quad{spacedim} <: Triangulation{2, spacedim}
	node::NTuple{4, Tuple{Int64, Point{spacedim}}}
	function Quad(node_ids::NTuple{4, Int64},
		coordinates::NTuple{4, Point{spacedim}}) where spacedim
		@assert spacedim >= 2
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end


"""
	Quad9{spacedim} <: Triangulation{2, spacedim}
2D quadratic (9 node) quadrilateral element
# Attributes
- `node::NTuple{9, Tuple{Int64, Point{spacedim}}}`
9-tuple of `(global_node_number, coordinates)`
"""
struct Quad9{spacedim} <: Triangulation{2, spacedim}
	node::NTuple{9, Tuple{Int64, Point{spacedim}}}
	function Quad9(node_ids::NTuple{9, Int64},
		coordinates::NTuple{9, Point{spacedim}}) where spacedim
		@assert spacedim >= 2
		node_data = zip(node_ids, coordinates)
		new{spacedim}(tuple(node_data...))
	end
end


"""
	Mesh{spacedim}
A mesh object in `spacedim` dimensional space.
# Attributes
- `elements::Dict{String, Array{Triangulation{dim, spacedim} where dim, 1}}`
Dictionary of triangulation objects accessed by keys that
indicate some logical tag for the elements, for example "Body",
"Surface", "Boundary", etc.
- `nodes::Dict{String, Array{Point{spacedim}, 1}}
Dictionary of point objects representing the nodes of the mesh, accessed by keys 
that indicate some logical tag for the nodes, for example "Body", "Boundary", etc.
"""
struct Mesh{spacedim}
	elements::Dict{String, Array{Triangulation{dim, spacedim} where dim, 1}}
	nodes::Dict{String, Array{Point{spacedim}, 1}}
	function Mesh(spacedim::Int64)
		elements = Dict{String, Array{Triangulation{dim, spacedim} where dim, 1}}()
		nodes = Dict{String, Array{Point{spacedim}, 1}}()
		new{spacedim}(elements, nodes)
	end
end

"""
	AssembleQuad(mesh::Mesh, mesh_data, tagToGroup::Dict{Int64, String})
Read `mesh_data[:cells]["quad"]` and construct `Quad` elements. Push these 
elements into `mesh.elements[group_name]` where `group_name` is obtained
from `tagToGroup`.
"""
function AssembleQuad(mesh::Mesh{spacedim}, mesh_data, tagToGroup::Dict{Int64, String}) where spacedim
	for i in 1:size(mesh_data[:cells]["quad"])[1]
		tag = mesh_data[:cell_data]["quad"]["gmsh:physical"][i]
		group_name = tagToGroup[tag]
		node_ids = tuple(mesh_data[:cells]["quad"][i,:]...)
		# gmsh uses 0-based indexing, so shift node numbers by 1
		node_ids = node_ids .+ 1
		coordinates = tuple([Point(tuple(mesh_data[:points][i,1:spacedim]...)) for i in node_ids]...)
		element = Quad(node_ids, coordinates)
		push!(mesh.elements[group_name], element)
	end
end


"""
	LoadMesh(mesh_data, spacedim::Int64)
Use the `mesh_data` `PyObject` to construct a `Mesh` object in
`spacedim` dimensions. `mesh_data` may be constructed using `PyCall`
and the `meshio` package. 
If `spacedim < 3`, drop the trailing coordinate values of nodal 
coordinates.
"""
function LoadMesh(mesh_data, spacedim::Int64)
	# Initialize the mesh object 
	mesh = Mesh(spacedim)

	# go through all the keys of field_data. These are the names of the physical
	# groups. Construct a dictionary tagToGroup such that
	# tagToGroup[tag::Int] = physical_group_name
	tagToGroup = Dict{Int64, String}()
	for key in keys(mesh_data[:field_data])
		tag = mesh_data[:field_data][key][1]
		tagToGroup[tag] = key
		# Initialize an empty array in the mesh object corresponding
		# to the current key
		mesh.elements[key] = []
	end

	for key in keys(mesh_data[:cells])
		if key == "quad"
			AssembleQuad(mesh, mesh_data, tagToGroup)
		end
	end
end

















# module geometry ends here
end
# module geometry ends here