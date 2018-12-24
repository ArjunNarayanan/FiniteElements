module geometry
# import base operators that will be overloaded. This is just
# to simplify arithmetic on Points (i.e. scaling and 
# translation)
import Base: *, +, -, ==, getindex

using StaticArrays

# Export types and methods
export Point, *, +, -, ==, getindex, Triangulation, Vertex, Line, Triangle, 
		Triangle, Quadrilateral, Mesh, loadMesh


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
	x::SVector{dim, Float64}
	function Point(x)
		new{1}((x,))
	end
	function Point(x, y)
		new{2}((x, y))
	end
	function Point(x, y, z)
		new{3}((x, y, z))
	end
	function Point(x::NTuple{N, T}) where {N,T}
		new{N}(x)
	end
	function Point(x::SVector{dim, Float64}) where dim
		new{dim}(x)
	end
end

"""
	getindex(p::Point, i::Int64)
Return the `i`th component of the coordinates of point `p`
i.e. `p.x[i]`.
"""
function getindex(p::Point, i::Int64)
	return p.x[i]
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
	==(p1::Point{dim}, p2::Point{dim}) where dim
Check if two points are equivalent by comparing their
corresponding tuples.
i.e. `p1 == p2` if and only if `p1.x == p2.x`
"""
function ==(p1::Point{dim}, p2::Point{dim}) where dim
	return (p1.x == p2.x)
end


"""
	Triangulation{N, dim, spacedim}
Abstract supertype for all `dim` dimensional geometric 
triangulations in `spacedim` dimensional space with `N` nodes.
# Examples
- 2 node `Line{2,spacedim}` is a linear line element that can be
embedded in `spacedim` dimensional space.
- 3 node `Line{3,spacedim}` is a quadratic line element that can be
embedded in `spacedim` dimensional space.
"""
abstract type Triangulation{P, dim, spacedim} end

"""
	Vertex{spacedim} <: Triangulation{1, 0, spacedim}
1-node, 0 dimensional point element.
# Attributes
- `points::Tuple{Point{spacedim}}
Coordinates of the point specifying the vertex.
- `nodes::Tuple{Int64}
Global node ID of the corresponding node
"""
struct Vertex{spacedim} <: Triangulation{1, 0, spacedim}
	nodes::Tuple{Int64}
	points::Tuple{Point{spacedim}}
	function Vertex(nodes::Tuple{Int64},
		points::Tuple{Point{spacedim}}) where spacedim
			new{spacedim}(nodes, points)
	end
end

"""
	Line{N, spacedim} <: Triangulation{N, 1, spacedim}
1D line element with `N` nodes in `spacedim` dimensional
space.
# Attributes
- `nodes::NTuple{N, Int64}` - the global node numbers
- `points::NTuple{N, Point{spacedim}}` - the corresponding coordinates
"""
struct Line{N, spacedim} <: Triangulation{N, 1, spacedim}
	nodes::NTuple{N, Int64}
	points::NTuple{N, Point{spacedim}}
	function Line{N}(nodes::NTuple{N, Int64},
		points::NTuple{N, Point{spacedim}}) where {N,spacedim}
		new{N,spacedim}(nodes, points)
	end
end


"""
	Triangle{N, spacedim} <: Triangulation{N, 2, spacedim}
2D triangular element with `N` nodes in `spacedim` dimensional
space.
# Attributes
- `nodes::NTuple{N, Int64}` - the global node numbers
- `points::NTuple{N, Point{spacedim}}` - the corresponding coordinates
"""
struct Triangle{N, spacedim} <: Triangulation{N, 2, spacedim}
	nodes::NTuple{N, Int64}
	points::NTuple{N, Point{spacedim}}
	function Triangle{N}(nodes::NTuple{N, Int64},
		points::NTuple{N, Point{spacedim}}) where {N,spacedim}
		new{N, spacedim}(nodes, points)
	end
end



"""
	Quadrilateral{N, spacedim} <: Triangulation{N, 2, spacedim}
2D quadrilateral element with `N` nodes in `spacedim` dimensional
space.
# Attributes
- `nodes::NTuple{N, Int64}` - the global node numbers
- `points::NTuple{N, Point{spacedim}}` - the corresponding coordinates
"""
struct Quadrilateral{N, spacedim} <: Triangulation{N, 2, spacedim}
	nodes::NTuple{N, Int64}
	points::NTuple{N, Point{spacedim}}
	function Quadrilateral{N}(nodes::NTuple{N, Int64},
		points::NTuple{N, Point{spacedim}}) where {N,spacedim}
		new{N, spacedim}(nodes, points)
	end
end


"""
	Mesh{spacedim}
A mesh object in `spacedim` dimensional space.
# Attributes
- `elements::Array{Triangulation{P, dim, spacedim} where {P, dim}, 1}`
Array of triangulation objects accessed by keys that
indicate some logical tag for the elements, for example "Body",
"Surface", "Boundary", etc.
- `nodes::Array{Point{spacedim}, 1}`
Array of point objects representing the nodes of the mesh, ordered consistently
with the global node numbering of `Triangulation` objects.
- `element_groups::Dict{String, Array{Int64, 1}}`
Dictionary whose keys are domain indicators like "Body", "Surface", "Boundary", etc.
The associated arrays contain the element numbers in these domains.
- `element_types::Dict{Symbol, Array{Int64, 1}}`
Dictionary whose keys are symbols representing element types, for example 
`:Quadrilateral{4}`, `Line{2}`, etc.
"""
struct Mesh{spacedim}
	elements::Array{Triangulation{P, dim, spacedim} where {P, dim}, 1}
	nodes::Array{Point{spacedim}, 1}
	element_groups::Dict{String, Array{Int64, 1}}
	element_types::Dict{Type{T} where T<:Triangulation, Array{Int64, 1}}
	function Mesh(spacedim::Int64)
		elements = Array{Triangulation{P, dim, spacedim} where {P, dim}, 1}()
		nodes = Array{Point{spacedim}, 1}()
		element_groups = Dict{String, Array{Int64, 1}}()
		element_types = Dict{Type{T} where T<:Triangulation, Array{Int64, 1}}()
		new{spacedim}(elements, nodes, element_groups, 
			element_types)
	end
end



elementTypes = Dict("vertex" => Vertex,
					"line" => Line{2},
					"line3" => Line{3},
					"triangle" => Triangle{3},
					"triangle6" => Triangle{6},
					"quad" => Quadrilateral{4},
					"quad9" => Quadrilateral{9})




"""
	AssembleQuad(mesh::Mesh, mesh_data, tagToGroup::Dict{Int64, String})
Read `mesh_data[:cells]["quad"]` and construct `Quad` elements. Push these 
elements into `mesh.elements[group_name]` where `group_name` is obtained
from `tagToGroup`.
"""
function AssembleElements(element_name::String, mesh::Mesh{spacedim}, 
	mesh_data, tagToGroup::Dict{Int64, String}) where spacedim
	element_type = elementTypes[element_name]
	mesh.element_types[element_type] = []
	for i in 1:size(mesh_data[:cells][element_name])[1]
		tag = mesh_data[:cell_data][element_name]["gmsh:physical"][i]
		group_name = tagToGroup[tag]
		node_ids = tuple(mesh_data[:cells][element_name][i,:]...)
		# gmsh uses 0-based indexing, so shift node numbers by 1
		node_ids = node_ids .+ 1
		coordinates = tuple([mesh.nodes[j] for j in node_ids]...)
		element = element_type(node_ids, coordinates)
		push!(mesh.elements, element)
		push!(mesh.element_groups[group_name], length(mesh.elements))
		push!(mesh.element_types[element_type], length(mesh.elements))
	end
end

"""
	LoadPoints(mesh::Mesh{spacedim}, mesh_data) where spacedim
Generate `Point` objects for each entry in `mesh_data[:points]` and
push into the `mesh.nodes` array.
"""
function LoadPoints(mesh::Mesh{spacedim}, mesh_data) where spacedim
	for i in 1:size(mesh_data[:points])[1]
		p = Point(tuple([c for c in mesh_data[:points][i,1:spacedim]]...))
		push!(mesh.nodes, p)
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
function loadMesh(mesh_data, spacedim::Int64)
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
		mesh.element_groups[key] = []
	end

	LoadPoints(mesh, mesh_data)

	for key in keys(mesh_data[:cells])
		AssembleElements(key, mesh, mesh_data, tagToGroup)
	end

	return mesh
end







# module geometry ends here
end
# module geometry ends here