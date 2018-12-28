module geometry

import Base: getindex

# Export types and methods
export Triangulation, Vertex, Line, Triangle, 
		Quadrilateral, Mesh





"""
	Triangulation{N, dim}
Abstract supertype for all `dim` dimensional geometric 
triangulations space with `N` nodes.
# Examples
- 2 node `Line{2}` is a 1D linear line element.
- 3 node `Line{3}` is a 1D quadratic line element.
- 4 node `Quadrilateral{4}` is a 2D linear quadrilateral element.
"""
abstract type Triangulation{N, dim} end

"""
	Vertex <: Triangulation{1, 0}
1-node, 0 dimensional point element.
"""
struct Vertex <: Triangulation{1, 0} end



"""
	Line{N} <: Triangulation{N, 1}
1D line element with `N` nodes.
"""
struct Line{N} <: Triangulation{N, 1} end


"""
	Triangle{N} <: Triangulation{N, 2}
2D triangular element with `N` nodes.
"""
struct Triangle{N} <: Triangulation{N, 2} end



"""
	Quadrilateral{N} <: Triangulation{N, 2}
2D quadrilateral element with `N` nodes.
"""
struct Quadrilateral{N} <: Triangulation{N, 2} end



elementTypes = Dict("vertex" => Vertex,
					"line" => Line{2},
					"line3" => Line{3},
					"triangle" => Triangle{3},
					"triangle6" => Triangle{6},
					"quad" => Quadrilateral{4},
					"quad9" => Quadrilateral{9})

"""
	Mesh{spacedim}
A mesh object in `spacedim` dimensional space.
# Attributes
- `data` : Dictionary with the following keys.

- `:nodes` : An array of size `(spacedim, number_of_nodes)` giving the 
corresponding nodal coordinates.
- `:elements` : A dictionary whose keys are `T::Type{<:Triangulation{N,dim}}`. 
Querying with the element type gives an array of size 
`(number_of_nodes_in_type, number_of_elements_of_type)`. This is the connectivity
matrix for this element type.
- `:node_groups` : A dictionary whose keys are strings `group_name` representing 
the names of node groups in the mesh. Querying this dictionary with the 
corresponding name gives an `Array{Int, 1}` with the node numbers of the nodes in 
this group.
- `:element_groups` : A dictionary whose keys are strings `group_name` representing
the names of element groups in the mesh. Querying this dictionary with the corresponding
name gives an `Array{Int, 1}` with the element numbers of the elements in this group.

# Constructor
	Mesh{spacedim}(data::Dict{Symbol, Any}) where spacedim
Store the contents of `data` in a `Mesh{spacedim}` object.
"""
struct Mesh{spacedim}
	data::Dict{Symbol, Any}
	function Mesh{spacedim}(data::Dict{Symbol, Any}) where spacedim
		new{spacedim}(data)
	end
end





# module geometry ends here
end
# module geometry ends here