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

- `:nodes` : An array of size `(number_of_nodes, spacedim)` giving the 
corresponding nodal coordinates.
- `:cells` : A dictionary whose keys are `T::Type{<:Triangulation{N,dim}}`. 
Querying with the element type gives an array of size 
`(number_of_elements_of_type, number_of_nodes_in_type)`. This is the connectivity
matrix.
- `:groups` : A dictionary whose keys are strings `physical_name` representing 
the names of physical groups in the mesh. Querying this dictionary with the 
corresponding name gives an Array of tuples `(element_type, index)`.

# Constructor
	function Mesh{spacedim}(mesh_data)
Construct a `Mesh{spacedim}` object from `mesh_data` (construct this
using the python `meshio` module).
"""
struct Mesh{spacedim}
	data::Dict{Symbol, Any}
	function Mesh{spacedim}(mesh_data) where spacedim

		data = Dict()

		data[:points] = mesh_data[:points][:,1:spacedim]

		data[:cells] = Dict{DataType, Array{Int64, 2}}()

		for key in keys(mesh_data[:cells])
			eType = elementTypes[key]
			# Shift the node ids by 1 because gmsh uses zero based
			# indexing
			data[:cells][eType] = mesh_data[:cells][key] .+ 1
		end

		tagToName = Dict{Int64, String}()

		data[:groups] = Dict{String, Array{Tuple{DataType, Int64}}}()

		for name in keys(mesh_data[:field_data])
			tag = mesh_data[:field_data][name][1]
			tagToName[tag] = name
			data[:groups][name] = Array{Tuple{DataType, Int64}, 1}()
		end


		for key in keys(mesh_data[:cell_data])
			eType = elementTypes[key]
			tags = mesh_data[:cell_data][key]["gmsh:physical"]
			for i in eachindex(tags)
				tag = tags[i]
				name = tagToName[tag]
				push!(data[:groups][name], (eType, i))
			end
		end

		new{spacedim}(data)
	end
end





# module geometry ends here
end
# module geometry ends here