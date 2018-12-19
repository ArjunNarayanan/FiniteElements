module maps

# We use Static Arrays to be able to use multiple
# dispatch on the size of arrays
using StaticArrays

using geometry, quadrature, elements

export Map



"""
	coordinates(arg::Symbol, data::Dict, nq::Int64, dim::Int64, spacedim::Int64)
Allocate static arrays to store the spatial coordinates of nq
quadrature points.
"""
function coordinates(arg::Symbol, data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64)
	data[arg] = [@SVector zeros(spacedim) for i in 1:nq]
end



"""
	Map{Triangulation, dim, spacedim}
Defines a map on a `Triangulation` from `dim` dimensions into 
`spacedim` dimensions. The map stores coordinates, and
jacobian transformation information.
# Keys
Currently supported keys
- `:coordinates` allocate memory to map the coordinates of 
quadrature points from the master element to the spatial element.
- `:derivatives` allocate memory to map derivatives of basis 
functions of quadrature points from the master element to the 
spatial element.
"""
struct Map{Triangulation, dim, spacedim}
	data::Dict{Symbol, Array}
	function Map(m::Master{T},
				::Type{<:Triangulation{P,dim,spacedim}},
				args::Vararg{Symbol}) where {T <: Triangulation{P,dim}}  where {P,dim,spacedim}
		data = Dict{Symbol, Array}()
		nq = length(m.quadrature.points)
		for arg in args
			eval(arg)(arg, data, nq, dim, spacedim)
		end
		new{T,dim,spacedim}(data)
	end
end



# module maps ends here
end
# module maps ends here