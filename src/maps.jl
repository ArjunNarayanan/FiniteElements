module maps

import Base: getindex

# We use StaticArrays to store mapping data 
# at the quadrature points. This allows to deploy
# specialized methods for inversion etc. based on 
# dimensions
using StaticArrays

using geometry, quadrature, master

export Map, reinit, getindex



"""
	coordinates(arg::Symbol, data::Dict, nq::Int64, dim::Int64, spacedim::Int64)
Allocate static arrays to store the spatial coordinates of `nq`
quadrature points given a specific instance of `Triangulation`.
"""
function coordinates(data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64)
	data[:coordinates] = [ @MVector zeros(spacedim) for i in 1:nq]
end



"""
	derivatives(arg::Symbol, data::Dict, nq::Int64, dim::Int64, spacedim::Int64)
Allocate static arrays to store the jacobian mapping, the inverse jacobian,
and the determinant of the jacobian.
After this call, `data` is modified as:
- `data[:jacobian]` - `nq` copies of zero static arrays of dimension `(spacedim,dim)`
This represents the following matrix
	[∂x/∂ξ 	  ∂x/∂η    ...
	 ∂y/∂ξ    ∂y/∂η    ...
	 .
	 .
	 						]
- `data[:inverse_jacobian] - `nq` copies of zero static arrays of dimension `(dim,spacedim)`
This represents the following matrix
	[∂ξ/∂x 	  ∂ξ/∂y    ...
	 ∂η/∂x    ∂η/∂y    ...
	 .
	 .
	 						]
- `data[:determinant] - `nq` zeros representing the determinant of the corresponding
entries of `data[:jacobian]`
"""
function derivatives(data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64)
	data[:jacobian] = [@MMatrix zeros(dim, spacedim) for i in 1:nq]
	data[:inverse_jacobian] = [@MMatrix zeros(spacedim, dim) for i in 1:nq]
	data[:determinant] = zeros(nq)
end




"""
	Map{Triangulation, dim, spacedim}
Defines a map on a `Triangulation` from `dim` dimensions into 
`spacedim` dimensions. The map stores coordinates, and
jacobian transformation information.
# Keys
Currently supported keys
- `:coordinates` allocate memory to store the spatial coordinates of 
quadrature points.
- `:derivatives` allocate memory to map derivatives of basis 
functions of quadrature points from the master element to the 
spatial element.
"""
struct Map{Triangulation, dim, spacedim}
	master::Master{T} where T
	data::Dict{Symbol, Array}
	args::NTuple{N, Symbol} where N
	function Map(master::Master{T}, spacedim::Int64,
				args::Vararg{Symbol}) where {T <: Triangulation{P,dim}}  where {P,dim}
		data = Dict{Symbol, Array}()
		nq = length(master.quadrature.points)
		for arg in args
			eval(arg)(data, nq, dim, spacedim)
		end
		new{T,dim,spacedim}(master, data, args)
	end
end


"""
	coordinates(mapping::Map{T}, master::Master{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
Compute the physical coordinates of each `master.quadrature.points` on the given `element`.
Store the result in `mapping.data[:coordinates]`
"""
function coordinates(mapping::Map{T}, master::Master{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
	for q in 1:length(master.quadrature.points)
		for i in 1:spacedim
			mapping.data[:coordinates][q][i] = sum([master[0][I,q]*element.points[I][i] for I in 1:length(master.basis.functions)])
		end
	end
end


"""
	invert(A::MArray{Tuple{2,2}}, B::MArray{Tuple{2,2}}, J::Float64)
Invert `A` and store it in `B`.
"""
function invert(A::MArray{Tuple{2,2}}, 
	B::MArray{Tuple{2,2}}, J::Float64)
	B[1,1] =  1.0/J*A[2,2]
	B[1,2] = -1.0/J*A[1,2]
	B[2,1] = -1.0/J*A[2,1]
	B[2,2] =  1.0/J*A[1,1]
end

"""
	determinant(A::MArray{Tuple{2,2}})
Return the determinant of `A`.
"""
function determinant(A::MArray{Tuple{2,2}})
	return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end


"""
	derivatives(mapping::Map{T}, master::Master{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
At each `master.quadrature.points`
- Compute the jacobian transformation from `master` to `element`. Store the result in `mapping.data[:jacobian]`.
- Compute the inverse of the jacobian transformation by calling `invert`.
Store the result in `mapping.data[:inverse_jacobian]`.
- Compute the determinant of the jacobian transformation by calling `determinant`.
Store the result in `mapping.data[:determinant]`.
"""
function derivatives(mapping::Map{T}, master::Master{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
	for q in 1:length(master.quadrature.points)
		for j in 1:spacedim
			for i in 1:dim
				mapping.data[:jacobian][q][i,j] = sum([master[1][I,q][i]*element.points[I][j] for I in 1:length(master.basis.functions)])
			end
		end
		mapping.data[:determinant][q] = determinant(mapping.data[:jacobian][q])
		invert(mapping.data[:jacobian][q], mapping.data[:inverse_jacobian][q], mapping.data[:determinant][q])
	end
end


"""
	reinit(mapping::Map{T}, master::Master{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
Reinitialize the map on the given element using the master element.
"""
function reinit(mapping::Map{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
	for arg in mapping.args
		eval(arg)(mapping, mapping.master, element)
	end
end



"""
	getindex(m::Mapping, v::Symbol)
Convenience shorthand for `m.data[v]`.
"""
function getindex(m::Map, v::Symbol)
	return m.data[v]
end







# module maps ends here
end
# module maps ends here