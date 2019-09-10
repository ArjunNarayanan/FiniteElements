module maps

import Base: getindex, values



using geometry, quadrature, master

export Map, getindex, coordinates, gradients, values



"""
	coordinates(data::Dict, nq::Int64, dim::Int64, spacedim::Int64)
Allocate arrays to store the spatial coordinates of `nq`
quadrature points in `spacedim` dimensional space.
"""
function coordinates(data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64, master_elmt::Master)
	data[:coordinates] = [zeros(spacedim) for i in 1:nq]
end



"""
	gradients(arg::Symbol, data::Dict, nq::Int64, dim::Int64, spacedim::Int64)
Allocate arrays to store the jacobian mapping, the inverse jacobian,
and the determinant of the jacobian.
After this call, `data` is modified as:
- `data[:jacobian]` - `nq` copies of zero arrays of dimension `(spacedim,dim)`
This represents the following matrix
	[∂x/∂ξ 	  ∂x/∂η    ...
	 ∂y/∂ξ    ∂y/∂η    ...
	 .
	 .
	 						]
- `data[:inverse_jacobian] - `nq` copies of zero arrays of dimension `(dim,spacedim)`
This represents the following matrix
	[∂ξ/∂x 	  ∂ξ/∂y    ...
	 ∂η/∂x    ∂η/∂y    ...
	 .
	 .
	 						]
- `data[:determinant]` - `nq` zeros representing the determinant of the corresponding
entries of `data[:jacobian]`
- `data[:dv]` - `nq` zeros to store the integrating measure `dx`.
- `data[:gradients]` - `(number_of_basis_funcs, nq)` array of `zeros(spacedim)`
to store the mapped gradients of basis functions in the physical element.
"""
function gradients(data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64, master_elmt::Master)
	data[:jacobian] = [zeros(spacedim, dim) for i in 1:nq]
	data[:inverse_jacobian] = [zeros(dim, spacedim) for i in 1:nq]
	data[:determinant] = zeros(nq)
	data[:dx] = zeros(nq)
	data[:gradients] = similar(master_elmt[:gradients])
	for i in eachindex(data[:gradients])
		data[:gradients][i] = zeros(spacedim)
	end
end


"""
	values(data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64, master_elmt::Master)
- Store `master_elmt[:values]` in `data[:values]`.
"""
function values(data::Dict, nq::Int64, dim::Int64,
					spacedim::Int64, master_elmt::Master)
	data[:values] = master_elmt[:values]
end


"""
	mapArgsToMasterArgs(args::Vararg{Symbol})
Convert the tuple of `args` supplied to `Map` into the corresponding
arguments required by the constructor of `Master`.
The arguments are mapped as:
- `:coordinates` -> `:values` since we need basis function values to map coordinates.
- `:gradients` -> `:gradients`
"""
function mapArgsToMasterArgs(args::Vararg{Symbol})
	masterArgsDict = Dict(:coordinates => :values,
						:gradients => :gradients,
						:values => :values)
	masterArgs = tuple([ haskey(masterArgsDict, arg) ? masterArgsDict[arg] : :nothing for arg in args]...)
	return masterArgs
end


"""
	Map{Triangulation, dim, spacedim}
Defines a map on a `Triangulation` from `dim` dimensions into
`spacedim` dimensions. The map stores coordinates, and
jacobian transformation information.
# Keys
Currently supported keys
- `:coordinates` Compute spatial coordinates of quadrature points upon reinitialization.
- `:values` Compute values of basis functions at quadrature points upon reinitialization.
- `:gradients` Compute gradients of basis functions at quadrature points upon reinitialization.
"""
struct Map{Triangulation, dim, spacedim}
	master::Master{T} where T
	data::Dict{Symbol, Array}
	args::NTuple{N, Symbol} where N
	function Map{T,spacedim}(quad_order::Int64,
				args::Vararg{Symbol}) where {T <: Triangulation{N,dim}} where {N,dim,spacedim}
		master_args = mapArgsToMasterArgs(args...)
		master_elmt = Master(T, quad_order, master_args...)

		data = Dict{Symbol, Array}()
		nq = length(master_elmt.quadrature.points)
		for arg in args
			eval(arg)(data, nq, dim, spacedim, master_elmt)
		end
		new{T,dim,spacedim}(master_elmt, data, args)
	end
	function Map{T,spacedim}(quadrature::Quadrature,
				args::Vararg{Symbol}) where {T <: Triangulation{N,dim}} where {N,dim,spacedim}
		master_args = mapArgsToMasterArgs(args...)
		master_elmt = Master(T, quadrature, master_args...)

		data = Dict{Symbol, Array}()
		nq = length(master_elmt.quadrature.points)
		for arg in args
			eval(arg)(data, nq, dim, spacedim, master_elmt)
		end
		new{T,dim,spacedim}(master_elmt, data, args)
	end
end


"""
	coordinates(mapping::Map{T,dim,spacedim},
	nodal_coordinates::Array{Float64, 2}) where {T,dim,spacedim}
Compute the physical coordinates of each `master_elmt.quadrature.points` on the given `element`.
Store the result in `mapping.data[:coordinates]`
"""
function coordinates(mapping::Map{T,dim,spacedim},
	nodal_coordinates::Array{Float64, 2}) where {T,dim,spacedim}
	for q in 1:length(mapping.master.quadrature.points)
		for i in 1:spacedim
			mapping.data[:coordinates][q][i] = sum([mapping.master[:values][I,q]*nodal_coordinates[i,I] for I in 1:length(mapping.master.basis.functions)])
		end
	end
end


"""
	values(mapping::Map{T,dim,spacedim},
	nodal_coordinates::Array{Float64, 2}) where {T,dim,spacedim}
Basis function values don't need to be updated.
This is just an empty function.
"""
function values(mapping::Map{T,dim,spacedim},
	nodal_coordinates::Array{Float64, 2}) where {T,dim,spacedim}

end




"""
	invert(A::Array{Float64, 2}, B::Array{Float64, 2},
	J::Float64, ::Type{<:Map{T, 2, 2}}) where T
Compute the inverse map of the jacobian. For a 2X2 matrix, this is the regular
matrix inverse.
"""
function invert(A::Array{Float64, 2},
				B::Array{Float64, 2}, J::Float64,
	 			::Type{<:Map{T, 2, 2}}) where T
	B[1,1] =  1.0/J*A[2,2]
	B[1,2] = -1.0/J*A[1,2]
	B[2,1] = -1.0/J*A[2,1]
	B[2,2] =  1.0/J*A[1,1]
end

"""
	invert(A::Array{Float64, 2}, B::Array{Float64, 2},
	J::Float64, ::Type{<:Map{T, 2, 2}}) where T
Compute the inverse map of the jacobian `A` and store it in `B`.
For a 2X1 matrix, this is the element-wise reciprocal.
"""
function invert(A::Array{Float64, 2},
				B::Array{Float64, 2}, J::Float64,
	 			::Type{<:Map{T, 1, 2}}) where T
	B[1,1] = 1.0/A[1,1]
	B[1,2] = 1.0/A[2,1]
end


"""
	determinant(A::Array{Float64, 2}, T::Type{<:Map{N,2,2}}) where N
Compute the 'determinant' of the jacobian mapping.
For a 2X2 matrix, this is the regular matrix determinant.
"""
function determinant(A::Array{Float64, 2}, ::Type{<:Map{T,2,2}}) where T
	return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end

"""
	determinant(A::Array{Float64, 2}, ::Type{<:Map{T,1,2}}) where T
Compute the 'determinant' of the jacobian mapping.
For a 2X1 matrix, this is the L2 norm.
"""
function determinant(A::Array{Float64, 2}, ::Type{<:Map{T,1,2}}) where T
	return sqrt(A[1,1]^2 + A[2,1]^2)
end



"""
	gradients(mapping::Map{T},
	nodal_coordinates::Array{Float64, 2}) where {T,dim,spacedim}
At each `master_elmt.quadrature.points`
- Compute the jacobian transformation from `master_elmt` to `element`. Store the result in `mapping.data[:jacobian]`.
- Compute the inverse of the jacobian transformation by calling `invert`.
Store the result in `mapping.data[:inverse_jacobian]`.
- Compute the determinant of the jacobian transformation by calling `determinant`.
Store the result in `mapping.data[:determinant]`.
"""
function gradients(mapping::Map{T,dim,spacedim},
	nodal_coordinates::Array{Float64, 2}) where {T<:Triangulation{N,dim}} where {N,dim,spacedim}
	for q in 1:length(mapping.master.quadrature.points)
		mapping.data[:jacobian][q] = sum([nodal_coordinates[:,I]*mapping.master[:gradients][I,q]' for I in 1:length(mapping.master.basis.functions)])
		mapping.data[:determinant][q] = determinant(mapping.data[:jacobian][q], typeof(mapping))
		invert(mapping.data[:jacobian][q],
			mapping.data[:inverse_jacobian][q],
			mapping.data[:determinant][q], typeof(mapping))
		mapping[:dx][q] = mapping[:determinant][q]*mapping.master.quadrature.weights[q]
		for I in 1:N
			mapping[:gradients][I,q] = mapping[:inverse_jacobian][q]'*mapping.master[:gradients][I,q]
		end
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
