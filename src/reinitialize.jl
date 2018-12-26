# Methods to reinitialize objects
module reinitialize

using geometry, maps, assembly

export reinit

"""
	reinit(my_array::Array{Array{Float64}, 2})
Reinitialize all entries to zero. Specifically for 
`element_matrix` and `element_rhs` arrays.
"""
function reinit(my_array::Array{Array{Float64}})
	for i in 1:length(my_array)
		fill!(my_array[i], 0.0)
	end
end


"""
	reinit(mapping::Map{T}, master::Master{T},
	element::Triangulation{P,dim,spacedim}) where {T <: Triangulation{P,dim}} where {P,dim,spacedim}
Reinitialize the map on the given element using the master element.
"""
function reinit(mapping::Map,
	nodal_coordinates::Array{Float64, 2})
	for arg in mapping.args
		eval(arg)(mapping, nodal_coordinates)
	end
end

"""
	reinit(assembler::Assembler)
Call `reinit` on `assembler.element_matrix` and `assembler.element_rhs`.
"""
function reinit(assembler::Assembler)
	reinit(assembler.element_matrix)
	reinit(assembler.element_rhs)
end


end