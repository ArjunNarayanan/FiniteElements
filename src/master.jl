module master

import Base: getindex

using geometry, basis, quadrature, ForwardDiff


export Basis, Master, getindex




"""
	values(basis::Basis, quadrature::Quadrature)
Compute the value of each `basis.functions` at each `quadrature.points`
and store it in an array of size
`(n_basis_functions, n_quadrature_points)`.
Return array.
"""
function values(basis::Basis, quadrature::Quadrature)
	vals = zeros(length(basis.functions), length(quadrature.points))
	for i in eachindex(quadrature.points)
		p = quadrature.points[i]
		for j in eachindex(basis.functions)
			vals[j,i] = basis.functions[j](p)
		end
	end
	return vals
end


"""
	gradients(basis::Basis, quadrature::Quadrature)
Compute the gradient of each `basis.functions` at each `quadrature.points`
store it in an array of size
`(n_basis_functions, n_quadrature_points)`.
Return array.
"""
function gradients(basis::Basis, quadrature::Quadrature)
	n = length(basis.functions)
	p = length(quadrature.points)
	grads = Array{Array{Float64, 1}, 2}(undef, n, p)
	for i in eachindex(quadrature.points)
		p = quadrature.points[i]
		for j in eachindex(basis.functions)
			grads[j,i] = ForwardDiff.gradient(basis.functions[j], p)
		end
	end
	return grads
end



"""
	Master{T <: Triangulation}
A type that stores all required basis function information like values,
gradients, etc.
# Attributes
	basis::Basis{T}
The associated `Basis` object.
	quadrature::Quadrature{T >: Triangulation}
The associated `Quadrature` object.
	data::Dict{Int, Array}
The required basis function derivative order evaluated at the quadrature points.
# Constructor
	Master(T::Type{<:Triangulation}, order::Int64, args::Vararg{Int})
`T` is the `Triangulation` object on which the basis and quadrature rules
are to be computed.
`order` is the order of quadrature rule to be used.
	args: The order of the basis function derivative to compute
- `:values` - compute basis function values
- `:gradients` - compute basis function gradients
"""
struct Master{T <: Triangulation}
	basis::Basis{T}
	quadrature::Quadrature
	data::Dict{Symbol, Array}
	function Master(T::Type{<:Triangulation}, order::Int64,
		args::Vararg{Symbol})
		basis = Basis(T)
		quadrature = Quadrature(T, order)
		data = Dict{Symbol, Array}()

		for arg in args
			if eval(arg) != nothing
				data[arg] = eval(arg)(basis, quadrature)
			end
		end

		new{T}(basis, quadrature, data)
	end
	function Master(T::Type{<:Triangulation},
		quadrature::Quadrature, args::Vararg{Symbol})
		basis = Basis(T)
		data = Dict{Symbol, Array}()

		for arg in args
			if eval(arg) != nothing
				data[arg] = eval(arg)(basis, quadrature)
			end
		end
		new{T}(basis, quadrature, data)
	end
	"""
		Master(surfaceElement::Type{<:Triangulation{N,2}},
			quad::Quadrilateral{Line,p}) where {N,p}
	special constructor for performing 1D quadrature on a 2D domain.
	This master element is used for performing 1D integration over an evolving
	interface, and is expected to be reinitialized several times. Therefore
	the struct is initialized to zero values.
	"""
	function Master(surfaceElement::Type{<:Triangulation{N,2}},
		quadrature::Quadrature{Line,p}, args::Vararg{Symbol}) where {N,p}
		basis = Basis(surfaceElement)
		data = Dict{Symbol, Array}()

		# for arg in args
		# 	if eval(arg) != nothing
		# 		data[arg] = eval(arg)(basis, quadrature)
		# 	end
		# end
		# new{T}(basis, quadrature, data)
	end
end


"""
	getindex(m::Master, v::Int)
Convenience shorthand for `m.data[v]`.
"""
function getindex(m::Master, v::Symbol)
	return m.data[v]
end




# module elements ends here
end
# module elements ends here
