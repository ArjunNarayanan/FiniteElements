module master

import Base: getindex

using geometry, quadrature, ForwardDiff


export Basis, Master, getindex

"""
	Basis{T <: Triangulation}
Stores the support points and basis functions associated with
a particular triangulation.
# Attributes
- `support_points::NTuple{N, Array{Float64, 1}} where N`
The `N` support points of the basis functions
- `functions::NTuple{N, Function} where N`
The `N` functions associated with this triangulation 
"""
struct Basis{Triangulation}
	support_points::NTuple{N, Array{Float64, 1}} where N
	functions::NTuple{N, Function} where N
	function Basis(T::Type{Line{2}})
		p1 = [-1.0]
		phi1(x) = (1.0 - x[1])/2.0
		p2 = [1.0]
		phi2(x) = (1.0 + x[1])/2.0
		new{T}((p1,p2),
					(phi1, phi2))
	end
	function Basis(T::Type{Line{3}})
		p1 = [-1.0]
		phi1(x) = x[1]*(x[1] - 1.0)/2.0
		p2 = [0.0]
		phi2(x) = (1.0 - x[1])*(1.0 + x[1])
		p3 = [1.0]
		phi3(x) = x[1]*(x[1] + 1.0)/2.0
		new{T}((p1,p2,p3),
					 (phi1,phi2,phi3))
	end
	function Basis(T::Type{Triangle{3}})
		p1 = [0.0,0.0]
		phi1(x) = 1.0 - x[1] - x[2]
		p2 = [1.0,0.0]
		phi2(x) = x[1]
		p3 = [0.0,1.0]
		phi3(x) = x[2]
		new{T}((p1,p2,p3),
						 (phi1,phi2,phi3))
	end
	function Basis(T::Type{Triangle{6}})
		p1 = [0.0, 0.0]
		phi1(x) = (1.0 - (x[1] + x[2]))*
				  (2*(1 - (x[1]+x[2])) - 1)
		p2 = [1.0, 0.0]
		phi2(x) = x[1]*(2*x[1] - 1)
		p3 = [0.0, 1.0]
		phi3(x) = x[2]*(2*x[2] - 1)
		p4 = [0.5,0.0]
		phi4(x) = 4.0*(1-(x[1]+x[2]))*x[1]
		p5 = [0.5,0.5]
		phi5(x) = 4.0*x[1]*x[2]
		p6 = [0.0,0.5]
		phi6(x) = 4.0*(1-(x[1]+x[2]))*x[2]
		new{T}((p1,p2,p3,p4,p5,p6),
						 (phi1,phi2,phi3,phi4,phi5,phi6))
	end
	function Basis(T::Type{Quadrilateral{4}})
		p1 = [-1.0,-1.0]
		phi1(x) = (1. - x[1])*(1. - x[2])/4
		p2 = [1.0,-1.0]
		phi2(x) = (1. + x[1])*(1. - x[2])/4
		p3 = [1.0,1.0]
		phi3(x) = (1. + x[1])*(1. + x[2])/4
		p4 = [-1.0,1.0]
		phi4(x) = (1. - x[1])*(1. + x[2])/4
		new{T}((p1,p2,p3,p4),
							  (phi1,phi2,phi3,phi4))
	end
	function Basis(T::Type{Quadrilateral{9}})
		p1 = [-1.0,-1.0]
		phi1(x) = x[1]*(x[1] - 1.)*x[2]*(x[2] - 1.)/4.
		p2 = [1.0,-1.0]
		phi2(x) = x[1]*(x[1] + 1.)*x[2]*(x[2] - 1.)/4.
		p3 = [1.0,1.0]
		phi3(x) = x[1]*(x[1] + 1.)*x[2]*(x[2] + 1.)/4.
		p4 = [-1.0,1.0]
		phi4(x) = x[1]*(x[1] - 1.)*x[2]*(x[2] + 1.)/4.
		p5 = [0.0,-1.0]
		phi5(x) = (1. - x[1])*(1. + x[1])*x[2]*(x[2] - 1.)/2.
		p6 = [1.0,0.0]
		phi6(x) = x[1]*(x[1] + 1.)*(1. - x[2])*(1. + x[2])/2.
		p7 = [0.0,1.0]
		phi7(x) = (1. - x[1])*(1. + x[1])*x[2]*(x[2] + 1.)/2.
		p8 = [-1.0,0.0]
		phi8(x) = x[1]*(x[1] - 1.)*(1. - x[2])*(1. + x[2])/2.
		p9 = [0.0, 0.0]
		phi9(x) = (1. - x[1]^2)*(1. - x[2]^2)
		new{T}((p1,p2,p3,p4,p5,p6,p7,p8,p9),
							(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9))
	end
end


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
			grads[j,i] = ForwardDiff.gradient(basis.functions[j], [x for x in p])
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
	quadrature::Quadrature{P} where {P >: T}
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