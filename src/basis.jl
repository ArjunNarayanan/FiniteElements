module basis

using geometry, ForwardDiff

export Basis, diameter, centroid, interpolate, gradient, evaluate

"""
	diameter(::Type{<:Line})
diameter of the reference `Line` element.
"""
function diameter(::Type{<:Line})
    return 2.0
end

"""
	centroid(::Type{<:Line})
centroid of the reference `Line` element.
"""
function centroid(::Type{<:Line})
	return [0.0]
end

"""
	centroid(::Type{<:Triangle})
centroid of the reference `Triangle` element.
"""
function centroid(::Type{<:Triangle})
	return [1.0/3.0, 1.0/3.0]
end

"""
	diameter(::Type{<:Triangle})
diameter of the reference `Triangle` element.
"""
function diameter(::Type{<:Triangle})
    return sqrt(2)
end

"""
	centroid(::Type{<:Quadrilateral})
centroid of the reference `Quadrilateral` element.
"""
function centroid(::Type{<:Quadrilateral})
	return [0.0, 0.0]
end

"""
	diameter(::Type{<:Quadrilateral})
diameter of the reference `Quadrilateral` element.
"""
function diameter(::Type{<:Quadrilateral})
    return 2*sqrt(2)
end


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
    interpolate(values::Array{Float64, 1}, xi::Array{Float64, 1},
        basis::Basis)
interpolate the nodal `values` on the master element at the coordinate `xi`
using the `basis`.
"""
function interpolate(values::Array{Float64, 1}, xi::Array{Float64, 1},
    basis::Basis)

    val = 0.0
    for I in eachindex(basis.functions)
        val += values[I]*basis.functions[I](xi)
    end
    return val
end

"""
    interpolate(values::Array{Float64, 2}, xi::Array{Float64, 1},
        basis::Basis)
interpolate the nodal `values` (treating each column as a nodal value),
on the master element at the coordinate `xi` using the `basis`.
A common application is to interpolate the point `xi` from the reference
element to the element in physical space. In this case, `values` are the
nodal coordinates of the element.
"""
function interpolate(values::Array{Float64, 2}, xi::Array{Float64, 1},
    basis::Basis)

    val = zeros(size(values)[1])
    for I in eachindex(basis.functions)
        val += values[:,I]*basis.functions[I](xi)
    end
    return val
end

"""
    interpolate(values::Array{Float64, 2}, xi::Array{Float64, 2},
        basis::Basis)
broadcast the `interpolate` operation onto the columns of `xi`.
"""
function interpolate(values::Array{Float64, 2}, xi::Array{Float64, 2},
    basis::Basis)

    val = zeros(size(values)[1], size(xi)[2])
    for J in 1:size(xi)[2]
        val[:,J] = interpolate(values, xi[:,J], basis)
    end
    return val
end

"""
    gradient(values::Array{Float64, 1}, xi::Array{Float64, 1},
        basis::Basis)
interpolate the gradient of the nodal `values` on the master element
at the coordinate `xi` using the `basis`.
"""
function gradient(values::Array{Float64, 1}, xi::Array{Float64, 1},
    basis::Basis)

    val = similar(xi)
    fill!(val, 0.0)
    for I in eachindex(basis.functions)
        val += values[I]*ForwardDiff.gradient(basis.functions[I], xi)
    end
    return val
end

"""
	outer(v1::AbstractArray, v2::AbstractArray)
outer product of `v1` and `v2` such that
`product[i,j] = v1[i]*v2[j]`
"""
function outer(v1::AbstractArray, v2::AbstractArray)
	N = length(v1)
	M = length(v2)
	product = zeros(N,M)
	for j in 1:M
		for i in 1:N
			product[i,j] = v1[i]*v2[j]
		end
	end
	return product
end

"""
	gradient(values::Array{Float64, 2}, xi::Array{Float64, 1},
		basis::Basis)
treat each column of `values` as a vector nodal value, and interpolate the
gradient of these values on the master element at the coordinate `xi` using
the given `basis`.
"""
function gradient(values::Array{Float64, 2}, xi::Array{Float64, 1},
	basis::Basis)

	num_components, num_basis = size(values)
	@assert num_basis == length(basis.support_points)
	dim = length(xi)
	grad = zeros(num_components,dim)
	for I in 1:num_basis
		basis_gradient = ForwardDiff.gradient(basis.functions[I], xi)
		grad += outer(values[:,I], basis_gradient)
	end
	return grad
end

"""
	evaluate(basis::Basis, coordinates::Array{Float64, 2})
return an array of size `(n_basis_functions, n_points)` representing the
value of the `basis` at each point in `coordinates`.
"""
function evaluate(basis::Basis, coordinates::Array{Float64, 2})
	val = zeros(length(basis.functions), size(coordinates)[2])
	for J in 1:size(coordinates)[2]
		for I in eachindex(basis.functions)
			val[I,J] = basis.functions[I](coordinates[:,J])
		end
	end
	return val
end

# module basis ends here
end
# module basis ends here
