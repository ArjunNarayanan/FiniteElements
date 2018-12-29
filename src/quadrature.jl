module quadrature
using geometry, basis, ForwardDiff, LinearAlgebra

# Import base operators to be overloaded
import Base: *, getindex, eachindex

export Quadrature, Integrate



"""
	domain(T::Type{<:Line})
returns the points whose convex combination determines the domain of
integration for a `T` type triangulation.
	return [-1.0 1.0]
"""
function domain(T::Type{<:Line})
	return [-1.0 1.0]
end

"""
	domain(T::Type{<:Triangle})
returns the points whose convex combination determines the domain of
integration for a `T` type triangulation.
	return [0.0   1.0    0.0
	        0.0   0.0    1.0]
"""
function domain(T::Type{<:Triangle})
	return [0.0   1.0    0.0
	        0.0   0.0    1.0]
end

"""
	domain(T::Type{<:Quadrilateral})
returns the points whose convex combination determines the domain of
integration for a `T` type triangulation.
	return [-1.0    1.0    1.0     -1.0
		    -1.0   -1.0    1.0      1.0]
"""
function domain(T::Type{<:Quadrilateral})
	return [-1.0    1.0    1.0     -1.0
		    -1.0   -1.0    1.0      1.0]
end


"""
	Quadrature{T <: Triangulation, P}
Defines an order `P` quadrature rule associated with a `Triangulation`.
# Attributes
- `points::NTuple{N, Array{Float64, 1}} where {N}`
Tuple of `N` quadrature points in `spacedim` dimensional space.
- `weights::NTuple{N, Float64} where N`
Associated weights of the quadrature points.
"""
struct Quadrature{Triangulation, order}
	points::NTuple{N, Array{Float64, 1}} where {N}
	weights::NTuple{N, Float64} where N
	function Quadrature(::Type{<:Line}, order::Int64)
		if order == 1
			p1 = [0.0]
			w1 = 2.0
			new{Line,order}((p1,), (w1,))
		elseif order == 2
			p1 = [-1.0/sqrt(3)]
			w1 = 1.0
			p2 = [1.0/sqrt(3)]
			w2 = 1.0
			new{Line,order}((p1,p2),
						(w1,w2))
		elseif order == 3
			p1 = [-sqrt(3/5)]
			w1 = 5.0/9.0
			p2 = [0.0]
			w2 = 8.0/9.0
			p3 = [sqrt(3/5)]
			w3 = 5.0/9.0
			new{Line,order}((p1,p2,p3),
						(w1,w2,w3))
		elseif order == 4
			p1 = [sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0))]
			w1 = (18.0 + sqrt(30))/36
			p2 = [-sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0))]
			w2 = (18.0 + sqrt(30))/36
			p3 = [sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0))]
			w3 = (18.0 - sqrt(30))/36
			p4 = [-sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0))]
			w4 = (18.0 - sqrt(30))/36
			new{Line,order}((p1,p2,p3,p4),
						(w1,w2,w3,w4))
		elseif order == 5
			p1 = [0.0]
			w1 = 128/225
			p2 = [1/3*sqrt(5 - 2*sqrt(10/7))]
			w2 = (322 + 13*sqrt(70))/900
			p3 = [-1/3*sqrt(5 - 2*sqrt(10/7))]
			w3 = (322 + 13*sqrt(70))/900
			p4 = [1/3*sqrt(5 + 2*sqrt(10/7))]
			w4 = (322 - 13*sqrt(70))/900
			p5 = [-1/3*sqrt(5 + 2*sqrt(10/7))]
			w5 = (322 - 13*sqrt(70))/900
			new{Line,order}((p1,p2,p3,p4,p5),
						(w1,w2,w3,w4,w5))
		else
			error("The requested quadrature order has not been implemented")
		end
	end
	function Quadrature(::Type{<:Triangle}, order::Int64)
		if order == 1
			p1 = [1/3, 1/3]
			w1 = 0.5
			new{Triangle,order}((p1,), (w1,))
		elseif order == 2
			p1 = [1/6, 1/6]
			w1 = 1/6
			p2 = [2/3, 1/6]
			w2 = 1/6
			p3 = [1/6, 2/3]
			w3 = 1/6
			new{Triangle,order}((p1,p2,p3),
							(w1,w2,w3))
		elseif order == 3
			p1 = [1/3,1/3]
			w1 = -27/96
			p2 = [1/5,1/5]
			w2 = 25/96
			p3 = [1/5,3/5]
			w3 = 25/96
			p4 = [3/5,1/5]
			w4 = 25/96
			new{Triangle,order}((p1,p2,p3,p4),
								(w1,w2,w3,w4))
		else
			error("The requested quadrature order has not been implemented")
		end
	end
"""
	Quadrature(::Type{<:Triangle{3}}, dim::Int64, order::Int64)
construct a quadrature of `dim = 1 < 2`. Quadrature is defined on only 
one edge of the `Triangle`, i.e. edge 1-2. Ensure that spatial elements
respect this correspondence while integrating.
"""
	function Quadrature(T::Type{<:Triangle{3}}, 
						dim::Int64, order::Int64)
		@assert dim == 1 "Expected dim = 1, got dim = "*string(dim)
		quad1D = Quadrature(Line{2}, order)
		basis = Basis(Line{2})
		nodes = domain(T)
		points = Array{Array{Float64, 1}, 1}()
		weights = Array{Float64, 1}()
		for i in eachindex(quad1D.points)
			p = quad1D.points[i]
			w = quad1D.weights[i]

			p_new = sum([ nodes[:,I]*basis.functions[I](p) for I in 1:2 ])
			dx = sum([nodes[:,I]*(ForwardDiff.gradient(basis.functions[I], p)[1]) for I in 1:2])
			w_new = w*norm(dx)

			push!(points, p_new)
			push!(weights, w_new)
		end
		points = tuple(points...)
		weights = tuple(weights...)
		new{Triangle, order}(points, weights)
	end
	function Quadrature(::Type{<:Quadrilateral}, order::Int64)
		q1 = Quadrature(Line,order)
		points = Array{Array{Float64, 1}, 1}()
		weights = Array{Float64, 1}()
		for i in 1:length(q1.points)
			p1 = q1.points[i]
			w1 = q1.weights[i]
			for j in 1:length(q1.points)
				p2 = q1.points[j]
				w2 = q1.weights[j]
				p_new = [p1[1], p2[1]]
				w_new = w1*w2
				push!(points, p_new)
				push!(weights, w_new)
			end
		end
		points = tuple(points...)
		weights = tuple(weights...)
		new{Quadrilateral, order}(points, weights)
	end
end


"""
	getindex(quad::Quadrature, i::Int64)
Return a tuple `(p,w)` where `p` and `w` are the `i`th
quadrature point and associated weight.
"""
function getindex(quad::Quadrature, i::Int64)
	return (quad.points[i], quad.weights[i])
end

"""
	eachindex(quad::Quadrature)
Return an iterator over `quad.points`.
"""
function eachindex(quad::Quadrature)
	return eachindex(quad.points)
end



"""
	Integrate(f::Function, q::Quadrature)
Integrate the function `f` using the quadrature rule `q`.
`f` is a function of ß`Point{spacedim}` where `spacedim` is 
consistent with the type of points in `q`.
"""
function Integrate(f::Function, q::Quadrature)
	I = 0.0
	for i in eachindex(q.points)
		I += f(q.points[i])*q.weights[i]
	end
	return I
end



# module quadrature ends here
end
# module quadrature ends here