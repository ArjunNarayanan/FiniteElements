module quadrature
using geometry

# Import base operators to be overloaded
import Base: *

export Quadrature, QLine1, QLine2, QLine3, QLine4, QLine5,
		QTriangle1, QTriangle3, QTriangle4,
		QQuad1, QQuad2, QQuad3, QQuad4, QQuad5

"""
	Quadrature{T <: Triangulation,N}
Defines a quadrature rule associated with a `Triangulation`
object with `N` points and associated weights.
# Attributes
- `points::NTuple{N, Point{spacedim}} where spacedim`
Tuple of `N` quadrature points in `spacedim` dimensional space.
- `weights::NTuple{N, Float64}`
Associated weights of the quadrature points.
"""
struct Quadrature{T <: Triangulation, N}
	points::NTuple{N, Point{spacedim}} where {spacedim}
	weights::NTuple{N, Float64}
end


###############################################################
# Defining specific 1D quadrature rules over lines
p1 = Point(0.0)
w1 = 2.0
QLine1 = Quadrature{Line,1}((p1,), (w1,))

p1 = Point(-1.0/sqrt(3))
w1 = 1.0
p2 = Point(1.0/sqrt(3))
w2 = 1.0
QLine2 = Quadrature{Line,2}((p1,p2),
					(w1,w2))


p1 = Point(-sqrt(3/5))
w1 = 5.0/9.0
p2 = Point(0.0)
w2 = 8.0/9.0
p3 = Point(sqrt(3/5))
w3 = 5.0/9.0
QLine3 = Quadrature{Line,3}((p1,p2,p3),
				    (w1,w2,w3))

p1 = Point(sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)))
w1 = (18.0 + sqrt(30))/36
p2 = Point(-sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)))
w2 = (18.0 + sqrt(30))/36
p3 = Point(sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)))
w3 = (18.0 - sqrt(30))/36
p4 = Point(-sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)))
w4 = (18.0 - sqrt(30))/36
QLine4 = Quadrature{Line,4}((p1,p2,p3,p4),
					(w1,w2,w3,w4))


p1 = Point(0.0)
w1 = 128/225
p2 = Point(1/3*sqrt(5 - 2*sqrt(10/7)))
w2 = (322 + 13*sqrt(70))/900
p3 = Point(-1/3*sqrt(5 - 2*sqrt(10/7)))
w3 = (322 + 13*sqrt(70))/900
p4 = Point(1/3*sqrt(5 + 2*sqrt(10/7)))
w4 = (322 - 13*sqrt(70))/900
p5 = Point(-1/3*sqrt(5 + 2*sqrt(10/7)))
w5 = (322 - 13*sqrt(70))/900
QLine5 = Quadrature{Line,5}((p1,p2,p3,p4,p5),
	 				(w1,w2,w3,w4,w5))
###############################################################



###############################################################
# Defining specific quadrature rules over triangles
p1 = Point(1/3, 1/3)
w1 = 1.0
QTriangle1 = Quadrature{Triangle,1}((p1,), (w1,))

p1 = Point(1/6, 1/6)
w1 = 1/3
p2 = Point(2/3, 1/6)
w2 = 1/3
p3 = Point(1/6, 2/3)
w3 = 1/3
QTriangle3 = Quadrature{Triangle,3}((p1,p2,p3),
									(w1,w2,w3))

p1 = Point(1/3,1/3)
w1 = -27/48
p2 = Point(1/5,1/5)
w2 = 25/48
p3 = Point(1/5,3/5)
w3 = 25/48
p4 = Point(3/5,1/5)
w4 = 25/48

QTriangle4 = Quadrature{Triangle,4}((p1,p2,p3,p4),
									(w1,w2,w3,w4))
###############################################################


###############################################################
# Quadrature rules for Quadrilaterals are easily obtained
# by tensor product construction from lines.

"""
	*(q1::Quadrature{Line}, q2::Quadrature{Line})
Return the tensor product of the two quadrature rules.
The result will be of type 
`Quadrature{Quadrilateral,N} where N`
"""
function *(q1::Quadrature{Line}, q2::Quadrature{Line})
	# compute the number of points in the new quadrature
	# rule
	N = length(q1.points)*length(q2.points)
	points = Array{Point{2}, 1}()
	weights = Array{Float64, 1}()
	for i in 1:length(q1.points)
		p1 = q1.points[i]
		w1 = q1.weights[i]
		for j in 1:length(q2.points)
			p2 = q2.points[j]
			w2 = q2.weights[j]
			p_new = Point(p1.x[1], p2.x[1])
			w_new = w1*w2
			push!(points, p_new)
			push!(weights, w_new)
		end
	end
	points = tuple(points...)
	weights = tuple(weights...)
	return Quadrature{Quadrilateral, N}(points, weights)
end


QQuad1 = QLine1*QLine1
QQuad2 = QLine2*QLine2
QQuad3 = QLine3*QLine3
QQuad4 = QLine4*QLine4
QQuad5 = QLine5*QLine5
###############################################################




# module quadrature ends here
end
# module quadrature ends here