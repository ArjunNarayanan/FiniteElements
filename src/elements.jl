module elements
using geometry, quadrature


export Basis

"""
	Basis{T <: Triangulation}
Stores the support points and basis functions associated with
a particular triangulation.
# Attributes
- `support_points::NTuple{N, Point} where N`
The `N` support points of the basis functions
- `functions::NTuple{N, Function} where N`
The `N` functions associated with this triangulation 
"""
struct Basis{Triangulation}
	support_points::NTuple{N, Point} where N
	functions::NTuple{N, Function} where N
	function Basis(::Type{Line{1}})
		p1 = Point(-1.0)
		phi1(x) = (1.0 - x[1])/2.0
		p2 = Point(1.0)
		phi2(x) = (1.0 + x[1])/2.0
		new{Line{1}}((p1,p2),
					(phi1, phi2))
	end
	function Basis(::Type{Line{2}})
		p1 = Point(-1.0)
		phi1(x) = x[1]*(x[1] - 1.0)/2.0
		p2 = Point(0.0)
		phi2(x) = (1.0 - x[1])*(1.0 + x[1])
		p3 = Point(1.0)
		phi3(x) = x[1]*(x[1] + 1.0)/2.0
		new{Line{2}}((p1,p2,p3),
					 (phi1,phi2,phi3))
	end
	function Basis(::Type{Triangle{1}})
		p1 = Point(0.0,0.0)
		phi1(x) = 1.0 - x[1] - x[2]
		p2 = Point(1.0,0.0)
		phi2(x) = x[1]
		p3 = Point(0.0,1.0)
		phi3(x) = x[2]
		new{Triangle{1}}((p1,p2,p3),
						 (phi1,phi2,phi3))
	end
	function Basis(::Type{Triangle{2}})
		p1 = Point(0.0, 0.0)
		phi1(x) = (1.0 - (x[1] + x[2]))*
				  (2*(1 - (x[1]+x[2])) - 1)
		p2 = Point(1.0, 0.0)
		phi2(x) = x[1]*(2*x[1] - 1)
		p3 = Point(0.0, 1.0)
		phi3(x) = x[2]*(2*x[2] - 1)
		p4 = Point(0.5,0.0)
		phi4(x) = 4.0*(1-(x[1]+x[2]))*x[1]
		p5 = Point(0.5,0.5)
		phi5(x) = 4.0*x[1]*x[2]
		p6 = Point(0.0,0.5)
		phi6(x) = 4.0*(1-(x[1]+x[2]))*x[2]
		new{Triangle{2}}((p1,p2,p3,p4,p5,p6),
						 (phi1,phi2,phi3,phi4,phi5,phi6))
	end
	function Basis(::Type{Quadrilateral{1}})
		p1 = Point(-1.0,-1.0)
		phi1(x) = (1. - x[1])*(1. - x[2])/4
		p2 = Point(1.0,-1.0)
		phi2(x) = (1. + x[1])*(1. - x[2])/4
		p3 = Point(1.0,1.0)
		phi3(x) = (1. + x[1])*(1. + x[2])/4
		p4 = Point(-1.0,1.0)
		phi4(x) = (1. - x[1])*(1. + x[2])/4
		new{Quadrilateral{1}}((p1,p2,p3,p4),
							  (phi1,phi2,phi3,phi4))
	end
	function Basis(::Type{Quadrilateral{2}})
		p1 = Point(-1.0,-1.0)
		phi1(x) = x[1]*(x[1] - 1.)*x[2]*(x[2] - 1.)/4.
		p2 = Point(1.0,-1.0)
		phi2(x) = x[1]*(x[1] + 1.)*x[2]*(x[2] - 1.)/4.
		p3 = Point(1.0,1.0)
		phi3(x) = x[1]*(x[1] + 1.)*x[2]*(x[2] + 1.)/4.
		p4 = Point(-1.0,1.0)
		phi4(x) = x[1]*(x[1] - 1.)*x[2]*(x[2] + 1.)/4.
		p5 = Point(0.0,-1.0)
		phi5(x) = (1. - x[1])*(1. + x[1])*x[2]*(x[2] - 1.)/2.
		p6 = Point(1.0,0.0)
		phi6(x) = x[1]*(x[1] + 1.)*(1. - x[2])*(1. + x[2])/2.
		p7 = Point(0.0,1.0)
		phi7(x) = (1. - x[1])*(1. + x[1])*x[2]*(x[2] + 1.)/2.
		p8 = Point(-1.0,0.0)
		phi8(x) = x[1]*(x[1] - 1.)*(1. - x[2])*(1. + x[2])/2.
		p9 = Point(0.0, 0.0)
		phi9(x) = (1. - x[1]^2)*(1. - x[2]^2)
		new{Quadrilateral{2}}((p1,p2,p3,p4,p5,p6,p7,p8,p9),
							(phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8,phi9))
	end
end









# module elements ends here
end
# module elements ends here