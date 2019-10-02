module basis

using geometry

export Basis, diameter

"""
	diameter(::Type{<:Triangle})
diameter of the reference `Triangle` element.
"""
function diameter(::Type{<:Triangle})
    return sqrt(2)
end

"""
	diameter(::Type{<:Quadrilateral})
diameter of the reference `Quadrilateral` element.
"""
function diameter(::Type{<:Quadrilateral})
    return 2*sqrt(2)
end

"""
	diameter(::Type{<:Line})
diameter of the reference `Line` element.
"""
function diameter(::Type{<:Line})
    return 2.0
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



















# module basis ends here
end
# module basis ends here
