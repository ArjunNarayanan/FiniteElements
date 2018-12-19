using FiniteElements, Test

###############################################################
QLine1 = Quadrature(Line,1)
QLine2 = Quadrature(Line,2)
QLine3 = Quadrature(Line,3)
QLine4 = Quadrature(Line,4)
QLine5 = Quadrature(Line,5)
###############################################################

###############################################################
# Check quadrature on line elements

# 1-point rule should integrate first order polynomial
function P1(p::Point{1})
	return 3*p.x[1] + 5
end

tol = 1e-10


I1 = Integrate(P1, QLine1)
@test abs(I1 - 10.0) < tol

# 2-point rule should be able to integrate cubics
function P3(p::Point{1})
	x = p.x[1]
	return 4x^3 - 3x^2 + 7x +2
end

I2 = Integrate(P3, QLine2)
@test abs(I2 - 2.0) < tol


# 3-point rule should be able to integrate quintics
function P5(p::Point{1})
	x = p.x[1]
	return 12x^5 - 10x^4 + 5x^3 + 6x^2 + 11x
end

I3 = Integrate(P5, QLine3)
@test abs(I3 - 0.0) < tol

# 4-point rule should integrate 7th order polynomial
function P7(p::Point{1})
	x = p.x[1]
	return 32x^7 + 21x^6 + 18x^5 - 20x^4 - 24x^3 - 12x^2 + x
end

I4 = Integrate(P7, QLine4)
@test abs(I4 - (-10.0)) < tol

# 5-point rule should integrate 9th order polynomial
function P9(p::Point{1})
	x = p.x[1]
	return 40x^9 - 45x^8 + 32x^7 + 21x^6 + 18x^5 - 20x^4 - 24x^3 - 12x^2 + x
end

I5 = Integrate(P9, QLine5)
@test abs(I5 - (-20.0)) < tol
###############################################################


###############################################################
# Check quadrature on triangles


###############################################################
QTriangle1 = Quadrature(Triangle,1)
QTriangle2 = Quadrature(Triangle,2)
QTriangle3 = Quadrature(Triangle,3)
###############################################################


# 1-point rule over triangles should integrate
# 2D linear
function P11(p::Point{2})
	return p.x[1] + p.x[2]
end

I1 = Integrate(P11, QTriangle1)
@test abs(I1 - 1/3) < tol

# 3 point rule over triangles should integrate
# xy, x^2, y^2
function P22(p::Point{2})
	x,y = p.x
	return x*y + x^2 + y^2
end

I2 = Integrate(P22, QTriangle2)
@test abs(I2 - 5/24) < tol

# 4 point rule should integrate 
# x^3, y^3, x^2y, xy^2
function P33(p::Point{2})
	x,y = p.x
	return x^3 + y^3 + x*y^2 + x^2*y
end

I3 = Integrate(P33, QTriangle3)
@test abs(I3 - 2/15) < tol
###############################################################


###############################################################
# Check quadrature on quadrilaterals


###############################################################
QQuad1 = Quadrature(Quadrilateral, 1)
QQuad2 = Quadrature(Quadrilateral, 2)
QQuad3 = Quadrature(Quadrilateral, 3)
QQuad4 = Quadrature(Quadrilateral, 4)
QQuad5 = Quadrature(Quadrilateral, 5)
###############################################################


# 1-point rule should integrate x, y
function PQ11(p::Point{2})
	return p.x[1] + p.x[2] + 5.0
end
I1 = Integrate(PQ11, QQuad1)

@test abs(I1 - 20.0) < tol

# 2nd order rule should integrate x^3, y^3, x^2*y, x*y^2
function PQ22(p::Point{2})
	x,y = p.x
	return x^3 + y^3 + 3x^2*y + 3*x*y^2
end

I2 = Integrate(PQ22, QQuad2)

@test abs(I2 - 0.0) < tol

# 3rd order rule should integrate x^5, y^5, x^4*y, y^4*x, x^2*y^2
function PQ33(p::Point{2})
	x,y = p.x
	return x^5 + y^5 + x^4*y + x*y^4 + 9*x^2*y^2
end

I3 = Integrate(PQ33, QQuad3)
@test abs(I3 - 4.0) < tol


# 4th order rule should integrate order 7
function PQ44(p::Point{2})
	x,y = p.x
	return 15*x^4*y^2 + 15*x^2*y^4 + 9*x^2*y^2 + x^7 + y*7 + x^6*y + x*y^6
end

I4 = Integrate(PQ44, QQuad4)
@test abs(I4 - 12.0) < tol

# 5th order rule should integrate order 9
function PQ55(p::Point{2})
	x,y = p.x
	return x^9 + y^9 + 25*x^4*y^4 + 15*x^4*y^2 + 15*x^2*y^4
end

I5 = Integrate(PQ55, QQuad5)
@test abs(I5 - 12.0) < tol
###############################################################






# If all tests pass
true