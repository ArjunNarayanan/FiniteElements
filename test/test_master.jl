using Test
using FiniteElements, LinearAlgebra



tol = 1e-10

function TestTriangulation(T,n,p)
	master = Master(T{n}, p, 0)
	l = length(master.basis.functions)
	b = master.basis.functions
	q = master.quadrature.points
	vals = master[0]
	@test norm([vals[i,j] - b[i](q[j]) for i in 1:l, j in 1:p]) < tol
end






TestTriangulation(Line,1,1)
TestTriangulation(Line,1,2)
TestTriangulation(Line,1,3)
TestTriangulation(Line,1,4)
TestTriangulation(Line,1,5)

TestTriangulation(Line,2,1)
TestTriangulation(Line,2,2)
TestTriangulation(Line,2,3)
TestTriangulation(Line,2,4)
TestTriangulation(Line,2,5)

TestTriangulation(Triangle,1,1)
TestTriangulation(Triangle,1,2)
TestTriangulation(Triangle,1,3)

TestTriangulation(Triangle,2,1)
TestTriangulation(Triangle,2,2)
TestTriangulation(Triangle,2,3)

TestTriangulation(Quadrilateral,1,1)
TestTriangulation(Quadrilateral,1,2)
TestTriangulation(Quadrilateral,1,3)
TestTriangulation(Quadrilateral,1,4)
TestTriangulation(Quadrilateral,1,5)

TestTriangulation(Quadrilateral,2,1)
TestTriangulation(Quadrilateral,2,2)
TestTriangulation(Quadrilateral,2,3)
TestTriangulation(Quadrilateral,2,4)
TestTriangulation(Quadrilateral,2,5)


true