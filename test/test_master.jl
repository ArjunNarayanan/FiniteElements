using Test
using geometry, master, LinearAlgebra



tol = 1e-10

function TestTriangulation(T,n,p)
	master = Master(T{n}, p, 0)
	l = length(master.basis.functions)
	b = master.basis.functions
	q = master.quadrature.points
	vals = master[0]
	@test norm([vals[i,j] - b[i](q[j]) for i in 1:l, j in 1:p]) < tol
end






TestTriangulation(Line,2,1)
TestTriangulation(Line,2,2)
TestTriangulation(Line,2,3)
TestTriangulation(Line,2,4)
TestTriangulation(Line,2,5)

TestTriangulation(Line,3,1)
TestTriangulation(Line,3,2)
TestTriangulation(Line,3,3)
TestTriangulation(Line,3,4)
TestTriangulation(Line,3,5)

TestTriangulation(Triangle,3,1)
TestTriangulation(Triangle,3,2)
TestTriangulation(Triangle,3,3)

TestTriangulation(Triangle,6,1)
TestTriangulation(Triangle,6,2)
TestTriangulation(Triangle,6,3)

TestTriangulation(Quadrilateral,4,1)
TestTriangulation(Quadrilateral,4,2)
TestTriangulation(Quadrilateral,4,3)
TestTriangulation(Quadrilateral,4,4)
TestTriangulation(Quadrilateral,4,5)

TestTriangulation(Quadrilateral,9,1)
TestTriangulation(Quadrilateral,9,2)
TestTriangulation(Quadrilateral,9,3)
TestTriangulation(Quadrilateral,9,4)
TestTriangulation(Quadrilateral,9,5)


true