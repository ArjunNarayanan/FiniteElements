using FiniteElements, LinearAlgebra, Test

tol = 1e-10



##################################################
# Test mapping on Triangle{1} with 1-point rule
master = Master(Triangle{3}, 1, 0, 1)

p1 = Point(1,1)
p2 = Point(4,2)
p3 = Point(2,4)
t1 = Triangle{3}((1,2,3), (p1,p2,p3))

mapping = Map(master, 2, :coordinates, :derivatives)

maps.reinit(mapping, t1)

exp_coords = [7/3, 7/3]
exp_jac = [3.0 1.0
		   1.0 3.0]

@test length(mapping[:coordinates]) == 1
@test length(mapping[:jacobian]) == 1
@test norm(mapping[:coordinates][1] - exp_coords) < tol
@test norm(mapping[:jacobian][1] - exp_jac) < tol
@test norm(mapping[:inverse_jacobian][1] - inv(exp_jac)) < tol
@test abs(mapping[:determinant][1] - det(exp_jac)) < tol
##################################################


##################################################
# Test mapping on Quadrilateral{1} with 3 point rule
master = Master(Quadrilateral{4}, 2, 0, 1)

p1 = Point(2,2)
p2 = Point(6,2)
p3 = Point(6,10)
p4 = Point(3,8)
q1 = Quadrilateral{4}((1,2,3,4), (p1,p2,p3,p4))

mapping = Map(master, 2, :coordinates, :derivatives)

exp_coords = [[]]

maps.reinit(mapping, q1)

@test length(mapping[:coordinates]) == 4
@test length(mapping[:jacobian]) == 4

p = [sum([master[0][i,j]*q1.points[i] for i = 1:4]) for j in 1:4]

@test norm([norm([p[j][i] - mapping[:coordinates][j][i] for i in 1:2]) for j in 1:4]) < tol

dxdξ = [sum([ master[1][i,j][1]*q1.points[i] for i in 1:4]) for j in 1:4]
dxdη = [sum([ master[1][i,j][2]*q1.points[i] for i in 1:4]) for j in 1:4]

for i in 1:4
	@test abs(dxdξ[i][1] - mapping[:jacobian][i][1,1]) < tol
	@test abs(dxdξ[i][2] - mapping[:jacobian][i][1,2]) < tol
	@test abs(dxdη[i][1] - mapping[:jacobian][i][2,1]) < tol
	@test abs(dxdη[i][2] - mapping[:jacobian][i][2,2]) < tol
end

for i in 1:4
	@test abs(det(mapping[:jacobian][i]) - mapping[:determinant][i]) < tol
	@test norm(inv(mapping[:jacobian][i]) - mapping[:inverse_jacobian][i]) < tol
end

##################################################


true