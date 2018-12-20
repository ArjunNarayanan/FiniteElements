using FiniteElements, LinearAlgebra, Test

tol = 1e-10



##################################################
# Test mapping on Triangle{1} with 1-point rule
master = Master(Triangle{1}, 1, :values, :gradients)

p1 = Point(1,1)
p2 = Point(4,2)
p3 = Point(2,4)
t1 = Triangle((1,2,3), (p1,p2,p3))

mapping = Map(master, Triangle{1,2}, :coordinates, :derivatives)

reinit(mapping, master, t1)

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
master = Master(Quadrilateral{1}, 2, :values, :gradients)

p1 = Point(2,2)
p2 = Point(6,2)
p3 = Point(6,10)
p4 = Point(3,8)

##################################################


true