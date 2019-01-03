using geometry, master, maps, reinitialize, 
		LinearAlgebra, Test, quadrature

tol = 1e-10



##################################################
# Test mapping on Triangle{3} with 3-point rule


t1 = [1.0  4.0  2.0
	  1.0  2.0  4.0]


mapping = Map{Triangle{3}, 2}(1, :coordinates, :gradients)

reinit(mapping, t1)

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
# Test edge mapping on Triangle{3} with 3 point rule

t1 = [1.0  4.0  2.0
	  1.0  2.0  4.0]

quad1D = Quadrature(Triangle{3}, 1, 3)
mapping = Map{Triangle{3},2}(quad1D, :coordinates, :gradients)
reinit(mapping, t1)

exp_coords = [ [1+3*(1-sqrt(3/5))/2, 1+(1-sqrt(3/5))/2],
			   [1+3/2, 1+ 1/2],
			   [1+3*(1+sqrt(3/5))/2, 1+(1+sqrt(3/5))/2] ]
@test norm([mapping[:coordinates][I] - exp_coords[I] for I = 1:3]) < tol






##################################################


##################################################
# Test mapping on Quadrilateral{1} with 3 point rule
master_elmt = Master(Quadrilateral{4}, 2, :values, :gradients)


q1 = [2.0  6.0  6.0   3.0
      2.0  2.0  10.0  8.0]


mapping = Map{Quadrilateral{4},2}(2, :coordinates, :gradients)

exp_coords = [[]]

reinit(mapping, q1)

@test length(mapping[:coordinates]) == 4
@test length(mapping[:jacobian]) == 4

p = [sum([master_elmt[:values][i,j]*q1[:,i] for i = 1:4]) for j in 1:4]

@test norm([norm([p[j][i] - mapping[:coordinates][j][i] for i in 1:2]) for j in 1:4]) < tol

dxdξ = [sum([ master_elmt[:gradients][i,j][1]*q1[:,i] for i in 1:4]) for j in 1:4]
dxdη = [sum([ master_elmt[:gradients][i,j][2]*q1[:,i] for i in 1:4]) for j in 1:4]

for i in 1:4
	@test abs(dxdξ[i][1] - mapping[:jacobian][i][1,1]) < tol
	@test abs(dxdξ[i][2] - mapping[:jacobian][i][2,1]) < tol
	@test abs(dxdη[i][1] - mapping[:jacobian][i][1,2]) < tol
	@test abs(dxdη[i][2] - mapping[:jacobian][i][2,2]) < tol
end

for i in 1:4
	@test abs(det(mapping[:jacobian][i]) - mapping[:determinant][i]) < tol
	@test norm(inv(mapping[:jacobian][i]) - mapping[:inverse_jacobian][i]) < tol
end

##################################################


true