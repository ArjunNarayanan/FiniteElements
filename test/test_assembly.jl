using FiniteElements, Test, TensorOperations, LinearAlgebra

const dofs = 2
# Construct an arbitrary Triangle{1} object
t1 = [1. 1
	  4 2
	  2 4]
T = Triangle{3}


element_matrix = elementMatrix(T, dofs)
element_rhs = elementRHS(T, dofs)

# Test that the element matrix is of correct size
@test size(element_matrix) == (3,3)
@test size(element_rhs) == (3,)

# Check that each node has the right sized matrix and vector
for m in element_matrix
	@test size(m) == (dofs,dofs)
end
for v in element_rhs
	@test size(v) == (dofs,)
end


# Check reinitialization
for m in element_matrix
	fill!(m, 1.0)
end
for v in element_rhs
	fill!(v, 1.0)
end
reinit(element_matrix)
reinit(element_rhs)

for m in element_matrix
	@test m == zeros(dofs,dofs)
end
for v in element_rhs
	@test v == zeros(dofs)
end


##################################################
# Test assembly for Triangle{3} element
reinit(element_matrix)
reinit(element_rhs)


mapping = Map{Triangle{3},2,2}(2, :coordinates, :gradients)
reinit(mapping, t1)

system_matrix = SystemMatrix(dofs)
system_rhs = SystemRHS(dofs)
KIJ = zeros(2,2)

for q in eachindex(mapping.master.quadrature.points)
	pq = mapping.master.quadrature.points[q]
	wq = mapping.master.quadrature.weights[q]
	for I in 1:size(t1)[1]
		∇ϕI = mapping[:inverse_jacobian][q]'*mapping.master[:gradients][I,q]
		for J in 1:size(t1)[1]
			∇ϕJ = mapping[:inverse_jacobian][q]'*mapping.master[:gradients][J,q]

			fill!(KIJ, 0.0)
			
			@tensor begin
				KIJ[i,j] = ∇ϕI[i]*∇ϕJ[j]
			end

			element_matrix[I,J] += KIJ*mapping[:determinant][q]*wq
		end
		element_rhs[I] += ∇ϕI*mapping[:determinant][q]*wq
	end
end


phi1 = [-0.25, -0.25]
phi2 = [0.375, -0.125]
phi3 = [-0.125, 0.375]

grads = [phi1, phi2, phi3]

for I = 1:3
	for J = 1:3
		KIJ = 8.0*0.5*grads[I]*grads[J]'
		@test norm(KIJ - element_matrix[I,J]) < 1e-15
	end
	FI = grads[I]*4.0
	@test norm(FI - element_rhs[I]) < 1e-15
end
##################################################




##################################################
# Test assembly for Quadrilateral{4}
T = Quadrilateral{4}

mapping = Map{T,2,2}(2, :coordinates, :gradients)

q1 = [2. 2
	  6 2
	  6 10
      3 8]


element_matrix = elementMatrix(T, dofs)
element_rhs = elementRHS(T, dofs)

# Test that the element matrix is of correct size
@test size(element_matrix) == (4,4)
@test size(element_rhs) == (4,)

reinit(element_matrix)
reinit(element_rhs)


reinit(mapping, q1)

system_matrix = SystemMatrix(dofs)
system_rhs = SystemRHS(dofs)
KIJ = zeros(2,2)

for q in eachindex(mapping.master.quadrature.points)
	pq = mapping.master.quadrature.points[q]
	wq = mapping.master.quadrature.weights[q]
	for I in 1:size(q1)[1]
		∇ϕI = mapping[:inverse_jacobian][q]'*mapping.master[:gradients][I,q]
		for J in 1:size(q1)[1]
			∇ϕJ = mapping[:inverse_jacobian][q]'*mapping.master[:gradients][J,q]

			fill!(KIJ, 0.0)
			
			@tensor begin
				KIJ[i,j] = ∇ϕI[i]*∇ϕJ[j]
			end

			element_matrix[I,J] += KIJ*mapping[:determinant][q]*wq
		end
		element_rhs[I] += ∇ϕI*mapping[:determinant][q]*wq
	end
end
##################################################




















true