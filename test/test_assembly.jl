using FiniteElements, Test, TensorOperations, LinearAlgebra

const dofs = 2
# Construct an arbitrary Triangle{1} object
p1 = Point(1,1)
p2 = Point(4,2)
p3 = Point(2,4)
T = Triangle{3}((1,2,3), (p1,p2,p3))


element_matrix = elementMatrix(typeof(T), dofs)
element_rhs = elementRHS(typeof(T), dofs)

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

master = Master(Triangle{3}, 2, 0, 1)
mapping = Map(master, 2, :coordinates, :derivatives)
reinit(mapping, T)

system_matrix = SystemMatrix(dofs)
system_rhs = SystemRHS(dofs)
KIJ = zeros(2,2)

for q in eachindex(master.quadrature.points)
	pq = master.quadrature.points[q]
	wq = master.quadrature.weights[q]
	for I in 1:length(T.nodes)
		∇ϕI = mapping[:inverse_jacobian][q]'*master[1][I,q]
		for J in 1:length(T.nodes)
			∇ϕJ = mapping[:inverse_jacobian][q]'*master[1][J,q]

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
master = Master(Quadrilateral{4}, 2, 0, 1)
mapping = Map(master, 2, :coordinates, :derivatives)

p1 = Point(2,2)
p2 = Point(6,2)
p3 = Point(6,10)
p4 = Point(3,8)
T = Quadrilateral{4}((1,2,3,4), (p1,p2,p3,p4))

element_matrix = elementMatrix(typeof(T), dofs)
element_rhs = elementRHS(typeof(T), dofs)

# Test that the element matrix is of correct size
@test size(element_matrix) == (4,4)
@test size(element_rhs) == (4,)

reinit(element_matrix)
reinit(element_rhs)


reinit(mapping, T)

system_matrix = SystemMatrix(dofs)
system_rhs = SystemRHS(dofs)
KIJ = zeros(2,2)

for q in eachindex(master.quadrature.points)
	pq = master.quadrature.points[q]
	wq = master.quadrature.weights[q]
	for I in 1:length(T.nodes)
		∇ϕI = mapping[:inverse_jacobian][q]'*master[1][I,q]
		for J in 1:length(T.nodes)
			∇ϕJ = mapping[:inverse_jacobian][q]'*master[1][J,q]

			fill!(KIJ, 0.0)
			
			@tensor begin
				KIJ[i,j] = ∇ϕI[i]*∇ϕJ[j]
			end

			element_matrix[I,J] += KIJ*mapping[:determinant][q]*wq
		end
		element_rhs[I] += ∇ϕI*mapping[:determinant][q]*wq
	end
end





















true