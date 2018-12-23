using TensorOperations

δ = [1.0 0.0
     0.0 1.0]

λ = 100.0
μ = 50.0

E = zeros(2,2,2,2)

@tensor begin
	E[i,j,k,l] = λ*δ[i,j]*δ[k,l] + 2*μ*δ[i,k]*δ[j,l]
end


# We define a Tensor type that is parameterized by its 
# rank and dimensionality
abstract type Tensor{rank} end

# We define the trial and test functions as subtypes
# of Tensor.
mutable struct TrialFunction{rank} <: Tensor{rank} 
	der::Int
	function TrialFunction{rank}() where rank
		new{rank}(0)
	end
end

mutable struct TestFunction{rank} <: Tensor{rank} 
	der::Int
	function TestFunction{rank}() where rank
		new{rank}(0)
	end
end


# The gradient operator has two effects -- it changes
# the rank and the arg field
function ∇(u::Tensor)
	u.arg == 0 ? 1 : error("Higher basis function derivatives are currently not implemented.")
	u.arg += 1
end




macro assemble(B)
	for arg in B.args
		println(arg)
	end
end



# How do we want to specify the bilinear and linear forms?
u = TrialFunction{1}()
v = TestFunction{1}()

B = :(∇(v)[i,j]*E[i,j,k,l]*∇(u)[k,l])

@assemble ∇(v)[i,j]*E[i,j,k,l]*∇(u)[k,l]





