# We define a Tensor type that is parameterized by its 
# rank and dimensionality
abstract type Tensor{rank, dim} end

# We define the trial and test functions as subtypes
# of Tensor.
struct TrialFunction{rank, dim} <: Tensor{rank, dim} 
	brank::Int64
	arg::Symbol
	function TrialFunction{rank, dim}() where {rank,dim}
		new{rank, dim}(0, :values)
	end
	function TrialFunction{rank, dim}(brank, arg) where {rank, dim}
		new{rank, dim}(brank, arg)
	end
end

struct TestFunction{rank, dim} <: Tensor{rank, dim} 
	brank::Int64
	arg::Symbol
	function TestFunction{rank, dim}() where {rank,dim}
		new{rank, dim}(0, :values)
	end
	function TestFunction{rank, dim}(brank, arg) where {rank, dim}
		new{rank, dim}(brank, arg)
	end
end


# The gradient operator has two effects -- it changes
# the rank and the arg field
function gradient(u::TrialFunction{rank,dim}) where {rank, dim}
	arg = u.arg == :values ? :gradients : error("Higher basis function derivatives are currently not implemented.")
	return TrialFunction{rank,dim}(u.brank+1, arg)
end

function gradient(u::TestFunction{rank,dim}) where {rank, dim}
	arg = u.arg == :values ? :gradients : error("Higher basis function derivatives are currently not implemented.")
	return TestFunction{rank,dim}(u.brank+1, arg)
end



# The divergence operator has a single effect. It lowers
# the rank 
function divergence(u::TrialFunction{rank, dim}) where {rank, dim}
	return TrialFunction{rank-1,dim}(u.coeff, u.arg)
end

function divergence(u::TestFunction{rank, dim}) where {rank, dim}
	return TestFunction{rank-1,dim}(u.coeff, u.arg)
end

function dot(u::Tensor{rank, dim}, v::Tensor{rank, dim}) where {rank, dim}
	return 
end




# How do we want to specify the bilinear and linear forms?
u = TrialFunction{1, 2}()
v = TestFunction{1, 2}()

B = :(dot(gradient(v), gradient(u)))



