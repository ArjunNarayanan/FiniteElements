#!/usr/bin/env julia

#Start Test Script
using Test

# Run tests

println("Test Geometry")
@time @test include("test_geometry.jl")
println("Test Quadrature")
@time @test include("test_quadrature.jl")
println("Test Basis")
@time @test include("test_basis.jl")
println("Test Master")
@time @test include("test_master.jl")
println("Test Map")
@time @test include("test_mapping.jl")
println("Test Assembly")
@time @test include("test_assembly.jl")
# println("Test transformation_2d Duplicate Nodes")
# @time @test include("test_duplicate_nodes.jl")
