#!/usr/bin/env julia
 
#Start Test Script
using geometry, quadrature
using Test
 
# Run tests
 
println("Test Geometry")
@time @test include("test_geometry.jl")
println("Test Quadrature")
@time @test include("test_quadrature.jl")
println("Test Basis")
@time @test include("test_basis.jl")