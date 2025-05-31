using Test: @test, @testset, @test_broken, @test_logs, @test_throws

using Matryoshka
using Aqua

#@testset "Precompile" begin
#    Base.ENV["JULIA_LOAD_PATH"] = "..:" * Base.ENV["JULIA_LOAD_PATH"]
#    include("../build/generate_precompile.jl")
#end

@testset "Configuration" begin
    include("configuration_test.jl")
end

@testset "Transform" begin
    include("transform_test.jl")
end

@testset "Criterion" begin
    include("criterion_test.jl")
end

@testset "Execution" begin
    include("execute_test.jl")
end

Aqua.test_all(Matryoshka; ambiguities = false, stale_deps = (; ignore = [:Aqua]), deps_compat = false, project_extras = false)
Aqua.test_ambiguities(Matryoshka)
