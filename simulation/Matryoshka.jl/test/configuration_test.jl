import Matryoshka
using JSON

f_helper(x) = x
f_helper(x::Vector) = [f_helper(e) for e in x]
f_helper(d::Dict) = Dict(Symbol(k) => f_helper(v) for (k, v) in d)
symbol_dict(d::Dict) = f_helper(d)

CONFIGURATION = nothing

@testset "Load configuration" begin
    CONFIGURATION = JSON.parsefile("configuration.json")

    @test CONFIGURATION["groups"] == 10

    CONFIGURATION = symbol_dict(CONFIGURATION)

    @test CONFIGURATION[:groups] == 10
end


@testset "Load internally" begin
    @test Matryoshka.configuration_load("configuration.json")

    @test get(Matryoshka.CONFIGURATION, :groups, 0) == 10
end
