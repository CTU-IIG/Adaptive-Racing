import Matryoshka

@testset "Set parameters for criterion" begin
    Matryoshka.criterion_init(
        Dict{Symbol, Any}(
            :mu => 1,
            :_mu => 0.4,
            :logfile => stdout
        )
    )

    @test Matryoshka.get_value(Matryoshka.P_cri, "mu") == 0.4
end
