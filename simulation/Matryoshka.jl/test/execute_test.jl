import Matryoshka

@testset "Terminal output" begin
    Matryoshka.configuration_set(Dict{Symbol, Any}(
        :criterion => "profile",
        :interpolator => "cubic_spline",
        :segmentator => "flood_fill",
        :selector => "uniform",
        :penalizer => "curvature",
        :penalizer_init => Dict(
            :k_max => 2
        ),
        :start_points => "ng_start_points_torino2.npy",
        :valid_points => "ng_valid_points_torino2.npy",
        :logging_verbosity => 2,
        :loops => 1,
        :groups => 10,
        :cascade => [
            Dict(
                :algorithm => "matryoshka",
                :budget => 10,
                :layers => 5,
            )
        ],
        :plot => true
    ))

    Matryoshka.execute()

    @test true
end

@testset "Terminal output 2" begin
    Matryoshka.configuration_set(Dict{Symbol, Any}(
        :criterion => "jazar_model",
        :interpolator => "cubic_spline",
        :segmentator => "flood_fill",
        :selector => "uniform_distance",
        :penalizer => "count",
        :start_points => "ng_start_points_torino2.npy",
        :valid_points => "ng_valid_points_torino2.npy",
        :logging_verbosity => 2,
        :loops => 1,
        :groups => 10,
        :cascade => [
            Dict(
                :algorithm => "matryoshka",
                :budget => 10,
                :layers => 5,
            )
        ],
        :plot => true
    ))

    Matryoshka.execute()

    @test true
end

@testset "Logfile output" begin
    Matryoshka.configuration_set(Dict{Symbol, Any}(
        :criterion => "profile",
        :interpolator => "cubic_spline",
        :segmentator => "flood_fill",
        :selector => "uniform",
        :penalizer => "count",
        :start_points => "ng_start_points_torino2.npy",
        :valid_points => "ng_valid_points_torino2.npy",
        :logging_verbosity => 2,
        :loops => 1,
        :groups => 10,
        :cascade => [
            Dict(
                :algorithm => "matryoshka",
                :budget => 10,
                :layers => 5,
            )
        ],
        :criterion_init => Dict(
            :_mu => 0.2,
            :mu => 0.5
        ),
        :prefix => "execute_test",
        :plot => true
    ))

    Matryoshka.execute()

    @test true
end

@testset "Variating test" begin
    Matryoshka.configuration_set(Dict{Symbol, Any}(
        :variate => "groups",
        :criterion => "profile",
        :interpolator => "cubic_spline",
        :segmentator => "flood_fill",
        :selector => "uniform",
        :penalizer => "count",
        :start_points => "ng_start_points_torino2.npy",
        :valid_points => "ng_valid_points_torino2.npy",
        :logging_verbosity => 2,
        :loops => 2,
        :groups => [9, 11],
        :cascade => [
            Dict(
                :algorithm => "matryoshka",
                :budget => 10,
                :layers => 5,
            )
        ],
        :hold_map => true,
        :prefix => "execute_var_test",
        :plot => true
    ))

    Matryoshka.execute()

    @test true
end

@testset "Execute algorithm" begin
    #Matryoshka.execute()
end
