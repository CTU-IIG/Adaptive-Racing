import Matryoshka
import NPZ: npzread
import BenchmarkTools: @btime, @benchmark
import Statistics: mean, minimum, median, maximum, std
import Printf: @printf

@testset "init" begin
    global _start_points, _valid_points

    _start_points = npzread(Matryoshka.CONFIGURATION[:start_points])
    _valid_points = npzread(Matryoshka.CONFIGURATION[:valid_points])

    Matryoshka.selector_init(Dict(:fixed_points => [[14.645512677269824, 2.543687937035088]]))
    Matryoshka.segmentator_init(_valid_points, Dict(:reserve_width => true, :range_limit => 0.5))
end

@testset "Matryoshka" begin
    _group_centers = Matryoshka.select(
        _start_points,
        Matryoshka.CONFIGURATION[:groups]
    )

    _groups = Matryoshka.segmentate(
        _valid_points,
        _group_centers
    )

    _glayers = Matryoshka.groups_border_obtain(
        _groups
    )

    _glayers = Matryoshka.groups_border_beautify(
        _glayers,
        400
    )


    _gcenters = Matryoshka.groups_center_compute(
        _groups
    )


    #b = @benchmark begin # Add $ before previously defined variables
    MATRYOSHKA = [Matryoshka.matryoshka_create(_glayers[_i], _gcenters[_i], 5) for _i in 1:length(_groups)]
    #end evals = 10

    #@printf "min: %s\nmed: %s\navg: %s\nmax: %s\nstd: %s\n" minimum(b) median(b) mean(b) maximum(b) std(b)
end
