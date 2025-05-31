import Printf: @printf


######################
# Utilities (segmentator)
######################

"""
    grid_compute(points::Matrix{Float64})::Float64
"""
function grid_compute(points::Matrix{Float64})::Float64
    minimum(minimum(diff(u) for u in [unique!(sort(c)) for c in eachcol(points)]))
end


"""
    map_create(points::Array{Float64,2}, origin=nothing, size=nothing, grid=nothing)::Tuple{Matrix{UInt8}, Vector{Float64}, Float64}
"""
function map_create(points::Array{Float64,2}, origin=nothing, size=nothing, grid=nothing)::Tuple{Matrix{UInt8}, Vector{Float64}, Float64}
    global MAP, MAP_ORIGIN, MAP_GRID, TRACK_BITMAP

    println("Creating map...")

    # Obtain grid size if not set
    _grid = grid !== nothing ? grid : grid_compute(points)
    println("\tGrid: ", _grid)

    # Obtain origin if not set
    _origin = origin !== nothing ? origin : vec(reshape(minimum(points, dims=1), (2, 1)))
    println("\tOrigin: ", _origin)

    # Obtain size if not set
    _size = size !== nothing ? size : vec(reshape(map(abs, maximum(points, dims=1) - minimum(points, dims=1)), (2, 1)))

    TRACK_BITMAP = TrackBitmap([minimum(points, dims = 1)...],
                               [maximum(points, dims = 1)...],
                               [_grid,_grid])

    println("\tMin: ", vec(reshape(minimum(points, dims=1), (2, 1))))
    println("\tMax: ", vec(reshape(maximum(points, dims=1), (2, 1))))
    println("\tDist: ", vec(reshape(maximum(points, dims=1), (2, 1)) - _origin))

    println("\tSize: ", _size)
    println("\tCell size: ", (_size ./ _grid) .+ 1, " ", vec(Int.(round.((_size ./ _grid) .+ 1))))

    _m = zeros(UInt8, Tuple(convert.(UInt64, round.((_size ./ _grid) .+ 1))))

    for _p in eachrow(points)
        index = Int.(round.((_p[1:2] - _origin) ./ _grid) .+ 1)
        _m[index[1], index[2]] = 100
        TRACK_BITMAP[_p] = true;
    end

    #=
    # req: using Gnuplot, Colors
    let tb = TRACK_BITMAP
        img = map(tb.bitmap') do p p ? colorant"#eed" : colorant"white" end
        @gp :map Gnuplot.palette(:gray1) "set size ratio -1" :-
        @gp :- :map recipe(img, "dx=$(tb.p_step[1]) dy=$(tb.p_step[2]) origin=($(tb.p_min[1]), $(tb.p_min[2]))")
    end
    =#

    MAP = _m
    MAP_ORIGIN = _origin
    MAP_GRID = _grid

    println("Map created.")

    return MAP, MAP_ORIGIN, MAP_GRID
end


"""
    point_to_map(points::Vector{Float64})::Vector{Int64}
"""
function point_to_map(points::Vector{Float64})::Vector{Int64}
    # global MAP_ORIGIN, MAP_GRID
    Int.(round.((points[1:2] - MAP_ORIGIN) ./ MAP_GRID) .+ 1)
end


"""
    points_to_map(points::Matrix{Float64})::Matrix{Int64}
"""
function points_to_map(points::Matrix{Float64})::Matrix{Int64}
    # global MAP_ORIGIN, MAP_GRID
    Int.(round.((points[:, 1:2] .- MAP_ORIGIN') ./ MAP_GRID) .+ 1)
end


######################
# Utilities (selector)
######################

"""
    point_distance(a::Vector{Float64}, b::Vector{Float64})::Float64
"""
function point_distance(a::Vector{Float64}, b::Vector{Float64})::Float64
    sqrt(sum([(b[i] - a[i])^2 for i in 1:min(length(a), length(b))]))
end


"""
    points_distance(points::Matrix{Float64})::Matrix{Float64}
"""
function points_distance(points::Matrix{Float64})::Matrix{Float64} # Should be vector?
    sqrt.(sum((circshift(points[:, 1:2], 1) .- points[:, 1:2]) .^ 2, dims=2))
end


"""
    trajectory_closest_index(points::Matrix{Float64}, reference::Vector{Float64}; from_left::Bool=false)::Int64
"""
function trajectory_closest_index(points::Matrix{Float64}, reference::Vector{Float64}; from_left::Bool=false)::Int64

    _distances = points[:, 1:2] .- reference[1:2]'
    index = argmin(hypot.(_distances[:, 1], _distances[:, 2]), dims=1)[1]

    if from_left == false
        return index
    else
        d1 = hypot(_distances[index, 1], distances[index, 2])
        d2 = hypot(_distances[index+1, 1], distances[index+1, 2])
        ds = point_distance(points[index, 1:2], points[index+1, 1:2])
        return (d1^2 - d2^2 + ds^2) / (2 * d1 * ds) > 0 ? index : index - 1
    end
end

######################
# Utilities (optimizer)
######################

"""
    trajectory_reduce(points::Matrix{Float64}, remain::Int64)::Matrix{Float64}
"""
function trajectory_reduce(points::Matrix{Float64}, remain::Int64)::Matrix{Float64}
    points[trunc.(Int, range(1, stop=size(points, 1), length=remain + 1)[1:end-1]), :]
end


"""
    trajectory_sort(points::Vector{Vector{Float64}}; verify_sort::Bool=false)::Matrix{Float64}
"""
function trajectory_sort(points::Vector{Vector{Float64}}; verify_sort::Bool=false)::Matrix{Float64}

    _points = points

    sorted_points = []
    push!(sorted_points, popfirst!(_points))


    while length(_points) > 0
        min_dist = 100000
        point = nothing

        for p in _points
            dist = point_distance(p, sorted_points[end])

            if dist < min_dist
                min_dist = dist
                point = p
            end
        end

        push!(sorted_points, point)
        filter!(e -> e != point, _points)
    end

    # Verify the sorting
    # Sometimes, when multiple points with same distance are present,
    # the result is not correct. This checks for points outside the
    # expected sqrt(2)*_grid distance and tries to move them to better
    # positions
    spoints = mapreduce(permutedims, vcat, sorted_points) # convert 2d vector to matrix

    if verify_sort == true
        # Obtain grid size
        _grid = minimum(abs.(minimum(spoints[2:end, :] - spoints[1:end-1, :])) for u in [unique(sort(c[:])) for c in eachcol(points)])

        while true

            # Get distances between consecutive points
            _dists = points_distance(spoints)

            # Find outliers
            _outliers = [ _i for (_i, _d) in enumerate(_dists) if _d > sqrt(2)*_grid ]#_dists[_dists.>sqrt(2)*_grid]

            if length(_outliers) == 1
                println("trajectorySort: Only one large jump in the trajectory found.")
                @printf("trajectorySort: points = %s\n", spoints)
                @printf("trajectorySort: dists = %s\n", _dists)
                println("trajectorySort: Continuing without dealing with outliers.")
                break
            # Continue only if outliers found
            elseif length(_outliers) > 0
                # TODO: Tohle v C verzi ted nebude fungovat.
                # Outlier indices
                _oi = findall(_dists .> sqrt(2) * _grid)

                # Find group sizes
                _groups = #(
                    [(_oi[_i][1] + 1, _oi[_i+1][1], _oi[_i+1][1] - _oi[_i][1]) for _i in 1:(length(_oi)-1)]
                #) # start id, end id (both inclusive), size of the group


                _groups = [_groups..., (_oi[end][1] + 1, _oi[1][1], _oi[1][1] + length(spoints) - _oi[end][1])]

                # Sort the groups in order to find the largest group
                _groups = sort(_groups, by=x->(x[end]), rev=true)

                # Delete outlier
                spoints = spoints[1:end.!=_groups[1][2]%length(spoints), :]
            else
                break
            end
        end
    end

    return spoints
end


if (abspath(PROGRAM_FILE) == @__FILE__)
    a = [[0.16433, 0.524746], [0.730177, 0.787651], [0.646905, 0.0135035], [0.796598, 0.0387711], [0.442782, 0.753235], [0.832315, 0.483352], [0.442524, 0.912381], [0.336651, 0.236891], [0.0954936, 0.303086], [0.459189, 0.374318]]
    b = [0.7589091211161472,
        0.8091539348190575,
        0.5256478329286531,
        0.41357337873861466]

    println(trajectory_sort(a; verify_sort=true))
end
