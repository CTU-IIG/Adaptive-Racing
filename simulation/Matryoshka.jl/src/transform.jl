import Dierckx: Spline2D
import Statistics: mean


"""
    matryoshka_create(layer0::Matrix{Float64}, layer0_center::Matrix{Float64}, layer_count::Int)::Vector{Dierckx.Spline2D}

Create interpolators for Matryoshka transformation on given Segment.

The Segment is described using its border `layer0` and center `layer0_center`. To improve the Matryoshka mapping,
`layer_count` layers are created on the Segment (by scaling the border).

Function returns two Spline2D objects, that converts transformed coordinates to a real coordinate.
"""
function matryoshka_create(layer0::Matrix{Float64}, layer0_center::Matrix{Float64}, layer_count::Int)::Vector{Spline2D}

    layers = group_layers_compute(layer0, layer0_center, layer_count)
    layer0_size = Int(size(layer0, 1))

    # Method, where we learn the interpolator all layers at once and use that information
    # for computing of the transformation
    # Note: This is much more precise.

    _tc = Array{Float64}(undef, 0, 2)
    _rc = Array{Float64}(undef, 0, 2)

    for _layer_index in 0:layer_count-1
        layer_params = layer_index_to_parameters(_layer_index, layer0_size, layer_count)
        indices = 1:layer_params[1]
        tc = indices_to_transformed_coordinates(indices, layer_params...)
        rc = indices_to_real_coordinates(indices, layers[Int(_layer_index)+1])

        _tc = vcat(_tc, tc)#mapreduce(permutedims, vcat, tc))
        _rc = vcat(_rc, rc)#mapreduce(permutedims, vcat, rc))
    end

    # transformedCoordinatesToRealCoordinates
    _ip2d = Vector{Spline2D}()

    for _d in axes(_rc, 2)
        push!(_ip2d, Spline2D(_tc[:, 1], _tc[:, 2], _rc[:, _d]; s=length(_tc[:, 1]) - sqrt(2 * length(_tc[:, 1]))))
    end

    return _ip2d
end


"""
    matryoshka_map(matryoshka::Vector{Spline2D}, coords)::Vector{Vector{Float64}}

Transform points using Matryoshka mapping.

See also [`matryoshka_create`](@ref).
"""
function matryoshka_map(matryoshka::Vector{Spline2D}, coords)::Vector{Vector{Float64}}
    _rcoords = Vector{Vector{Float64}}()
    for _c in coords
        _dims = Vector{Float64}()
        for _interpolator in matryoshka
            push!(_dims, _interpolator(_c[1], _c[2]))
        end
        push!(_rcoords, _dims)
    end
    return _rcoords
end


"""
    layer_index_to_parameters(layer_index::Int, layer0_size::Int, layer_count::Int)::Tuple{Int64, Float64}

Compute layer parameters.

Computes size of a layer and its scale using `layer_index`. The layer size is scaled down accordingly to its index
and the number of layers.

Examples:
1 layer :
    layer 0 : 100% points, scale 1

2 layers:
    layer 0 : 100% points, scale 1
    layer 1 : 50% points, scale 0.5

3 layers:
    layer 0 : 100% points, scale 1
    layer 1 : 66% points, scale 0.66
    layer 2 : 33% points, scale 0.33
"""
function layer_index_to_parameters(layer_index::Int, layer0_size::Int, layer_count::Int)::Tuple{Int64, Float64}
    (trunc(Int, layer0_size - ((layer0_size / layer_count) * layer_index)), 1 - (1.0 / layer_count) * layer_index)
end


"""
    indices_to_transformed_coordinates(indices::UnitRange{Int64}, layer_size::Int, scale::Float64)::Matrix{Float64}

Convert `indices` of a layer to transformed coordinates.

Computed transformed coordinates are uniformly placed alongside the layer. Scale is provided as a size of a square size
of this coordinate frame.
"""
function indices_to_transformed_coordinates(indices::UnitRange{Int64}, layer_size::Int, scale::Float64)::Matrix{Float64}
    _coords = Vector{Vector{Float64}}()
    l = Float64(layer_size)
    for _i in indices
        _point = Vector{Float64}()
        for _d in [1, -1]
            push!(_point, (
                max(
                    min(
                        (abs(-((((l / 2.0) + _i + _d * (l / 8.0)) % l) - (l / 2.0))) / (l / 4.0)) - 0.5, 1),
                    0
                ) * scale + ((1 - scale) / 2.0)
            ))
        end
        push!(_coords, _point)
    end
    return reduce(vcat, _coords')
end


"""
Convert `indices` of a layer to real coordinates.

Given `points` are from the same layer and ordered.
"""
function indices_to_real_coordinates end


"""
    indices_to_real_coordinates(indices::UnitRange{Int64}, points::Matrix{Float64})::Matrix{Float64}
"""
@inline function indices_to_real_coordinates(indices::UnitRange{Int64}, points::Matrix{Float64})::Matrix{Float64}
    return points[indices, :]
end


"""
    indices_to_real_coordinates(indices::Vector{Float64}, points::Matrix{Float64})::Matrix{Float64}
"""
function indices_to_real_coordinates(indices::Vector{Float64}, points::Matrix{Float64})::Matrix{Float64}
    _rcoords = Vector{Vector{Float64}}()
    for _i in indices
        if trunc(Int, _i) == _i
            push!(_rcoords, points[trunc(Int, _i), :])
        else
            push!(_rcoords,
                points[trunc(Int, _i), :] .+ (_i - trunc(Int, _i)) .* (points[(trunc(Int, _i)+1)%size(points)[1], :] .- points[trunc(Int, _i), :])
            )
        end
    end
    return reduce(vcat, _rcoords')
end


"""
    points_filter(points::Matrix{Float64}, grid=nothing)::Vector{Vector{Float64}}

Filter out points that are not necessary in the trajectory planning.

# Notes
It is expected that points are aligned to a square grid.

# TODO
Do not use Vector of Vectors.
"""
function points_filter(points::Matrix{Float64}, grid=nothing)::Vector{Vector{Float64}}
    # _points = Array{Float64}(undef, 0, 2)
    _points = Vector{Float64}[]
    _grid = grid !== nothing ? grid : minimum(minimum(u[2:length(u)] - u[1:length(u)-1] for u in [unique(sort(c[:])) for c in eachcol(points)]))
    _cells = [[trunc(Int, (round.(_p[_d] / _grid))) for _d = 1:ndims(points)] for _p in eachrow(points)]
    _cells_copy = [[trunc(Int, (round.(_p[_d] / _grid))) for _d = 1:ndims(points)] for _p in eachrow(points)]
    points_l = points

    for _p in collect.(eachrow(points_l))
        _cell = [trunc(Int, (round.(_p[_d] / _grid))) for _d = 1:ndims(points)]

        _xr = -1:1
        _yr = -1:1

        x = sum([any([[_cell[1] + _x, _cell[2] + _y] ∈ _cells for _x in _xr]) for _y in _yr])
        y = sum([any([[_cell[1] + _x, _cell[2] + _y] ∈ _cells for _y in _yr]) for _x in _xr])

        if (x < 3) && (y < 3)
            # Return nearby points back to the loop
            for _xr in -1:1
                for _yr in -1:1
                    if _xr == _yr == 0
                        continue
                    end

                    _nearbyc = [_cell[1] + _xr, _cell[2] + _yr]

                    if _nearbyc in _cells_copy && points_l[findfirst(item -> item == _nearbyc, _cells_copy)] in _points
                        _nearbyp = points_l[findfirst(item -> item == _nearbyc, _cells_copy)]
                        filter!(e -> e ≠ _nearbyp, _points)
                        push!(points_l, _nearbyp)
                    end
                end
            end
            filter!(e -> e ≠ _cell, _cells)
        else
            push!(_points, _p)
        end
    end
    return _points
end


"""
    groups_border_obtain(groups::Vector{Matrix{Float64}}, grid::Float64=nothing)::Vector{Matrix{Float64}}

Obtain border points of the groups.

# Notes
It is expected that the points are aligned to a square grid.
Border points are not ordered.
"""
function groups_border_obtain(groups::Vector{Matrix{Float64}}, grid=nothing)::Vector{Matrix{Float64}}

    _borders = []

    for (_i, _g) in enumerate(groups)
        _border = Array{Float64}(undef, 0, 2)

        # Obtain grid size if not set
        _grid = grid !== nothing ? grid : minimum(minimum(u[2:end] - u[1:end-1] for u in [unique(sort(_g[:, d])) for d in 1:size(_g)[2]]))

        # Go through all dimensions
        for _d in 1:ndims(_g)

            # Find unique values in that dimension
            for _u in unique(sort(_g[:, _d]))
                _v = _g[findall(_g[:, _d] .== _u), :]

                # Append points with max / min values in another dimension
                _border = vcat(_border, minimum(_v, dims=1))
                _border = vcat(_border, maximum(_v, dims=1))

                # Append inner borders
                # Obtain values in the dimension
                _v = _g[findall(_g[:, _d] .== _u), :]

                # Sort them
                # Enforce the axis as otherwise it is not sorted in ascending order everytime.
                _v = _v[sortperm(_v[:, 3 - _d]), :]

                # Find distances between neighbouring points
				_dists = map(maximum, eachrow(_v[2:end, :] .- _v[1:end-1, :]))
				# Find points in the distance larger than 1.5x _grid
                _bords = findall(_dists .> (_grid * 1.5))

				for _b in _bords
                    _border = vcat(_border, _v[[_b], :])
                    _border = vcat(_border, _v[[_b+1], :])
                end
            end
        end
        push!(_borders, _border)
    end
    [unique(sortslices(b, dims=1, by=x -> (x[1], x[2])), dims=1) for b in _borders]
end


"""
    groups_border_beautify(borders::Vector{Matrix{Float64}}, border_length::Int)::Vector{Matrix{Float64}}

Filter, sort and interpolate the borders to get smoother points.
"""
function groups_border_beautify(borders::Vector{Matrix{Float64}}, border_length::Int)::Vector{Matrix{Float64}}
    bborders = Vector{Matrix{Float64}}()

    for (group_i, border) in enumerate(borders)
        # FIXME: Temporarily hidden as we are working with 0.02 map in Stage.
        border_filtered = points_filter(border)#0.05

        border_sorted = trajectory_sort(border_filtered, verify_sort=true)

        border_interpolated = interpolate(border_sorted, int_size=border_length)

        push!(bborders, border_interpolated)
    end

    return bborders
end


"""
    groups_center_compute(_groups::Vector{Matrix{Float64}})::Vector{Matrix{Float64}}
"""
@inline function groups_center_compute(_groups::Vector{Matrix{Float64}})::Vector{Matrix{Float64}}
    return mean.(_groups, dims=1)
end


"""
    group_layers_compute(layer0::Matrix{Float64}, layer0_center::Matrix{Float64}, layer_count::Int)::Vector{Matrix{Float64}}

Compute `layer_count` layers of a group in real coordinates.

When computing the layers of a segment, points are moved towards the `layer0_center`. Number of points in each layer is lowered.

See also [`layer_index_to_parameters`](@ref).
"""
function group_layers_compute(layer0::Matrix{Float64}, layer0_center::Matrix{Float64}, layer_count::Int)::Vector{Matrix{Float64}}
    layers_size = size(layer0, 1)
    points = [(layer0[:, 1:2] .- layer0_center) .* (1 - (1 / layer_count) * layer_index) .+ layer0_center for layer_index in 0:layer_count-1]
    remains = (trunc.(Int, layers_size - (layers_size / layer_count) * layer_index for layer_index in 0:layer_count-1))

    [trajectory_reduce(points[i], remains[i]) for i in eachindex(remains)]
end
