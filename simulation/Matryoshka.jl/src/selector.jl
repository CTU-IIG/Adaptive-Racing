# using .ParameterListClass.ParameterClass
# using .ParameterListClass

P_sel = ParameterList()
add_parameter!(P_sel, Parameter("sampling_distance", 1.0, 1.0, float, "[m] Distance of super-sampling before the interpolation, skipped when 0.", "init"))
add_parameter!(P_sel, Parameter("distance", 0.05, 0.05, float, "[m] Distance between the individual points, ignored when 0, used when requesting negative number of points.", "init"))
#=
add_parameter!(P_sel, Parameter("rotate", 0, 0, float, "Parameter for rotating the input path. 0 is not rotated. <0, 1)", "init"))
add_parameter!(P_sel, Parameter("fixed_points", [], [], Array, "Points to be used in the selection upon calling 'select'.", "init"))
=#
add_parameter!(P_sel, Parameter("overlap", 100, 0, float, "Size of the trajectory overlap. 0 disables this.", "init"))


SELECTOR = "uniform_time"


"""
    selector_init(kwargs)

Select and initialize selector with given parameters.
"""
function selector_init(kwargs)
    global SELECTOR

    #TODO: rotate
    update_all!(P_sel, kwargs)

    _selector = get(kwargs, :selector, SELECTOR)

    if _selector in ["curvature", "uniform", "curvature_sample", "fixed", "curvature2"]
        write(kwargs[:logfile], "Passed unsupported selector: $(_selector). Using $(SELECTOR) instead.\n")
    elseif _selector in ["uniform_distance", "uniform_time"]
        SELECTOR = _selector
    else
        throw(ArgumentError("Unknown selector: " * String(_selector)))
    end
end


"""
    select(points::Matrix{Float64}, remain::Int; overflown...)::Matrix{Float64}

Select points from `points` according to the selected selector and `remain`.

When `remain` is negative number of points is selected by the selector.
(In case that it supports negative selection.)
"""
function select(points::Matrix{Float64}, remain::Int; overflown...)::Matrix{Float64}
    global SELECTOR

    if remain < 0 && getValue(P, "distance") <= 0
        # Raise an exception, as we cannot proceed without further information.
        throw(ArgumentError("Negative selection requires set 'distance' parameter for 'uniform_distance' selector."))
    end

    if SELECTOR == "uniform_distance"
        rpoints = trajectory_resample(points, remain)

        # !Force number of points
        if remain > 0 && size(rpoints, 1) != remain
            return trajectory_resample(points, remain - 1)
        end
    elseif SELECTOR == "uniform_time"
        rpoints = trajectory_time_resample(points, remain)
    else
        throw(ArgumentError("Unsupported selector: " * String(SELECTOR)))
    end

    return rpoints
end


"""
    path_point_distance_avg(points::Matrix{Float64})::Float64

Get average distance between consecutive points.
"""
function path_point_distance_avg(points::Matrix{Float64})::Float64
    path_length(points) / length(points)
end


"""
    factor_compute(points::Matrix{Float64}, resolution::Float64)::Float64

Compute a factor for modifying the number of line points in order to obtain `resolution`.
"""
function factor_compute(points::Matrix{Float64}, resolution::Float64)::Float64
    path_point_distance_avg(points) / resolution
end


"""
    resolution_estimate(points::Matrix{Float64}, resolution::Float64)::Int64

Estimate the number of points of a line to obtain `resolution`.
"""
function resolution_estimate(points::Matrix{Float64}, resolution::Float64)::Int64
    trunc(Int, length(points) * factor_compute(points, resolution))
end


"""
    trajectory_resample(points::Matrix{Float64}, remain::Int64)::Matrix{Float64}

Resample `points` to obtain sequence of points of length `remain`.

# Note
Resampling is affected by parameters:
- fixed_points (currently disabled)
- sampling_distance
- distance
"""
function trajectory_resample(points::Matrix{Float64}, remain::Int64)::Matrix{Float64}

    # Throw away repeated point
    if points[1, :1] == points[end, :1]
        points = points[1:end-1, :]
    end

    # Keep fixed points local
    raw_fixed_points = []#copy(get_value(P_sel, "fixed_points"))

    # Result
    rpoints = []
    # Intermediate results
    fixed_points = Vector{Float64}()
    upoints = []

    # Other values
    #rotate = typeof(get_value(P_sel, "rotate")) != Vector ? [get_value(P_sel, "rotate") for _ in range(1, stop=max(1, length(get_value(P_sel, "fixed_points"))))] : copy(get_value(P_sel, "rotate"))

    while true
        # Rotate to get to the first fixed point
        _points = circshift(points, length(raw_fixed_points) > 0 ? -trajectory_closest_index(points, popfirst!(raw_fixed_points)) : 0)

        # Resample if requested
        if get_value(P_sel, "sampling_distance") != 0.0
            _points = interpolate(_points[:, 1:2], int_size=resolution_estimate(_points, get_value(P_sel, "sampling_distance")))
        end

        # Select points equidistantly
        if remain < 0
            _rpoints = interpolate(_points[:, 1:2], int_size=resolution_estimate(_points, get_value(P_sel, "distance")))
        # Select 'remain' points
        else
            _rpoints = interpolate(_points[:, 1:2], int_size=remain)
        end

        # Rotate when required
        #=
        if rotate[1] > 0.0
            #TODO
        else
            popfirst!(rotate)
        =#

        # Create fpoints with a set factor to allow concatenating
        _fpoints = interpolate(_points[:, 1:2], int_size=10 * length(_rpoints))

        push!(fixed_points, _fpoints[1])
        push!(upoints, _fpoints)
        #=
        end
        =#

        push!(rpoints, _rpoints)

        if length(raw_fixed_points) <= 0
            break
        end
    end

    #=
    if length(rpoints) == 1
    =#
    return rpoints[1]
    #=
    # This is a part for fixed points. But it is not working properly.
    else
        # Build up the new path waypoints
        result = nothing

        for _i in range(1, stop=length(rpoints))

            # 1) Take the fixed point from the next segment
            _p = fixed_points[((_i+1)%length(rpoints))+1, :]

            # 2) Find it in current path (rotated, full)
            _cpi = trajectory_closest_index(upoints[_i], _p; from_left=true)

            # 3) Loop through the selection to find all points that are more to the left
            _max_i = 0

            while _max_i + 1 < length(rpoints[_i]) && trajectory_closest_index(upoints[_i], rpoints[_i][_max_i+1, :]', from_left=true) < _cpi
                _max_i += 1
            end

            # 4) Append them to the result
            if _max_i >= 1
                if result == Nothing
                    result = rpoints[_i][1:_max_i+1, :]
                else
                    result = vcat(result, rpoints[_i][1:_max_i+1, :])
                end
            end
        end
        _l = result .+ 100
        @printf "trajectory_resample: %s %s out: %s\n" typeof(points) typeof(remain) typeof(_l)
        return _l
    end
    =#
end


"""
    trajectory_time_resample(points::Matrix{Float64}, remain::Int64)::Matrix{Float64}

Resample the given points to obtain ordered sequence of `remain` points that are isochronal.

# Note
Negative selection is not supported.
Resampling is affected by parameters:
- fixed_points (currently disabled)
- sampling_distance
- distance
- overlap
"""
function trajectory_time_resample(points::Matrix{Float64}, remain::Int64)::Matrix{Float64}
    if remain < 0
        throw(
            ArgumentError("Negative selection is not supported by 'uniform_time' selector.")
        )
    end

    if get_value(P_sel, "distance") <= 0
        throw(
            ArgumentError("Selector 'uniform_time' requires 'distance' parameter to be set '>0.0'.")
        )
    end

    if get_value(P_sel, "overlap") <= 0
        println(stderr, "Warning: Consider setting 'overlap' parameter, as otherwise, the initial conditions affect the results.")
    end


    resampled_trajectory = trajectory_resample(vcat(points), -1)

    _, _, _t = jazar_profile_compute(resampled_trajectory, get_value(P_sel, "overlap"))

    # TODO: Rotate is not supported as it needs to be done in 'trajectory_resample' first.
    #function trajectoryRotate(points, npi, rotation = 0.0)
	#    # Tohle stejně není v té jl implementaci podporované.
	#    return circshift(points, -floor(npi * rotation))
    #end
    #rotated_points = trajectoryRotate(resample, next_point_index, 0.0)
    rotated_points = resampled_trajectory

    # Time sample
    equidistant_trajectory = hcat(
        [ 
            rotated_points[argmin(broadcast(abs, _t .- i))[1], :] 
            for i in range(0, _t[end], remain+1)[1:end-1]
        ]...
    )'

    return vcat(
        hcat(
            [
                points[trajectory_closest_index(hcat(points), resampled_point), :]
                for resampled_point in collect.(eachrow(equidistant_trajectory))
            ]...
        )'
    )
end


if (abspath(PROGRAM_FILE) == @__FILE__)
    using NPZ
    START_POINTS = npzread("configuration/ng_start_points_torino2.npy")
    res = trajectory_resample(vcat(START_POINTS), 12)
    println(res)
    println(size(res, 1))
end
