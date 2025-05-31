P_pen = ParameterList()
# TODO: Consider using only curvature (without count subroutine).
add_parameter!(P_pen, Parameter("k_max", 1.5, 1.5, float, "(curvature) [m^-1] Maximum allowed curvature in abs.", ""))

PENALIZER = "curvature"


"""
    penalizer_init(kwargs)

Select and initialize penalizer with given parameters.
"""
function penalizer_init(kwargs)
    global PENALIZER

    update_all!(P_pen, kwargs)

    _penalizer = get(kwargs, :penalizer, PENALIZER)

    if _penalizer in ["segment", "centerline", "borderlines"]
        write(kwargs[:logfile], "Passed unsupported selector: $(_penalizer). Using $(PENALIZER) instead.\n")
    elseif _penalizer in ["count", "curvature", "none"]
        PENALIZER = _penalizer
    else
        throw(ArgumentError("Unknown penalizer: " * String(_penalizer)))
    end
end


"""
    penalize(points::Matrix{Float64}, valid_points::Array{Float64, 2}, grid::Float64, penalty::Real = 100.0; overflown...)::Float64

Count penalty for candidate solution `points` using their position in the map `valid_points`.
"""
function penalize(points::Matrix{Float64}, valid_points::Array{Float64, 2}, grid::Float64, penalty::Real = 100; overflown...)::Float64
    global MAP, MAP_ORIGIN, PENALIZER

    update_all!(P_pen, overflown)

    invalid = 0

    try
        invalid = sum(map(eachrow(points)) do p
                          #index = Int.(round.((p[1:2] - MAP_ORIGIN) ./ grid) .+ 1)
                          #println(index, p)
                          #MAP[index[1], index[2]] == 100 ? 0 : 1 #TRACK_BITMAP[p[1:2]] ? 0 : 1
                          MAP[round(Int, ((p[1] - MAP_ORIGIN[1]) / grid) + 1),
                              round(Int, ((p[2] - MAP_ORIGIN[2]) / grid) + 1)] == 100 ? 0 : 1
                      end)
    catch e
        return 100 * penalty
    end

    if PENALIZER == "count"
        return invalid * penalty

    elseif PENALIZER == "curvature"

        if invalid == 0
            _k = get_value(P_pen, "k_max")

            invalid_1 = sum(map(eachrow(points)) do p
                          p[3] > _k ? p[3] : 0
                      end) / 100

            invalid_2 = sum(map(eachrow(points)) do p
                          p[3] < -_k ? p[3] : 0
                      end) / 100

            invalid = invalid_1 - invalid_2
        end

        return invalid * penalty * 10

    elseif PENALIZER == "none"
        return 0.0

    else
        throw(ArgumentError("Unsupported penalizer: " * String(PENALIZER)))
    end
end

if (abspath(PROGRAM_FILE) == @__FILE__)
    a = [ 0.16433    0.524746;
        0.730177   0.787651;
        0.646905   0.0135035;
        0.796598   0.0387711;
        0.442782   0.753235;
        0.832315   0.483352;
        0.442524   0.912381;
        0.336651   0.236891;
        0.0954936  0.303086;
        0.459189   0.374318]
    b = a .+ 0.1
    println(penalize(a, b, 0.03999999910593033))
end
