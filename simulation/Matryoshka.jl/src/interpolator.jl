import Dierckx: Spline1D, derivative

INTERPOLATOR = "cubic_spline"


"""
    interpolator_init(kwargs)::Nothing

Select and initialize interpolator. Currently only checking that 'cubic_spline' is selected.
"""
function interpolator_init(kwargs)::Nothing
    global INTERPOLATOR

    if haskey(kwargs, :interpolator)
        if kwargs[:interpolator] == "none"
            write(kwargs[:logfile], "Interpolator 'none' does not currently return curvature.\n")
            INTERPOLATOR = "none"
            nothing  # For some jl reason this needs to be here, otherwise
            # MethodError: Cannot `convert` an object of type String to an object of type Nothing
        elseif kwargs[:interpolator] != "cubic_spline"
            write(kwargs[:logfile], "Passed unsupported interpolator: $(kwargs[:interpolator]). Defaulting to cubic_spline.\n")
        end
    end
end


"""
    interpolate(points::Matrix{Float64}; int_size::Int=400, overflown...)::Matrix{Float64}

Interpolate `points` with a selected interpolator to get `int_size` points.
"""
function interpolate(points::Matrix{Float64}; int_size::Int=400, overflown...)::Matrix{Float64}
    global INTERPOLATOR

    if INTERPOLATOR == "none"
        return hcat(points, zeros(size(points, 1), 1))
    end

    _points = vcat(points, points[1:1, :])
    x, y = eachcol(_points)
    distance = vec(cumsum(sqrt.(sum(diff(_points, dims=1) .^ 2, dims=2)), dims=1))
    distance = insert!(distance, 1, 0) ./ last(distance)

    alpha = range(0, stop=1, length=int_size + 1)[1:end-1]

    splinex = Spline1D(distance, x; periodic=true)
    spliney = Spline1D(distance, y; periodic=true)

    ipolx = splinex(alpha)
    ipol2x = derivative(splinex, alpha; nu=2)
    ipol1x = derivative(splinex, alpha; nu=1)

    ipoly = spliney(alpha)
    ipol2y = derivative(spliney, alpha; nu=2)
    ipol1y = derivative(spliney, alpha; nu=1)

    xx = (ipol1x .* ipol2y) .- (ipol1y .* ipol2x)
    yy = sqrt.((ipol1x .^ 2 .+ ipol1y .^ 2) .^ 3)

    return hcat(ipolx, ipoly, xx ./ yy)
end

if (abspath(PROGRAM_FILE) == @__FILE__)
    using Plots
    a = [0.16433 0.524746
        0.730177 0.787651
        0.646905 0.0135035
        0.796598 0.0387711
        0.442782 0.753235
        0.832315 0.483352
        0.442524 0.912381
        0.336651 0.236891
        0.0954936 0.303086
        0.459189 0.374318]
    b = interpolate(a)
    plot(b[:, 1], b[:, 2])
    savefig("myplot.png")
end
