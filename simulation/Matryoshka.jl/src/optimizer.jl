#import Statistics: mean
#import Dierckx: Spline2D
#import Printf: @printf
import Evolutionary
#using Metaheuristics
#import Gnuplot: @gp
#using PyCall
#import Colors: @colorant_str

#f_helper(x) = x
#f_helper(x::Vector) = [f_helper(e) for e in x]
#f_helper(d::Dict) = Dict(Symbol(k) => f_helper(v) for (k, v) in d)
#symbol_dict(d::Dict) = f_helper(d)

# Global variables
OPTIMIZER = nothing
MATRYOSHKA = nothing
VALID_POINTS = nothing
CRITERION_ARGS = nothing
INTERPOLATOR_ARGS = nothing
SEGMENTATOR_ARGS = nothing
SELECTOR_ARGS = nothing
PENALIZER_INIT = nothing
PENALIZER_ARGS = nothing
LOGFILE = nothing
VERBOSITY = 3
FILELOCK = ReentrantLock()
HOLDMAP = nothing
GRID = nothing
PENALTY = nothing
FIGURE = nothing
PLOT = nothing
BUDGET = nothing
NUM_WORKERS = Sys.CPU_THREADS # Currently not used for Evolutionary.jl

mutable struct TrackBitmap
    const p_min::Vector{Float64} # const here needs Julia 1.8
    const p_max::Vector{Float64}
    const p_step::Vector{Float64}
    bitmap::BitArray{2}
    TrackBitmap(p_min, p_max, p_step) = new(p_min, p_max, p_step,
                                            falses((cld.(p_max .- p_min, p_step) .|> Int) + [1, 1] |> Tuple))
end

function index(tb::TrackBitmap, point)
    idx = ((point .- tb.p_min) .÷ tb.p_step .|> Int) + [1, 1]
end

Base.setindex!(tb::TrackBitmap, x, pt) = tb.bitmap[index(tb, pt)...] = x
function Base.getindex(tb::TrackBitmap, pt)
    idx = index(tb, pt)
    if any(idx .< (1,1)) || any(idx .> size(tb.bitmap))
        return false
    end
    return tb.bitmap[idx...]
end

TRACK_BITMAP = nothing

P = ParameterList()

add_parameter!(P, Parameter("budget", 100, 100, Int, "Budget parameter for the genetic algorithm.", "init (general)"))
add_parameter!(P, Parameter("groups", 8, 8, Int, "Number of groups to segmentate the track into.", "init (general)"))
add_parameter!(P, Parameter("workers", "Sys.CPU_CORES", "Sys.CPU_CORES", Int, "Number threads for the genetic algorithm.", "init (general)"))
add_parameter!(P, Parameter("penalty", 100, 100, Float64, "Constant used for increasing the penalty criterion.", "init (general)"))
add_parameter!(P, Parameter("criterion_args", Dict(), Dict(), Dict, "Arguments for the criterion function.", "init (general)"))
add_parameter!(P, Parameter("interpolator_args", Dict(), Dict(), Dict, "Arguments for the interpolator function.", "init (general)"))
add_parameter!(P, Parameter("segmentator_args", Dict(), Dict(), Dict, "Arguments for the segmentator function.", "init (general)"))
add_parameter!(P, Parameter("selector_args", Dict(), Dict(), Dict, "Arguments for the selector function.", "init (general)"))
add_parameter!(P, Parameter("penalizer_init", Dict(), Dict(), Dict, "Arguments for the init part of the penalizer function.", "init (general)"))
add_parameter!(P, Parameter("penalizer_args", Dict(), Dict(), Dict, "Arguments for the penalizer function.", "init (general)"))
add_parameter!(P, Parameter("logging_verbosity", 2, 2, Int, "Index for verbosity of the logger.", "init (general)"))
add_parameter!(P, Parameter("hold_matryoshka", false, false, Bool, "Whether the transformation should be created only once.", "init (Matryoshka)"))
add_parameter!(P, Parameter("plot", false, false, Bool, "Whether a graphical representation should be created.", "init (viz.)"))
add_parameter!(P, Parameter("grid", "computed by default", "computed by default", Vector, "X-size and y-size of the grid used for points discretization.", "init (Matryoshka)"))
add_parameter!(P, Parameter("plot_mapping", false, false, Bool, "Whether a grid should be mapped onto the track (to show the mapping).", "init (viz.)"))

######################
# Utils
######################

"""
    optimizer_init(kwargs)

Initialize optimizer for Matryoshka mapping.
"""
function optimizer_init(kwargs)"""points,
    group_centers,
    group_centerline,
    budget::Int=10,
    layers::Int=5,
    groups::Int=12,
    workers::Int=Sys.CPU_THREADS,
    penalty=100,
    criterion_args::Dict=Dict(),
    interpolator_args::Dict=Dict(),
    segmentator_args::Dict=Dict(),
    selector_args::Dict=Dict(),
    penalizer_init::Dict=Dict(),
    penalizer_args::Dict=Dict(),
    logfile=stdout,
    logging_verbosity::Int=2,
    hold_matryoshka::Bool=false,
    plot::Bool=false,
    grid::Vector=[],
    figure=nothing,
    #penalize = nothing,
    kwargs)"""

    points = get(kwargs, :points, [])
    group_centers = get(kwargs, :group_centers, [])
    group_centerline = get(kwargs, :group_centerline, [])
    budget = get(kwargs, :budget, 10)
    layers = get(kwargs, :layers, 5)
    groups = get(kwargs, :groups, 12)
    workers = 1#Sys.CPU_THREADS,
    penalty = get(kwargs, :penalty, 100)
    criterion_args = get(kwargs, :criterion_args, Dict())
    interpolator_args = get(kwargs, :interpolator_args, Dict())
    segmentator_args = get(kwargs, :segmentator_args, Dict())
    selector_args = get(kwargs, :selector_args, Dict())
    penalizer_init = get(kwargs, :penalizer_init, Dict())
    penalizer_args = get(kwargs, :penalizer_args, Dict())
    logfile = get(kwargs, :logfile, stdout)
    logging_verbosity = get(kwargs, :logging_verbosity, 2)
    hold_matryoshka = get(kwargs, :hold_matryoshka, false)
    plot = get(kwargs, :plot, false)
    grid = get(kwargs, :grid, [])
    figure = get(kwargs, :figure, nothing)

    global OPTIMIZER, MATRYOSHKA, VALID_POINTS, LOGFILE, VERBOSITY, HOLDMAP, GRID, PENALTY, FIGURE, PLOT, BUDGET, NUM_WORKERS
    global CRITERION, CRITERION_ARGS, INTERPOLATOR, INTERPOLATOR_ARGS, SEGMENTATOR, SEGMENTATOR_ARGS, SELECTOR, SELECTOR_ARGS, PENALIZER, PENALIZER_INIT, PENALIZER_ARGS

    # Local to global variables
    CRITERION_ARGS = symbol_dict(criterion_args)
    INTERPOLATOR_ARGS = symbol_dict(interpolator_args)
    SEGMENTATOR_ARGS = symbol_dict(segmentator_args)
    SELECTOR_ARGS = symbol_dict(selector_args)
    PENALIZER_INIT = symbol_dict(penalizer_init)
    PENALIZER_ARGS = symbol_dict(penalizer_args)
    LOGFILE = logfile#stdout#logfile
    VERBOSITY = logging_verbosity
    _holdmatryoshka = hold_matryoshka
    PENALTY = penalty
    FIGURE = figure
    PLOT = plot
    BUDGET = budget
    NUM_WORKERS = workers

    if haskey(kwargs, :algorithm)
        if kwargs[:algorithm] != "matryoshka"
            write(logfile, "Passed unsupported optimizer: $(kwargs[:algorithm]). Using matryoshka instead.\n")
        end
    end

    VALID_POINTS = points
    if MATRYOSHKA === nothing || _holdmatryoshka == false
        # Note: In version <=1.3.0 the group_centerline passed to the SELECTOR was sorted using
        #       ng_trajectory.interpolators.utils.trajectorySort, but it sometimes rotated the
        #       already sorted centerline; interestingly, the result was counterclockwise at all
        #       times (or at least very very often).

        group_centers = select(group_centerline, groups; SELECTOR_ARGS...)

        if plot == true
            # TODO: plot
        end

        # Matryoshka construction
        _groups = segmentate(points, group_centers; SEGMENTATOR_ARGS...)

        grouplayers = groups_border_obtain(_groups)
        grouplayers = groups_border_beautify(grouplayers, 400)

        if plot == true
            # TODO: plot
        end

        layers_center = groups_center_compute(_groups)
        layers_count = [layers for _ in 1:length(grouplayers)]

        MATRYOSHKA = [matryoshka_create(grouplayers[_i], layers_center[_i], layers_count[_i]) for _i in 1:length(_groups)]
        # TODO: plot

        println("Matryoshka mapping constructed.")

        if GRID === nothing
            GRID = grid_compute(points)
        end
    end
end

#=
using VideoIO
import FileIO

WRITER = nothing

function plot_population_gp(population, value; video=false)
    n = length(MATRYOSHKA)
    @gp tit="Best value: $(@sprintf "%5.3f" value)" "set size ratio -1" :-
    #@gp :- VALID_POINTS[:, 1] VALID_POINTS[:, 2] "w p pt 1 lc rgbcolor '0xeeeeee' notitle" :-
    #let tb = TRACK_BITMAP
    #    img = map(tb.bitmap') do p p ? colorant"#eed" : colorant"white" end
    #let tb = VALID_POINTS
    #    img = foreach(tb) do p p ? colorant"#eed" : colorant"white" end
    #    @gp :- Gnuplot.palette(:gray1) "set size ratio -1" :-
    #    @gp :- recipe(img, "dx=$(tb.p_step[1]) dy=$(tb.p_step[2]) origin=($(tb.p_min[1]), $(tb.p_min[2]))") :-
    #end

    foreach(population) do p
        points01 = reshape(p, (n, 2))
        points = zeros(n, 2)
        for (i, p) in enumerate(eachrow(points01))
            points[i, :] .= matryoshka_map(MATRYOSHKA[i], [p])[1]
        end
        _points = interpolate(points)
        @gp :- _points[:, 1] _points[:, 2] "w l notitle" :-
        @gp :- points[:, 1] points[:, 2] "w p pt 1 lc -1 notitle" :-
    end
    @gp :- ""
    save(term="pngcairo size 1280, 720 fontscale 0.8", output="frame.png")
    #write(WRITER, FileIO.load("frame.png"))
end

function plot_population(population, value; video=false)
    n = length(MATRYOSHKA)
    plot(plot_title="Best value: $(@sprintf "%5.3f" value)")
    video && plot!(size=(1280, 720))
    let tb = TRACK_BITMAP
        img = map(tb.bitmap') do p
            p ? colorant"#eed" : colorant"white"
        end
        plot!(xlims=[tb.p_min[1], tb.p_max[1]], ylims=[tb.p_min[2], tb.p_max[2]])
        plot!([tb.p_min[1], tb.p_max[1]], [tb.p_min[2], tb.p_max[2]],
              img, yflip=false)
    end

    foreach(population) do p
        points01 = reshape(p, (n, 2))
        points = zeros(n, 2)
        for (i, p) in enumerate(eachrow(points01))
            points[i, :] .= matryoshka_map(MATRYOSHKA[i], [p])[1]
        end
        _points = interpolate(points)
        plot!(_points[:, 1], _points[:, 2], legend=false)
        plot!(points[:, 1], points[:, 2], seriestype=:scatter, shape=:+, mc=:black)
    end

    if video
        Plots.png("frame.png")
        #write(WRITER, FileIO.load("frame.png"))
    else
        plot!(show=true)
    end
end
=#

#=
function Evolutionary.trace!(record::Dict{String,Any}, objfun, state, population, method, options)
    #plot_population(population, value(state), video=false)
end

function Evolutionary.after_while!(objfun, state, method, options)
    global STATE = state # allow investigating the internal state
                         # after the optimization ends
end
=#


"""
    optimize_evolutionary()

See also [`optimize`](@ref) and [`opt`](@ref).
"""
@inline function optimize_evolutionary()
    global OPTIMIZER, MATRYOSHKA, LOGFILE, FILELOCK, VERBOSITY, INTERPOLATOR, INTERPOLATOR_ARGS, FIGURE, PLOT, PENALIZER, PENALIZER_ARGS
    n = length(MATRYOSHKA)
    constr = Evolutionary.BoxConstraints(zeros(2n), ones(2n))
    x0 = fill(0.5, 2n)
    #method = Evolutionary.GA(populationSize=30, selection=uniformranking(10), mutation=gaussian(0.1), crossover=TPX)
    method = Evolutionary.CMAES(sigma0=0.1, c_1=0.01, c_mu=0.001, c_sigma=0.02)
    #method = Evolutionary.CMAES(sigma0=0.1)
    res = Evolutionary.optimize(opt, constr, x0, method,
                                Evolutionary.Options(iterations=BUDGET,
                                                     #parallelization=:thread,
                                                     show_trace=false,
                                                     store_trace=true,
                                                     #reltol=1e-6,              # Comment out for notolerance run
                                                     successive_f_tol=10000,    # Uncomment for notolerance run
                                                     ))
    write(LOGFILE, string(res))
    write(LOGFILE, "\n")
    global RESULT = res
    trace = Evolutionary.trace(res)
    #@gp :trace value.(trace)[10:end] "w l" "set grid"

    return points01 = reshape(Evolutionary.minimizer(res), (n, 2))
end

#=
function optimize_metaheuristics()
    global OPTIMIZER, MATRYOSHKA, LOGFILE, FILELOCK, VERBOSITY, INTERPOLATOR, INTERPOLATOR_ARGS, FIGURE, PLOT, PENALIZER, PENALIZER_ARGS
    x0 = [0.5 for _ in 1:length(MATRYOSHKA)*2]
    bounds = repeat([0.0, 1.0], 1, length(MATRYOSHKA) * 2)
    ga = Metaheuristics.GA(;
                           crossover=Metaheuristics.OrderCrossover(),
                           mutation=Metaheuristics.SlightMutation())
    points01 = Metaheuristics.optimize(_opt, bounds, ga)
end

function optimize_nevergrad()
    global OPTIMIZER, MATRYOSHKA, LOGFILE, FILELOCK, VERBOSITY, INTERPOLATOR, INTERPOLATOR_ARGS, FIGURE, PLOT, PENALIZER, PENALIZER_ARGS

    num_rows = length(MATRYOSHKA)

    nevergrad = pyimport("nevergrad")
    # concurrent = pyimport("concurrent")
    # multiprocessing = pyimport("multiprocessing")

    # Optimizer definition
    instrum = nevergrad.Instrumentation(nevergrad.var.Array(num_rows, 2).bounded(0, 1))
    OPTIMIZER = nevergrad.optimizers.DoubleFastGADiscreteOnePlusOne(instrumentation=instrum, budget=BUDGET)
    # OPTIMIZER = nevergrad.optimizers.DoubleFastGADiscreteOnePlusOne(instrumentation=instrum, budget=BUDGET, num_workers=NUM_WORKERS)

    # @pywith concurrent.futures.ProcessPoolExecutor(max_workers=OPTIMIZER.num_workers, mp_context=multiprocessing.get_context("fork")) as executor begin
    #     recommendation = OPTIMIZER.minimize(_opt, executor=executor, batch_mode=false)
    # end
    recommendation = OPTIMIZER.minimize(_opt, batch_mode=false)
    points01 = convert(Array{Float64,2}, recommendation.args[1])
end
=#


"""
    optimize()

Run the optimization using Matryoshka mapping and Evolutionary.jl.
"""
function optimize()
    global OPTIMIZER, MATRYOSHKA, LOGFILE, FILELOCK, VERBOSITY, INTERPOLATOR, INTERPOLATOR_ARGS, FIGURE, PLOT, PENALIZER, PENALIZER_ARGS

    #encoder_options = (crf=23, preset="medium")

    #global WRITER
    #WRITER = open_video_out("video.mp4", RGB{N0f8}, (720, 1280),
    #                        framerate=10, encoder_options=encoder_options,
    #                        target_pix_fmt=VideoIO.AV_PIX_FMT_YUV420P)

    points01 = optimize_evolutionary()

    points = [matryoshka_map(MATRYOSHKA[i], [p])[1] for (i, p) in enumerate(eachrow(points01))]

    PENALIZER_ARGS[:optimization] = false
    final = _opt(points01)
    #plot_population_gp([points01], final)

    # Interpolate received points
    # It is expected that they are unique and sorted.
    _points = interpolate(mapreduce(permutedims, vcat, points))


    #close_video_out!(WRITER)

    lock(FILELOCK) do
        if VERBOSITY > 0
            write(LOGFILE, "solution:$(points)\n")
            write(LOGFILE, "final:$(final)\n")
        end
    end

    return (final, points, points01, _points)
end

prepare_points(points::Array{Float64, 2}, matr_len) = convert(Array{Float64,2}, points) # Nevergrad
prepare_points(points::Vector{Float64}, matr_len)::Matrix{Float64} = reshape(points, (matr_len, 2)) # Evolutionary


"""
    opt(points::Vector{Float64})::Float64

Function executed at every optimization step.

See also [`_opt`](@ref).

# TODO
Resolve variable type confict, as `points` can be Vector or Matrix.
"""
function opt(points)::Float64
    #@time _opt(points)
    _opt(points)
end


"""
    _opt(points)::Float64
"""
function _opt(points)::Float64
    global VALID_POINTS, CRITERION_ARGS, INTERPOLATOR_ARGS, PENALIZER_ARGS
    global MATRYOSHKA, LOGFILE, FILELOCK, VERBOSITY, GRID, PENALTY, penalizeZ

    points = prepare_points(points, length(MATRYOSHKA))

    # Transform points
    points = [matryoshka_map(MATRYOSHKA[i], [p])[1] for (i, p) in enumerate(eachrow(points))]
    _points = interpolate(mapreduce(permutedims, vcat, points); INTERPOLATOR_ARGS...)

    # Check the correctness of the points and compute penalty
    #penalty = penalizeZ(_points, VALID_POINTS, GRID, PENALTY; PENALIZER_ARGS...)
    # Předat py funkci jde, ale pak je to děsně pomalé.
    penalty = penalize(_points, VALID_POINTS, GRID, PENALTY; PENALIZER_ARGS...)

    if penalty != 0
        if VERBOSITY > 0
            lock(FILELOCK) do
                if VERBOSITY > 2
                    write(LOGFILE, "pointsA:$(string(points))\n")
                    write(LOGFILE, "pointsT:$(string(_points))\n")
                end
                if VERBOSITY > 1
                    write(LOGFILE, "penalty:$(penalty)\n")
                end
                #flush(LOGFILE)
            end
        end
        return Float64(penalty)
    end

    _c = compute(_points; CRITERION_ARGS...)
    lock(FILELOCK) do
        if VERBOSITY > 2
            write(LOGFILE, "pointsA:$(string(points))\n")
            write(LOGFILE, "pointsT:$(string(_points))\n")
        end
        if VERBOSITY > 1
            write(LOGFILE, "correct:$(_c)\n")
        end
        #flush(LOGFILE)
    end
    return _c
end


if (abspath(PROGRAM_FILE) == @__FILE__)
    # grouplayers = interpolate(a, 400)
    # x = [2 * i for i in 0:19]
    # y = [i for i in 0:19]
    # z = (x .- 1) .^ 2 .+ y .^ 2
    # spline = Spline2D(x, y, z, s=1e-4)
    # # println("result: ", spline(a[:, 1], a[:, 2]))
    # for _c in eachrow(a)
    #     println(spline(_c[1], _c[2]))
    # end
end
