using JSON # print, parse, parsefile, json
import NPZ: npzread
using Printf # @printf, @sprintf, Format, format
#using Plots
using Gnuplot # @gp, palette # For some reason cannot be imported this way.
import PkgVersion: Version
import Colors: @colorant_str

f_helper(x) = x
f_helper(x::Vector) = [f_helper(e) for e in x]
f_helper(d::Dict) = Dict(Symbol(k) => f_helper(v) for (k, v) in d)
symbol_dict(d::Dict) = f_helper(d)

CONFIGURATION = nothing


"Loads a configuration stored in 'filename'."
function configuration_load(filename)

    global CONFIGURATION

    conf = nothing

    try
        conf = symbol_dict(JSON.parsefile(filename))

        if get(conf, :_version, 1) < 2
            println(stderr, "Unsupported version of the configuration file.")
        else
            CONFIGURATION = conf
        end
    catch e
        println(e)
        return false
    end

    if get(conf, :logging_verbosity, 1) > 1
        JSON.print(CONFIGURATION)
    end

    return true
end


function configuration_set(conf)
    global CONFIGURATION = conf
    return true
end


function cascade_run(track, fileformat, notification, loop_i, loop_output, conf)

    # Cascade step timing
    step_time = @elapsed begin

    # Get configuration for current step
    _alg = merge(conf, loop_i[2])

    # Rename output from previous stages
    fitness, rcandidate, tcandidate, result = loop_output

    # Check for logging file
    # If exists, we continue to next round.
    if fileformat !== ""
        try
            _stored = Dict{String, Any}()

            for line in eachline(replace(fileformat, "%d" => loop_i[1]) * "-$(get(_alg, :algorithm, "")).log")
                if startswith(line, "#")
                    name, value = split(line, ":")
                    merge!(_stored, Dict(name => JSON.parse(value)))
                end
            end

            if _stored["#fitness"] < fitness
                return (_stored["#fitness"],
                    [ Float64.(v) for v in _stored["#tcandidate"] ],
                    reduce(hcat, [ Float64.(v) for v in _stored["#rcandidate"] ])',
                    reduce(hcat, [ Float64.(v) for v in _stored["#trajectory"] ])'
                )
            else
                return loop_output
            end
        catch e
        end
    end

    # Open up logging file
    if fileformat !== ""
        LOGFILE = open(replace(fileformat, "%d" => loop_i[1]) * "-$(get(_alg, :algorithm, "")).log", "w")
        #write(LOGFILE, JSON.json(_alg))
        write(LOGFILE, "\nRunning Matryoshka.jl version " * string(Version(Matryoshka)) * "\n")
        flush(LOGFILE)
    else
        LOGFILE = stdout
    end


    ## Initialization
    # Get optimizers etc.
    # Note: This is not done here and is moved to `*_init`. Mostly because Julia has all algorithms in one file.

    # Show up current progress
    write(LOGFILE, replace(notification, "%d" => loop_i[1]) * " $(get(_alg, :algorithm, "")) with $(get(_alg, :criterion, "")) criterion (penalized by $(get(_alg, :penalizer, ""))), int. by $(get(_alg, :interpolator, ""))\n")
    flush(LOGFILE)

    # Initialize parts
    selector_init(merge(_alg, get(_alg, :selector_init, Dict()), Dict(:logfile => LOGFILE)))
    interpolator_init(merge(_alg, get(_alg, :interpolator_init, Dict()), Dict(:logfile => LOGFILE)))
    segmentator_init(track, merge(_alg, get(_alg, :segmentator_init, Dict()), Dict(:logfile => LOGFILE)))
    criterion_init(merge(_alg, get(_alg, :criterion_init, Dict()), Dict(:logfile => LOGFILE)))
    penalizer_init(merge(_alg, get(_alg, :penalizer_init, Dict()), Dict(:logfile => LOGFILE)))
    optimizer_init(merge(Dict(:points=>track, :group_centers=>rcandidate, :group_centerline=>result, :logfile=>LOGFILE), _alg))


    ## Optimization
    _fitness, _rcandidate, _tcandidate, _result = optimize()


    ## Plot the solution
    # Note: Dynamic plotting is not supported yet. This is the stable configuration I was using before.
    if get(_alg, :plot, false)
        img = map(TRACK_BITMAP.bitmap') do p p ? colorant"#eed" : colorant"white" end
        @gp :map Gnuplot.palette(:gray1) "set size ratio -1" :-
        @gp :- recipe(img, "dx=$(TRACK_BITMAP.p_step[1]) dy=$(TRACK_BITMAP.p_step[2]) origin=($(TRACK_BITMAP.p_min[1]), $(TRACK_BITMAP.p_min[2]))") :-
        @gp :- append!(_result[:, 1], _result[1, 1]) append!(_result[:, 2], _result[1, 2]) "with lines linecolor 1 title ''" :-
        @gp :- [ _rc[1] for _rc in _rcandidate ] [ _rc[2] for _rc in _rcandidate ] "linecolor 1 linetype 1 pointsize 1.5 title ''" :-

        if fileformat !== ""
            save(term="pngcairo size 1280, 720 fontscale 0.8", output=replace(fileformat, "%d" => loop_i[1]) * "-$(get(_alg, :algorithm, "")).png")
        else
            # Figure out a way to pause the execution here, until the figure is closed.
            @gp :-
        end
    end

    @show _fitness
    ## End parts
    if fileformat !== ""
        # Show all results of optimize function (log only)
        write(LOGFILE, "#fitness:$(_fitness)\n")
        write(LOGFILE, "#rcandidate:$(string(_rcandidate))\n")
        __tcandidate = [ _tcandidate[i, :] for i in range(1, length(_tcandidate[:, 1])) ]
        write(LOGFILE, "#tcandidate:$(string(__tcandidate))\n")
        __result = [ _result[i, :] for i in range(1, length(_result[:, 1])) ]
        write(LOGFILE, "#trajectory:$(string(__result))\n")
    end
    end #@elapsed
    # Show up time elapsed
    write(LOGFILE, "time:$(step_time)\n")
    write(LOGFILE, "==============\n")

    # Close file if opened
    if fileformat !== ""
        close(LOGFILE)
    end


    # Store only better solution for next steps of the cascade
    if _fitness < fitness
        return (_fitness, _rcandidate, _tcandidate, _result)
    else
        return loop_output
    end
end

function loop_cascade_run(fileformat, notification, loop_i, loop_output, conf)
    # Using global variables to avoid segfaults
    global _VALID_POINTS, _START_POINTS

    track = _VALID_POINTS
    initline = _START_POINTS

    # Cascade timing
    cascade_time = @elapsed begin

    # Initial solution
    fitness = 10000000
    result = initline
    rcandidate = initline
    tcandidate = [0.5 for _ in 1:size(initline, 1)*2]

    # Update logging file
    _fileformat = ""
    if fileformat !== ""
        # TODO: Missing leading zeroes.
        # Printf.format(Printf.Format(fileformat), loop_i) * @sprintf("-%%0%dd", length(string(length(conf[:cascade]))))
        _fileformat = replace(fileformat, "%d" => loop_i) * "-%d"
    end

    # Printf.format(Printf.Format(notification), loop_i) * @sprintf(" Running step %%d/%d", length(conf[:cascade]))
    notification = replace(notification, "%d" => loop_i) * " Running step %d/$(length(conf[:cascade]))"


    # Run cascade
    cascade_output = nothing
    for i in 1:length(conf[:cascade])
        cascade_output = cascade_run(
            track,
            _fileformat,
            notification,
            (i, conf[:cascade][i]),
            cascade_output === nothing ? (fitness, rcandidate, tcandidate, result) : cascade_output,
            conf
        )
    end

    end #@elapsed


    if fileformat !== ""
        open(replace(fileformat, "%d" => loop_i) * ".log", "w") do logfile
            write(logfile, "timeA:$(cascade_time)\n")
        end
    else
        @info "timeA:$(cascade_time)"
    end

    if loop_output === "" || cascade_output[1] < loop_output[1]
        return cascade_output
    else
        return loop_output
    end
end

function variate_run(fileformat, notification, loop_i, loop_output, conf)

    # Group timing
    variate_time = @elapsed begin

    # Local variables
    _i = loop_i[1]
    _param = loop_i[2][1]
    _value = loop_i[2][2]

    # Update logging file
    if fileformat !== ""
        # Fill group count, and add format for number of loops
        # TODO: Missing leading zeroes.
        #fileformat = Printf.format(Printf.Format(fileformat), _value) * @sprintf("-%%0%dd", length(string(CONFIGURATION[:loops])))
        fileformat = replace(fileformat, "%d" => _value) * "-%d"
    end

    # Update notification
    # Fill loop index, group count and prepare loops progress
    # TODO: Missing leading zeroes.
    # notification = Printf.format(Printf.Format(notification), _i, _value, _param) * @sprintf(" [%%d / %d]", CONFIGURATION[:loops])
    notification = replace(notification, "%dcurrent" => _i, "%svalue" => _value, "%sparam" => _param) * " [%d / $(CONFIGURATION[:loops])]"

    cascade_output = ""

    ## Loop cascade
    for i in 1:conf[:loops]
        cascade_output = loop_cascade_run(
            fileformat,
            notification,
            i,
            cascade_output,
            merge(conf, Dict(Symbol(_param) => _value))
        )
    end

    end #@elapsed

    @printf("Variating %s %s finished in %fs.\n", _param, _value, variate_time)

    if loop_output === "" || cascade_output[1] < loop_output[1]
        return cascade_output
    else
        return loop_output
    end
end

function execute(START_POINTS=nothing, VALID_POINTS=nothing)

    global CONFIGURATION
    global _START_POINTS, _VALID_POINTS

    println("Starting Matryoshka.jl version ", Version(Matryoshka))

    # Overall time
    overall_time = @elapsed begin

    # Load data about the track
    if START_POINTS === nothing
        START_POINTS = npzread(CONFIGURATION[:start_points])
    end
    if VALID_POINTS === nothing
        VALID_POINTS = npzread(CONFIGURATION[:valid_points])
    end

    _START_POINTS = START_POINTS
    _VALID_POINTS = VALID_POINTS

    # Logging file format
    fileformat = ""
    if haskey(CONFIGURATION, :prefix)
        fileformat = string(CONFIGURATION[:prefix])
    end

    # Notification about progress
    notification = ""

    solution = nothing

    # Identify and prepare variating variable
    if haskey(CONFIGURATION, :variate) && haskey(CONFIGURATION, Symbol(CONFIGURATION[:variate]))
        # Currently unsupported in C.
        param = CONFIGURATION[:variate]
        values = CONFIGURATION[Symbol(param)]

        # Force list
        if !(values isa Array{<:Number})
            values = [values]
        end

        # Convert to tuples
        tvalues = [(param, value) for value in values]

        # Add variate to the file format
        if fileformat !== ""
            if all(typeof(_value) == Int for _value in values)
                # TODO: Missing leading zeroes.
                fileformat = fileformat * "-%d" # @sprintf("-%%0%dd", length(string(maximum(values))))
            else
                fileformat = fileformat * "-%s"
            end
        end

        # ... and also to the notification
        # notification = notification * @sprintf("{%%d / %d (%%s %%s)}", length(values))
        notification = notification * "{%dcurrent / $(length(values)) (%svalue %sparam)}"

        ## And variate the parameter
        solution = ""
        for i in 1:length(tvalues)
            solution = variate_run(
                fileformat,
                notification,
                (i, tvalues[i]),
                solution,
                CONFIGURATION
            )
        end

    else
        # Skip to the loop
        # Update logging file
        if fileformat !== ""
            # TODO: Missing leading zeroes.
            fileformat = fileformat * "-%d" #@sprintf("-%%0%dd", length(string(CONFIGURATION[:loops])))
        end

        # Update notification
        notification = notification * "[%d / $(CONFIGURATION[:loops])]" #@sprintf("[%%d / %d]", CONFIGURATION[:loops])

        ## Loop cascade
        solution = ""
        for i in 1:CONFIGURATION[:loops]
            solution = loop_cascade_run(
                fileformat,
                notification,
                i,
                solution,
                CONFIGURATION
            )
        end
    end

    end #@elapsed

    @info "Optimization finished in $(overall_time)s."

    return solution
end

if (abspath(PROGRAM_FILE) == @__FILE__)
    println(execute())
end
