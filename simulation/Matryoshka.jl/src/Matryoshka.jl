module Matryoshka

#using Revise

include("parameter.jl")
include("utils.jl")
include("transform.jl")
include("selector.jl")
include("interpolator.jl")
include("segmentator.jl")
include("criterions.jl")
include("optimizer.jl")
include("penalizer.jl")
include("run.jl")

#=
export trajectory_resample
export P_sel
export profile_compute
export trajectory_closest_index
export trajectory_time_resample
export jazar_profile_compute
export P_cri
export optimizer_init
export set_maps
=#

export configuration_load
export execute

export run_ma

Base.@ccallable function run_ma(start::Ptr{Cdouble}, valid::Ptr{Cdouble})::Cint
    @show start
    return 0
end


export load_configuration

Base.@ccallable function load_configuration(filename::Cstring)::Cint
    @show filename
    @show unsafe_string(filename)

    ret = configuration_load(unsafe_string(filename))

    return ret ? 0 : 1
end

export execute_ma

Base.@ccallable function execute_ma()::Cint
    k = execute()
    @show k
    return 0
end

export execute_ma_points

Base.@ccallable function execute_ma_points(line::Ptr{Cdouble}, line_length::Cint, track::Ptr{Cdouble}, track_length::Cint)::Cint
    @show line
    @show line_length
    @show track
    @show track_length

    _line = reshape(unsafe_wrap(Vector{Float64}, line, Base.cconvert(Int, line_length)*2; own=false), :, 2)
    @show _line

    _track = reshape(unsafe_wrap(Vector{Float64}, track, Base.cconvert(Int, track_length)*2; own=false), :, 2)

    @show execute(_line, _track)

    return 0
end

end # module Matryoshka
