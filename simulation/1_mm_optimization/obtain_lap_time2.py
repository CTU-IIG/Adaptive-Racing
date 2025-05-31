#!/usr/bin/env python3
# obtain_lap_time2.py
"""Script for obtaining lap time for ng_trajectory files.

Inspired and partially reused from 'ng_plot' and 'obtain_lap_time.py'.
This reads data from the log and interprets them in a smarter way.
"""
######################
# Imports & Globals
######################

import argparse
import ast
import os

import ng_trajectory.interpolators.cubic_spline as cspline
import ng_trajectory.criterions.profile.main as profile
import ng_trajectory.criterions.profile.profiler as profiler
from ng_trajectory.segmentators.utils import mapCreate

import numpy
import json

# Progress
import tqdm


# Global variables
VALID_POINTS = None


######################
# Argparse
######################

PARSER = argparse.ArgumentParser(
    prog = "obtain_lap_time",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description = """
Script for obtaining lap time for ng_trajectory files.

Inspired by 'ng_plot'.
    """,
)

PARSER.add_argument(
    "input_file",

    # nargs = 1,
    help = "Path to the ng_trajectory file.",
    type = argparse.FileType("r"),
)

PARSER.add_argument(
    "-d",

    dest = "display",
    help = "When given, show the result instead of storing it into a file.",
    action = "store_true",
)

PARSER.add_argument(
    "-v",

    dest = "verbose",
    help = "Give more output.",
    action = "store_true",
)

# Dynamicaly build up the parameters for profiler
_PROFILE_SETTINGS = {
    "overlap": 100,
    "_mu": 0.3,
    "v_lim": 8.0,
    "a_acc_max": 10.0,
    "a_break_max": 1.8
}

_PROFILE_PARAMETERS = [
    "overlap", "_mu", "_g", "_m", "_ro", "_A", "_cl", "v_0", "v_lim",
    "a_acc_max", "a_break_max", "_lf", "_lr"
]

PROFILER_GROUP = PARSER.add_argument_group(
    "profiler arguments",
    "Parameters used for the profiler."
)

for name, parameter in profile.P.parameters.items():
    if name not in _PROFILE_PARAMETERS:
        continue

    PROFILER_GROUP.add_argument(
        "--" + name,

        dest = name,
        help = parameter.description + " (%(type)s, default: %(default)s)",
        default = (
            _PROFILE_SETTINGS.get(name) if name in _PROFILE_SETTINGS
            else parameter.default
        ),
        type = parameter.type,
        metavar = "VALUE",
    )

    if name not in _PROFILE_SETTINGS:
        _PROFILE_SETTINGS.setdefault(name, parameter.default)


######################
# Utilities
######################

def write_result(value):
    """Write the result into a file or to the terminal."""
    if args.display:
        print(value)
    else:
        with open(result_file, "w") as f:
            f.write(value)



# python 2
# https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists
def gen_dict_extract(key, var):
    if hasattr(var, 'iteritems'):  # hasattr(var,'items') for python 3
        for k, v in var.iteritems():  # var.items() for python 3
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result

    if hasattr(var, 'items'):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result


######################
# Main
######################

if __name__ == "__main__":

    # Obtain arguments
    args = PARSER.parse_args()

    if args.verbose:
        print ("Received arguments:")
        print (
            "\n".join([
                "\t%s: %s"
                % (key, value) for key, value in args._get_kwargs()
            ])
        )

    # Check for result
    result_file = args.input_file.name + ".lap_time"

    if os.path.exists(result_file):
        print (open(result_file, "r").readline())
        exit(0)


    # Obtain configuration
    try:
        cfg = ast.literal_eval(
            args.input_file.readline()
        )
    except Exception:
        if args.verbose:
            print ("Unable to detect in-place configuration.")

        # Try to do it ourselves
        for f in os.scandir(os.path.dirname(args.input_file.name)):
            if f.name.endswith(".json"):
                cfg_path = os.path.join(
                    os.path.dirname(args.input_file.name),
                    f.name
                )

                if args.verbose:
                    print (
                        "Loading configuration from '%s'." % cfg_path
                    )

                cfg = json.load(
                    open(cfg_path, "r")
                )

                break

    map_file = os.path.join(
        os.path.dirname(args.input_file.name),
        list(gen_dict_extract("valid_points", cfg))[0]
    )

    optimizer = list(gen_dict_extract("algorithm", cfg))[0]
    criterion = list(gen_dict_extract("criterion", cfg))[0]

    if args.verbose:
        print ("Interpreted configuration:")
        print ("valid_points:", map_file)
        print ("criterion:", criterion)
        print ("optimizer:", optimizer)


    # Obtain solutions from logs
    solutions = []
    trajectory = None
    valid = False
    final = None

    for line in tqdm.tqdm(args.input_file, desc = "Loading", leave = False):
        if line.startswith("correct:"):
            valid = True
        if line.startswith("solution:"):
            solutions.append(
                json.loads(line[9:].rstrip())
            )
        if line.startswith("final:"):
            final = float(line[6:].rstrip())
        if line.startswith("#trajectory:") and optimizer == "matryoshka_julia":
            # Julia results are different when interpolated
            # So by default, use the trajectory data instead of solution
            trajectory = numpy.asarray(
                json.loads(line[len("#trajectory:"):].rstrip())
            )

    if args.verbose:
        print ("Found solutions: %d" % len(solutions))

    if len(solutions) < 1 or not valid:
        write_result("nan")
        exit(0)

    elif criterion == "profile":
        write_result("%f" % final)
        exit(0)

    # Obtain map from the file
    VALID_POINTS = numpy.load(map_file)
    mapCreate(VALID_POINTS)

    # Prepare criterion
    if True:
        for param in _PROFILE_SETTINGS.keys():
            _PROFILE_SETTINGS[param] = getattr(args, param)

        if args.verbose:
            print ("Setting profiler to:")
            print (_PROFILE_SETTINGS)

        profile.init(**_PROFILE_SETTINGS)

    if trajectory is None:
        solution = solutions[0]

        line = numpy.asarray(solution)
        iline = cspline.interpolate(line)

    else:
        iline = trajectory

    _v, _a, _t = profiler.profileCompute(iline, 400, lap_time = True)

    if args.verbose:
        print ("Estimated lap time: %f" % _t[-1])

        if not args.display:
            write_result("%f" % _t[-1])
    else:
        write_result("%f" % _t[-1])
