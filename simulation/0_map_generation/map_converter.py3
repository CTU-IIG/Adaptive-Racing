#!/usr/bin/env python3
# map_converter.py
"""Convert bitmap map data into race trajectory '.csv' file.

This file is based on:
https://github.com/CL2-UWaterloo/Raceline-Optimization/map_converter.ipynb
"""
######################
# Imports & Globals
######################

import argparse
import os
import sys
import numpy
from PIL import Image
import yaml
import scipy.ndimage
import skimage.morphology
import matplotlib.pyplot


# Increase recursion limit for DFS
sys.setrecursionlimit(20000)


######################
# Argparse
######################

PARSER = argparse.ArgumentParser(
    prog = "map_converter",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description = """
Convert bitmap map data into race trajectory '.csv' file.
    """,
)

PARSER.add_argument(
    "map_name",

    help = (
        "Name of the map. Expecting files "
        "'MAP_NAME.pgm' and 'MAP_NAME.yaml'."
    ),
    type = str,
)

PARSER.add_argument(
    "--threshold",

    help = "Threshold value for edt.",
    default = 0.17,
    type = float,
)

PARSER.add_argument(
    "-v",

    dest = "verbose",
    help = "Give more output.",
    default = 0,
    action = "count",
)


######################
# Utilities
######################

def printv(*a, **kw) -> None:
    """Print only when verbose is set."""
    if args.verbose:
        print (*a, **kw)


def savevv(nd, filename: str) -> None:
    """Save an image from ndarray when verbose is set to >=2."""
    if args.verbose >= 2:
        matplotlib.pyplot.figure()
        matplotlib.pyplot.imshow(nd, cmap = "gray", origin = "lower")
        matplotlib.pyplot.savefig(filename)


######################
# Map converter
######################

def map_converter(map_name, map_image, map_yaml):
    """Convert the map into '.csv' file.

    Arguments
    ---------
    map_name: str
        base name of the map
    map_image: str
        path to the image of the track
    map_yaml: str
        path to the configuration of the track
    """
    # 1. Load configuration
    cfg = yaml.safe_load(
        open(map_yaml, "r").read()
    )
    printv("Loaded configuration:", cfg)

    # 2. Load the image
    raw_map_image = Image.open(map_image).transpose(Image.FLIP_TOP_BOTTOM)
    printv("Loaded image:", raw_map_image)
    raw_map_img = numpy.array(raw_map_image).astype(numpy.float64)
    printv("Numpy array shape:", raw_map_img.shape)
    # printv("Color statistics:", numpy.unique(raw_map_img, return_counts=True))

    # 3. Convert the image to grayscale
    map_img = raw_map_img.copy()
    map_img[map_img <= 255. * (1. - cfg["free_thresh"])] = 0
    map_img[map_img > 0.0] = 1
    printv("Free space ratio:", numpy.sum(map_img > 0.0) / numpy.size(map_img))
    savevv(map_img, "grayscale.png")

    # 4. Calculate Euclidean Distance Transform
    dist_transform = scipy.ndimage.distance_transform_edt(map_img)
    savevv(dist_transform, "edt.png")

    # 5. Threshold edt to create binary image
    centers = dist_transform > args.threshold * dist_transform.max()
    savevv(centers, "edt_centers.png")

    # 6. Obtain centerline
    centerline = skimage.morphology.skeletonize(centers)
    savevv(centerline, "centerline.png")

    # 7. Obtain centerline with distance / track width
    centerline_dist = numpy.where(centerline, dist_transform, 0)
    savevv(centerline_dist, "centerline_dist.png")

    # 8. Select starting point for the raceline data
    NON_EDGE = 0.0
    left_start_y = map_img.shape[0] // 2 - 120
    left_start_x = 0

    while (centerline_dist[left_start_y][left_start_x] == NON_EDGE):
        left_start_x += 1

    printv(
        "Starting position for left edge: %d %d"
        % (left_start_x, left_start_y)
    )

    # 9. Run DFS to sort the raceline
    visited = {}
    centerline_points = []
    track_widths = []
    DIRECTIONS = [
        (0, -1), (-1, 0), (0, 1), (1, 0), (-1, 1), (-1, -1), (1, 1), (1, -1)
    ]
    starting_point = (left_start_x, left_start_y)

    def dfs(point):
        if (point in visited):
            return

        visited[point] = True
        centerline_points.append(numpy.array(point))
        track_widths.append(
            numpy.array([
                centerline_dist[point[1]][point[0]],
                centerline_dist[point[1]][point[0]]
            ])
        )

        for direction in DIRECTIONS:
            _x = point[0] + direction[0]
            _y = point[1] + direction[1]

            if (
                centerline_dist[_y][_x] != NON_EDGE
                and (_x, _y) not in visited
            ):
                dfs((_x, _y))

    dfs(starting_point)

    # 10. Convert the data
    widths = numpy.array(track_widths)
    waypoints = numpy.array(centerline_points)
    printv("Track widths shape:", widths.shape)
    printv("Waypoints shape:", waypoints.shape)

    data = numpy.concatenate((waypoints, widths), axis = 1)
    data *= cfg["resolution"]
    data += numpy.array([
        cfg["origin"][0], cfg["origin"][1], 0, 0
    ])

    # 11. Save the data
    with open(map_name + ".csv", "wb") as f:
        numpy.savetxt(
            f, data, fmt = "%0.4f", delimiter = ",",
            header = "x_m,y_m,w_tr_right_m,w_tr_left_m"
        )


######################
# Main
######################

if __name__ == "__main__":

    # Obtain arguments
    args = PARSER.parse_args()


    # Check files
    map_name = args.map_name
    map_image = map_name + ".pgm"
    map_yaml = map_name + ".yaml"

    if not os.path.exists(map_image):
        raise IOError(
            "Map not found! File '%s' does not exist." % map_image
        )

    if not os.path.exists(map_yaml):
        raise IOError(
            "Configuration not found! File '%s' does not exist." % map_yaml
        )

    map_converter(map_name, map_image, map_yaml)
