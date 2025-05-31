#!/usr/bin/env python
# process.py
# Process bag file and generate graph for the RAL paper.
# Source: /media/klapajar/MyBook/julia_experiment/experiment/ral_final

import rospy
import rosbag
import sys

import matplotlib.pyplot

import numpy, math

# GridCells
from nav_msgs.msg import GridCells
# an array of cells in a 2D grid
#: Header header
#: float32 cell_width
#: float32 cell_height
#: geometry_msgs/Point[] cells

# Odometry
from nav_msgs.msg import Odometry
# This represents an estimate of a position and velocity in free space. 
# The pose in this message should be specified in the coordinate frame given by header.frame_id.
# The twist in this message should be specified in the coordinate frame given by the child_frame_id
#: Header header
#: string child_frame_id
#: geometry_msgs/PoseWithCovariance pose
#: geometry_msgs/TwistWithCovariance twist

# Trajectory
from autoware_auto_msgs.msg import Trajectory
#: std_msgs/Header header
#: autoware_auto_msgs/TrajectoryPoint[] points

# TrajectoryPoint
from autoware_auto_msgs.msg import TrajectoryPoint
# Point of a Trajectory
# Zero heading corresponds to the positive x direction in the given frame.
# All fields should default to 0.
# As noted in the original message, it is possible that more fields will
# be added overtime (e.g., jerk).
# Note: We use type duration instead of std_msgs/Duration to make the API
# consistent. The only difference is, that 'time_from_start.nsecs' is int32,
# whereas it is uint32 in ROS2.
#: duration time_from_start                # ROS duration, .secs and .nsecs
#: geometry_msgs/Pose pose
#: float32 longitudinal_velocity_mps       # [m.s^-1]
#: float32 lateral_velocity_mps            # [m.s^-1]
#: float32 acceleration_mps2               # [m.s^-2]
#: float32 heading_rate_rps                # [rad.s^-1]
#: float32 front_wheel_angle_rad           # [rad]
#: float32 rear_wheel_angle_rad            # [rad]

def setPlot():
    matplotlib.pyplot.xlim([-3.5, 9])
    matplotlib.pyplot.ylim([-2, 7])
    matplotlib.pyplot.plot([-3.5, 9, 9, -3.5], [-2, -2, 7, 7], color = "black")
    matplotlib.pyplot.axis("equal")


def savePlot(filename):
    matplotlib.pyplot.savefig(filename, bbox_inches="tight", dpi=300, pad_inches=0)
    matplotlib.pyplot.close()


import colorsys
def scale_lightness(rgb, scale_l):
    """Scale the lightness of a color.

    Source:
    https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib

    Taken from 'ng_graph'.
    """
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)


if len(sys.argv) < 2:
    print >>sys.stderr, "Missing a name of the bag."
    exit()

bag = rosbag.Bag(sys.argv[1])

GRID = None

odom_x = []
odom_y = []

traj_x = []
traj_y = []
traj_xy = []
traj_t = []
traj_v = []

trajectory_i = 0

for topic, msg, t in bag.read_messages(topics = ["/map/map_grid_cells", "/odom", "/trajectory/autoware"]):
    if topic == "/odom":
        odom_x.append(msg.pose.pose.position.x)
        odom_y.append(msg.pose.pose.position.y)
    elif topic == "/trajectory/autoware":
        trajectory_i += 1
        traj_x = [ tpoint.pose.position.x for tpoint in msg.points ]
        traj_y = [ tpoint.pose.position.y for tpoint in msg.points ]
        traj_t = [ tpoint.time_from_start.to_sec() for tpoint in msg.points ]
        traj_v = [ numpy.hypot(tpoint.longitudinal_velocity_mps, tpoint.lateral_velocity_mps) for tpoint in msg.points ]

        #if trajectory_i == 19:
        #    odom_x = []
        #    odom_y = []
        #elif trajectory_i > 19:
        #    break
    else:
        GRID = msg.cells

#print GRID
print(traj_t)


traj_x.append(traj_x[0])
traj_y.append(traj_y[0])
traj_xy = numpy.vstack((traj_x, traj_y)).T.reshape(-1, 1, 2)

# Inspired by:
# https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm, LogNorm, LinearSegmentedColormap

base_color = matplotlib.colors.ColorConverter.to_rgb("#00ff00")#'#2ca02c')
tabgreen_cmap = LinearSegmentedColormap.from_list("tabgreen", [ (0.0, scale_lightness(base_color, 0.4)), (1.0, scale_lightness(base_color, 1.6)) ])

norm = matplotlib.pyplot.Normalize(min(traj_v),max(traj_v))
#norm = LogNorm(min(traj_v), max(traj_v))
lc = LineCollection(numpy.concatenate([traj_xy[:-1], traj_xy[1:]], axis=1),
    #colors = [ scale_lightness(base_color, ((_v - min(traj_v))/(max(traj_v) - min(traj_v)))) for _v in traj_v ]
    cmap=tabgreen_cmap,#"viridis",
    norm=norm,
    linewidths = numpy.asarray(traj_v)*1.7 # Norm: Ale nevypada to tak dobre numpy.asarray(norm(traj_v))*5
)
lc.set_array(numpy.asarray(traj_v))
#lc.set_linewidth(1.2)

matplotlib.pyplot.gca().add_collection(lc)
#matplotlib.pyplot.colorbar(matplotlib.pyplot.gca().add_collection(lc))
setPlot()
savePlot(sys.argv[1] + ".intensity.png")

#exit()


if False:
    isopoints = []

    for _tp in numpy.linspace(0, 13.276, 13, endpoint = False):
        first = None
        second = None

        for _i, _t in enumerate(traj_t):
            if _t == _tp:
                isopoints.append((traj_x[_i], traj_y[_i]))
                break
            elif _t < _tp:
                first = _i
            else:
                second = _i
                _dt = _tp - traj_t[first]
                _dv = traj_v[second] - traj_v[first]
                _ds = numpy.hypot(traj_x[second] - traj_x[first], traj_y[second] - traj_y[first])
                _dx = traj_x[second] - traj_x[first]
                _dy = traj_y[second] - traj_y[first]

                _k = (2 * traj_v[first] * _dt) / (2 * _ds - _dv * _dt)


                _vk = traj_v[first] + _dv * _k
                _alpha = math.atan2(_dy, _dx)
                _sk = ( (traj_v[first] + _vk) / 2 ) * _dt
                isopoints.append(
                    (
                        traj_x[first] + _sk * numpy.cos(_alpha),
                        traj_y[first] + _sk * numpy.sin(_alpha)
                    )
                )

                break
            
    matplotlib.pyplot.scatter([ p[0] for p in isopoints ], [ p[1] for p in isopoints ], color = "black")
    setPlot()
    savePlot(sys.argv[1] + ".isopoints.png")

    exit()

#_x = [ point.x for point in GRID ]
#_y = [ point.y for point in GRID ]


#matplotlib.pyplot.scatter(_x, _y, color = "lightgray")
#matplotlib.pyplot.plot(odom_x, odom_y)
#matplotlib.pyplot.plot(traj_x, traj_y)
#setPlot()
#savePlot(sys.argv[1] + ".track.png")


matplotlib.pyplot.plot(odom_x, odom_y, color = "tab:blue")
setPlot()
savePlot(sys.argv[1] + ".odometry.png")


matplotlib.pyplot.plot(traj_x, traj_y, ":", color = "tab:orange")
setPlot()
savePlot(sys.argv[1] + ".trajectory.png")



exit()
bag = rosbag.Bag("final_ftg.bag")

ftg_x = []
ftg_y = []

for topic, msg, t in bag.read_messages(topics = ["/odom"]):
    if topic == "/odom":
        ftg_x.append(msg.pose.pose.position.x)
        ftg_y.append(msg.pose.pose.position.y)


matplotlib.pyplot.plot(ftg_x, ftg_y, color = "tab:green")
setPlot()
savePlot("ftg_old_map.png")

