i=1 rostopic echo -p -b "$1".bag "/odom" | sed 1d | cut -d, -f6,7 > "$1".odom.csv
