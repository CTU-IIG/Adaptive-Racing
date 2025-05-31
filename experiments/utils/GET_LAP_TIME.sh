rostopic echo -p -b "$1".bag "/lap_time" | sed 1d  > "$1".lap_time.csv
