rostopic echo -p -b "$1".bag "/trajectory/tracking_error" | sed 1d | cut -d, -f2- > "$1".tracking_error.csv
