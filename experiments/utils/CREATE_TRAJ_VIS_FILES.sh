i=1 rostopic echo -p -b "$1".bag "/trajectory" | sed 1d | { while read line; do echo -n "$line" | cut -d, -f5- | sed 's/,/\n/14;P;D' | cut -d, -f2,3 > "$1".trajectory."$i".csv; i=$((i+1)); done }
