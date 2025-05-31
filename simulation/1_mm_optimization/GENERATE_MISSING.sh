find . -name "*ia.log" -or -name "*ka.log" -or -name "*in.log" | xargs -I'{}' bash -c "test -f '{}.lap_time' || python3 obtain_lap_time2.py '{}'"
