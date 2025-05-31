ls -1 | { while read line; do cd "$line"; make; cd ..; done }
