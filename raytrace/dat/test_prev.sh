#!/bin/sh
echo 'Apply AA test too'
for file in *.dat 
do 
    echo "RT file: $file, press enter"
    ls -l "$file"
    read z
    if [ "$z" = "q" ] || [ "$z" = "Q" ]
	then
	    exit 1
	fi
    ../rays.fast -i "$file" -f -J
done

