#!/bin/sh
cd src
for f in $(ls ../graphs) ; do
    echo "-------------------processing $f---------------------"
    python3.7 -m app.main -g ../graphs/$f
done
cd ..
