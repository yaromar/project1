#!/bin/bash

if [ "$#" -eq 1 ]; then
	curl -O files.rcsb.org/view/$1.pdb
elif [ "$#" -eq 2 ]; then
	curl -o $2/$1.pdb files.rcsb.org/view/$1.pdb
fi
