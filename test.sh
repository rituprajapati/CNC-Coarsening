#!/bin/bash

for threads in 1 2 4 8 16 32 64
do
	export CNC_NUM_THREADS=$threads
	for max_level in 20
	do
		for thresh in 1e-16
		do
			./main $max_level $thresh $threads
		done
	done
done
