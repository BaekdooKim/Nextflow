#!/bin/bash

rm input_list.txt > /dev/null 2>&1
for i in $(pwd)/*; do if [ -f "$i" ] && [ "$(basename ${i##*.})" = "gtf" ]; then printf "$(basename ${i%.*})\t$i\n" >> input_list.txt; fi done