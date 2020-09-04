#!/bin/bash

#Script not tried yet but should be close to working
grep "event" coo400b-9.log > file1   #Enter GBD log file name and desired output file name
grep "event" coo400b-9.log > file3
grep "event" coo400b-9.log > file3

/home/tim/bin/scripts/kon_combine3 file1 file2 file3 500 > konmulti

awk -F, 'length($0) != 1 { print }' konmulti > tout_clean
