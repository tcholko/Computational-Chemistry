#!/bin/bash

grep "event" COOe11.log > t_COOe11   #Enter GBD log file name and desired output file name

/home/tim/bin/scripts/kon t_COOe11 105 > kon_COOe11   #Enter number of binders from GBD log
