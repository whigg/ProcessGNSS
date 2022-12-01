#!/bin/bash

echo The following script converts .T00 to .dat to Rinex

for FILE in *; do ./teqc -O.obs L1+l2+ca+p2+s1+s2 $FILE > $FILE".obs" ; done