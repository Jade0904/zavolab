#!/bin/bash

awk '/^ATOM/ {split($0, a, " "); if (a[5]=="A") print}' original_files/4i43.pdb > chain_files