#!/bin/bash

awk '/^ATOM/ {split($0, a, " "); if (a[5]=="A") print}'