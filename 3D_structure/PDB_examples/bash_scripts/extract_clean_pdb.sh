#!/bin/bash

directory="original_files"

for file in "$directory"/*.pdb; do
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta.binary.linux.release-371/main/tools/protein_tools/scripts/clean_pdb.py "$file" 