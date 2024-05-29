#!/bin/bash

cd /scicore/home/zavolan/zhu0006/3D_structure/deletion_mutants/alphafold_res

del_folders=("del2" "del3" "del5" "del50" "del51" "del52" "del62" "del64" "del66" "del67" "del68" "del69" "del70" "del71" "del72")

for folder in "${del_folders[@]}"; do
    echo "Processing $folder"
    cd "$folder"
    basename=$(head -n 1 "lowest.tag")
    lowest="${basename}.pdb"
    score_perRes="${basename}_perRes.sc"
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/per_residue_energies.static.linuxgccrelease -in:file:s "$lowest" -out:file:silent "$score_perRes"
    cd ..
done

echo "Done!"