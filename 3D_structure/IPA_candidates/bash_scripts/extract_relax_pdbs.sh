#!/bin/bash

cd /scicore/home/zavolan/zhu0006/3D_structure/IPA_candidates/alphafold_res

folders=("bad_10_iso" "bad_10_ref" "bad_11_iso" "bad_11_ref" "bad_13_iso" "bad_13_ref" "good_0_iso" "good_0_ref" "good_1_iso" "good_1_ref" "good_2_iso" "good_2_ref" "good_3_iso" "good_3_ref" "good_4_iso" "good_4_ref" "good_5_iso" "good_5_ref" "good_6_iso" "good_6_ref" "good_7_iso" "good_7_ref" "good_8_iso" "good_8_ref" "good_9_iso" "good_9_ref")

for folder in "${folders[@]}"; do
    echo "Processing $folder"
    cd "$folder"
    relax_out="relax.out"
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta.binary.linux.release-371/main/source/bin/extract_pdbs.static.linuxgccrelease -in:file:silent "$relax_out" -in:file:tagfile lowest.tag
    cd ..
done

echo "Done!"
