#!/bin/bash

cd /scicore/home/zavolan/zhu0006/3D_structure/indel/indel_examples_pdb/alphafold_res

folders=("4pe5_A" "7mpi_L" "1vso_A" "5fxj_A" "6k4y_F" "6tac_A" "8e31_C" "1yae_A" "2z2a_A" "1c8r_A" "8e3b_C" "6b2f_A" "5vot_A" "5l1h_A" "6hoo_A" "4awh_A" "4ilg_A" "6hb8_A" "6ta2_A" "2z2b_A" "8sar_B" "7b7d_L" "4i43_A" "8saq_B" "1c8s_A" "6n4c_F" "6bhl_A" "4e5e_A")

for folder in "${folders[@]}"; do
    echo "Processing $folder"
    cd "$folder"
    relax_out="relax.out"
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta.binary.linux.release-371/main/source/bin/extract_pdbs.static.linuxgccrelease -in:file:silent "$relax_out" -in:file:tagfile lowest.tag
    cd ..
done

echo "Done!"
