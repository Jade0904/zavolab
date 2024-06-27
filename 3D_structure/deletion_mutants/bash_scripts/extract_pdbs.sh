#!/bin/bash

# to extract pdb files out of relax outputs

cd /scicore/home/zavolan/zhu0006/3D_structure/deletion_mutants/alphafold_res

del_folders=("del2" "del3" "del5" "del50" "del51" "del52" "del62" "del64" "del66" "del67" "del68" "del69" "del70" "del71" "del72")

for folder in "${del_folders[@]}"; do
    echo "Processing $folder"
    cd "$folder"
    relax_out="af_${folder}.out"
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta.binary.linux.release-371/main/source/bin/extract_pdbs.static.linuxgccrelease -in:file:silent "$relax_out" -in:file:tagfile lowest.tag
    cd ..
done

echo "Done!"
