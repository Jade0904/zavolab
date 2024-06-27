#!/bin/bash 
#SBATCH --job-name=del-Relax_Rosetta
#SBATCH --time=14:00:00
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=8 
#SBATCH --partition=a100 
#SBATCH --gres=gpu:1
 
cd /scicore/home/zavolan/zhu0006/3D_structure/deletion_mutants/alphafold_res

for dir in del*/; do
    echo "Processing $dir"
    cd "$dir"

    ls ranked_?.pdb > "${dir%/}_models.list"
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/relax.static.linuxgccrelease -in:file:l "${dir%/}_models.list" -relax:constrain_relax_to_start_coords 1 -nstruct 10 -out:file:silent "af_${dir%/}.out" -out:file:scorefile "af_${dir%/}.sc" > "af_${dir%/}.log"

    cd ..
done