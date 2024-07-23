#!/bin/bash 
#SBATCH --job-name=relax
#SBATCH --time=7-00:00:00
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=16
#SBATCH --partition=a100 
#SBATCH --gres=gpu:1

# this script is an example
cd /scicore/home/zavolan/zhu0006/3D_structure/indel/indel_examples_pdb/alphafold_res

ids=("6tac_A")

for id in "${ids[@]}"; do
    echo "Processing $id"
    cd "$id"
    ls ranked_?.pdb > models.list
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta.binary.linux.release-371/main/source/bin/relax.static.linuxgccrelease -in:file:l models.list -relax:constrain_relax_to_start_coords 1 -nstruct 10 -out:file:silent relax.out -out:file:scorefile relax.sc > relax.log
    cd ..
done


echo "Done!"