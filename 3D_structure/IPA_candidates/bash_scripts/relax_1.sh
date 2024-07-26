#!/bin/bash 
#SBATCH --job-name=relax_1
#SBATCH --qos=1day
#SBATCH --time=1-00:00:00
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=16
#SBATCH --partition=a100 
#SBATCH --gres=gpu:1

cd /scicore/home/zavolan/zhu0006/3D_structure/IPA_candidates/alphafold_res

ids=("bad_13_iso" "bad_13_ref")

for id in "${ids[@]}"; do
    echo "Processing $id"
    cd "$id"
    ls ranked_?.pdb > models.list
    /scicore/home/zavolan/zhu0006/3D_structure/Rosetta/rosetta.binary.linux.release-371/main/source/bin/relax.static.linuxgccrelease -in:file:l models.list -relax:constrain_relax_to_start_coords 1 -nstruct 5 -out:file:silent relax.out -out:file:scorefile relax.sc > relax.log
    cd ..
done


echo "Done!"