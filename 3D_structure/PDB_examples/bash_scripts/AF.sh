#!/bin/bash 
#SBATCH --job-name=AF
#SBATCH --time=06:00:00 
#SBATCH --mem=64G 
#SBATCH --cpus-per-task=8 
#SBATCH --partition=a100 
#SBATCH --gres=gpu:1

# this script is an example
module load AlphaFold/2.2.0
# db_preset can be changed by updating the environment variable: DB_PRESET (defaults to full_dbs)
# model_preset can be changed by updating the environment variable: MODEL_PRESET_MONOMER (Defaults to monomer)

run_alphafold.sh --fasta_paths=/scicore/home/zavolan/zhu0006/3D_structure/indel/indel_examples_pdb/chain_files/4i43_A.fasta --output_dir=/scicore/home/zavolan/zhu0006/3D_structure/indel/indel_examples_pdb/alphafold_res --max_template_date=9999-12-31 --use_gpu_relax --uniref90_database_path=$UNIREF90_DATABASE_PATH --mgnify_database_path=$MGNIFY_DATABASE_PATH --uniclust30_database_path=$UNICLUST30_DATABASE_PATH --bfd_database_path=$BFD_DATABASE_PATH --pdb70_database_path=$PDB70_DATABASE_PATH --template_mmcif_dir=$TEMPLATE_MMCIF_DIR --obsolete_pdbs_path=$OBSOLETE_PDBS_PATH --model_preset=monomer
