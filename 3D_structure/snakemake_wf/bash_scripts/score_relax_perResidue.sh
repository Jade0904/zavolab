#!/bin/bash

cd ${snakemake_params[current_working_dir]}
basename=$(head -n 1 "lowest.tag")
lowest="${basename}.pdb"
${snakemake_params[per_residue_energies_path]} -in:file:s "$lowest" -out:file:silent "${snakemake_output[0]}"
cd ${snakemake_params[working_dir]}