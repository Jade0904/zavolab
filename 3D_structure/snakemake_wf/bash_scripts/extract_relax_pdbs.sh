#!/bin/bash

cd ${snakemake_params[current_working_dir]}
${snakemake_params[extract_pdbs_path]} -in:file:silent ${snakemake_params[relax_out]} -in:file:tagfile "${snakemake_input[0]}"
cd ${snakemake_params[working_dir]}