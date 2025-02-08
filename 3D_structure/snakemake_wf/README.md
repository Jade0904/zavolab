```
snakemake -p --executor slurm --profile profile --use-envmodules --keep-going
```

Every input is an amino acid sequence, output is the per-residue score of the structure with lowest energy after AlphaFold+RosettaRelax.