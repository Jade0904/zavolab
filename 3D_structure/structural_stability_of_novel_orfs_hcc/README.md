**extract_fasta.ipynb:**

Extract sequences from .xlsx file and save in fasta files.

**pairs_csv.py:**

Integrate csv file for per-residue metrics in each protein pair (readthrough vs. canonical).

*possible options*
- --sl: support level (Proteomic\_Evidence, RPF\_Evidence, No\_Evidence)
- --pid: Protein ID, [ProteinName]\_[Index]

*use example*

```
python pairs_csv.py --sl Proteomics_Evidence --pid ARHGEF12_0
```

**generate_array_cmd.ipynb:**

Generate cmd batch file to use pairs\_csv.py.

**generate_wfconfig.ipynb:**

Generate "config.yaml" for snakemake workflow.

**stability_metrics.ipynb:**

Generate OVERALL metrics for each protein pair, one row per protein ID. Per-residue metrics were extracted instead by pairs\_csv.py and generate\_array\_cmd.ipynb using Slurm array.

**plots.ipynb:**

Generate plots from stability metrics, for both of overall structures and of per-residue scores.
