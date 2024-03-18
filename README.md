IRworkflow.py is used to extract the count of reads that supports (1) intron retention events (2) splicing at spliced sites.

```
python IRworkflow.py --a <gtf file> --sj <SJ.out.tab> --bam <bam file> --out <output folder>
```

In which:

(1) gtf file: annotation file in gtf format.

(2) SJ.out.tab & bam file: output from STAR.

(3) output folder: for example "xx/xx/out", without "/" at the end.

This will result in six files in the output folder:

(1) All the spliced sites from both gtf file and SJ.out.tab (merged.bed);

(2) The filtered bed file from bam file (filteredReads.bed);

(3) The filtered bed flle "groupby" the position of the reads, using the "name" column to save the count of the same reads (filteredNameAsCount.bed);

(4,5) The output of "bedtools intersect" (intersect.bed & intersect.log);

(6) The result table containing the count of the reads of both events (result.csv).
