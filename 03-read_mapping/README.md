## Mapping

First need to index our reference database

```bash
ln ../01-simulate_aDNA/tb_genomes/tb_genomes.fa .
bwa index tb_genomes.fa
```


```bash




bwa aln -l 1000 -n 0.1 -t 8 $ref $path/Sample_"${NAME}"/"${NAME}"_trimmed_merged.fastq.gz > $path/Sample_"${NAME}"/$outdir/"${NAME}"_trimmed_merged.sai
```
