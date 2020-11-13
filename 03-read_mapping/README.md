## Mapping

First need to index our reference database

```bash
cd 03-read_mapping/
ln ../01-simulate_aDNA/tb_genomes/tb_genomes.fa .
bwa index tb_genomes.fa
```

Map samples to full TB genome reference file

```bash
ls *fq.gz | sed 's/.fq.gz//' | while read line; do bwa aln -l 1000 -n 0.1 -t 8 tb_genomes.fa ../02-mock_community_generation/spiked_in/$line.fq.gz > $line.sai 
```

Convert sai to sam, sort

```bash
ls *sai | sed 's/.sai//' | while read line; do bwa samse tb_genomes.fa $line.sai $line.fq.gz > $line.sam; done
ls *sam | sed 's/.sam//' | while read line; do samtools view -bSu $line.sam | samtools sort - $line.sort.sam; done
```
