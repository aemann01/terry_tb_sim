## Mapping

First need to index our reference database

```bash
cd 03-read_mapping/
ln ../01-simulate_aDNA/tb_genomes/tb_genomes.fa .
bwa index tb_genomes.fa
```

Map samples to full TB genome reference file

```bash
ls *fq.gz | sed 's/.fq.gz//' | while read line; do bwa aln -l 1000 -n 0.1 -t 8 tb_genomes.fa ../02-mock_community_generation/spiked_in/$line.fq.gz > $line.sai ; done
```

Convert sai to sam, sort

```bash
ls *sai | sed 's/.sai//' | while read line; do bwa samse tb_genomes.fa $line.sai $line.fq.gz > $line.sam; done
ls *sam | sed 's/.sam//' | while read line; do samtools view -bSu $line.sam | samtools sort - > $line.sort.sam; done
```

Convert to bam, filter by mapping quality score (>30)

```bash
ls *.sort.sam | sed 's/.sort.sam//' | while read line; do samtools view -bS $line.sort.sam > $line.bam; done
ls *.bam | sed 's/.bam//' | while read line; do samtools view -h -q 30 $line.bam > $line.filt.bam; done
rm *sam
```

Convert to bed file, pull stats

```bash
mkdir results
ls *filt.bam | sed 's/.filt.bam//' | parallel 'bedtools bamtobed -i {}.filt.bam > results/{}.bed'
cd results
awk -F"\t" '{print $1, "\t", $4}' 30bp.01per.bed | sed 's/:.*//' | sort | uniq -c | sed 's/^[ \t]*//g' | sed 's/ /\t/g' # do for each file
```
