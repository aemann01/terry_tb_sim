# Synthetic and Real Dataset Processing

Allison E. Mann (2020)

## Set Up

Create working folders

```bash
mkdir adapterremoval spiked_in
```

## Read Adapter Removal, Merging and Quality Filtering

Adapter removal, quality filter, merge reads.

```bash
cd adapterremoval
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/tb_sim_30bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/tb_sim_30bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_30bp
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/tb_sim_50bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/tb_sim_50bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_50bp
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/tb_sim_75bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/tb_sim_75bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_75bp
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/tb_sim_100bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/tb_sim_100bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_100bp

# bacterial taxa
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/mock_oral_30bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/mock_oral_sim_30bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_30bp
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/mock_oral_sim_50bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/mock_oral_sim_50bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_50bp
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/mock_oral_sim_75bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/mock_oral_sim_75bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_75bp
AdapterRemoval --file1 ../../01-simulate_aDNA/sim/mock_oral_sim_100bp_s1.fq.gz --file2 ../../01-simulate_aDNA/sim/mock_oral_sim_100bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_100bp
```

## Dereplication

```bash
#keeping collapsed, singletons, and high qualtiy forward reads
ls *pair1*gz | parallel 'gzip -d {}' &
ls *singleton*gz | parallel 'gzip -d {}' &
ls *collapsed*gz | parallel 'gzip -d {}' 
ls | grep -v "gz" | grep -v "settings" | parallel 'fastq_to_fasta -i {} -o {}.fna -Q33 -z'
ls *.collapsed.fna.gz | sed 's/.collapsed.fna.gz//' | while read line; do cat $line.collapsed.fna.gz $line.collapsed.truncated.fna.gz $line.pair1.fna.gz $line.pair1.truncated.fna.gz $line.singletons.fna.gz $line.singletons.truncated.fna.gz > $line.fa.gz; done 
rm *fna
ls *fa.gz | parallel 'gzip -d {}'
```
And then remove unnecessary sequencing duplicates (not 30bp or you'll lose all of your data)

```bash
ls *fa | gzip -v "309bp" | sed 's/.fa//' | parallel 'vsearch --derep_fulllength {}.fa --output ../{}.uniq.fa'
mv *30bp* ..
cd ..
rm -r adapterremoval
```

## Dataset Spiking

Add sequence counts to headers

```bash
mv mock_oral_30bp.fa mock_oral_30bp.uniq.fa

sed 's/:.*//' mock_euks.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_euks.trim.fa
sed 's/:.*//' mock_oral.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_oral.trim.fa
sed 's/:.*//' mock_wheat.trim.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp.mock
mv temp.mock mock_wheat.trim.fa
```

Next split out by species

```bash
cat mock_wheat.trim.fa mock_euks.trim.fa > temp.seqs
mv temp.seqs mock_euks.trim.fa
awk '{print $2}' euk.rename | parallel 'grep {} -A 1 mock_euks.trim.fa > temp/{}.fa'
echo "How many dietary reads post quality filter and dereplication?"
grep ">" temp/*fa -c
```

Now we can make the mock oral commmunities

```bash
#5k depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4995000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_5k.fa
#500 depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999500 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_500.fa
#50 depth
cat mock_oral.trim.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 4999950 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50.fa
```

Mock oral files are now generated. Should be: 4995000, 499500, 4999950

```bash
grep ">" *oral*5*.fa -c
```

Mock Eukaryotic community preparation

```bash
cd temp
ls *fa | sed 's/.fa$//' > diet.ids
#5000 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 5000 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.5k.fa; done &
#500 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 500 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.500.fa; done &
#50 subsample
cat diet.ids | while read line; do cat $line.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 50 | awk '{printf("%s\n%s\n",$1,$2)}' > $line.50.fa; done

Finally, concatenate appropriate sets, and assign sample names

```bash
cat diet.ids | while read line; do cat $line.5k.fa ../mock_oral_5k.fa > ../samples/$line.5k.spike.fa; done &
cat diet.ids | while read line; do cat $line.500.fa ../mock_oral_500.fa > ../samples/$line.500.spike.fa; done &
cat diet.ids | while read line; do cat $line.50.fa ../mock_oral_50.fa > ../samples/$line.50.spike.fa; done

```bash
echo "Generate mock samples, should all be 5000000"
grep ">" samples/*spike*fa -c
```
And cleanup

```
rm -r temp
mv neanderthal.trim.fa *calc*trim* modHuman.trim.fa samples/
```
