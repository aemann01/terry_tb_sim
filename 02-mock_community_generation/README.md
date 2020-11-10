# Synthetic and Real Dataset Processing

Allison E. Mann (2020)

## TO DO: relative paths

## Set Up

Create working folders

```bash
mkdir adapterremoval spiked_in
```

## Read Adapter Removal, Merging and Quality Filtering

Adapter removal, quality filter, merge reads.

```bash
cd adapterremoval
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.30bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.30bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_30bp
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.50bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.50bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_50bp
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.75bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.75bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_75bp
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.100bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.100bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_100bp

# bacterial taxa
cd /Volumes/histolytica/terry_tb_sim
AdapterRemoval --file1 mock_oral_sim.30bp_s1.fq.gz --file2 mock_oral_sim.30bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_30bp
AdapterRemoval --file1 mock_oral_sim.50bp_s1.fq.gz --file2 mock_oral_sim.50bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_50bp
AdapterRemoval --file1 mock_oral_sim.75bp_s1.fq.gz --file2 mock_oral_sim.75bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_75bp
AdapterRemoval --file1 mock_oral_sim.100bp_s1.fq.gz --file2 mock_oral_sim.100bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_100bp
```

## Dereplication

```bash
#keeping collapsed, singletons, and high qualtiy forward reads
ls *pair1*gz | parallel 'gzip -d {}' &
ls *singleton*gz | parallel 'gzip -d {}' &
ls *collapsed*gz | parallel 'gzip -d {}' 
ls | grep -v "gz" | grep -v "settings" | parallel 'fastq_to_fasta -i {} -o {}.fna.gz -Q33 -z'
ls *.collapsed.fna.gz | sed 's/.collapsed.fna.gz//' | while read line; do cat $line.collapsed.fna.gz $line.collapsed.truncated.fna.gz $line.pair1.fna.gz $line.pair1.truncated.fna.gz $line.singletons.fna.gz $line.singletons.truncated.fna.gz > $line.fa.gz; done 
rm *fna.gz *collapsed* *discarded* *pair* *singleton* *settings*
ls *fa.gz | parallel 'gzip -d {}'
```
And then remove unnecessary sequencing duplicates (not 30bp or you'll lose all of your data)

```bash
ls *fa | grep -v "30bp" | sed 's/.fa//' | parallel 'vsearch --derep_fulllength {}.fa --output ../{}.uniq.fa'
mv *30bp* ..
cd ..
rm -r adapterremoval
```

## Dataset Spiking

Add sequence counts to headers

```bash
mv mock_oral_30bp.fa mock_oral_30bp.uniq.fa
mv tb_30bp.fa tb_30bp.uniq.fa
sed 's/:.*//' mock_oral_30bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp mock_oral_30bp.uniq.fa
sed 's/:.*//' mock_oral_50bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp mock_oral_50bp.uniq.fa
sed 's/:.*//' mock_oral_75bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp mock_oral_75bp.uniq.fa
sed 's/:.*//' mock_oral_100bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp mock_oral_100bp.uniq.fa
sed 's/:.*//' tb_30bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp tb_30bp.uniq.fa
sed 's/:.*//' tb_50bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp tb_50bp.uniq.fa
sed 's/:.*//' tb_75bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp tb_75bp.uniq.fa
sed 's/:.*//' tb_100bp.uniq.fa | awk '/>/{print $0=$0"_"(++i)}!/>/' > temp
mv temp tb_100bp.uniq.fa
```

Prepare mock oral communities

```bash
# 19800000 1%
cat mock_oral_30bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19800000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_30bp.1per.fa
cat mock_oral_50bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19800000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50bp.1per.fa
 cat mock_oral_75bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19800000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_75bp.1per.fa
cat mock_oral_100bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19800000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_100bp.1per.fa
# 19900000 0.5%
cat mock_oral_30bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19900000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_30bp.05per.fa
cat mock_oral_50bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19900000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50bp.05per.fa


 cat mock_oral_75bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19900000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_75bp.05per.fa
 cat mock_oral_100bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19900000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_100bp.05per.fa
 # 19980000 0.1%
cat mock_oral_30bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19980000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_30bp.01per.fa
cat mock_oral_50bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19980000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_50bp.01per.fa
 cat mock_oral_75bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19980000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_75bp.01per.fa
 cat mock_oral_100bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19980000 | awk '{printf("%s\n%s\n",$1,$2)}' > mock_oral_100bp.01per.fa
```

Mock oral files are now generated. Should be: 19800000, 19900000, 19980000

```bash
grep ">" *oral*per*.fa -c
```

Prepare TB for spike in

```bash
# 200000 1%
cat tb_30bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 200000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_30bp.1per.fa
cat tb_50bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 200000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_50bp.1per.fa
cat tb_75bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 200000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_75bp.1per.fa
cat tb_100bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 200000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_100bp.1per.fa
# 100000 0.5%
cat tb_30bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 100000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_30bp.05per.fa
cat tb_50bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 100000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_50bp.05per.fa
cat tb_75bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 100000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_75bp.05per.fa
cat tb_100bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 100000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_100bp.05per.fa
# 20000 0.1%
cat tb_30bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 20000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_30bp.01per.fa
cat tb_50bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 20000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_50bp.01per.fa
cat tb_75bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 20000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_75bp.01per.fa
cat tb_100bp.uniq.fa | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 20000 | awk '{printf("%s\n%s\n",$1,$2)}' > tb_100bp.01per.fa
```

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
