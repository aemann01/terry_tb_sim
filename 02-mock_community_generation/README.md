# Synthetic Dataset Processing

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
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.30bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.30bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_30bp &
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.50bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.50bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_50bp &
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.75bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.75bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_75bp &
AdapterRemoval --file1 ../../01-simulate_aDNA/samples/tb_sim_s1.100bp.fq.gz --file2 ../../01-simulate_aDNA/samples/tb_sim_s2.100bp.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename tb_100bp

# bacterial taxa
cd /Volumes/histolytica/terry_tb_sim
AdapterRemoval --file1 mock_oral_sim.30bp_s1.fq.gz --file2 mock_oral_sim.30bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_30bp &
AdapterRemoval --file1 mock_oral_sim.50bp_s1.fq.gz --file2 mock_oral_sim.50bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_50bp &
AdapterRemoval --file1 mock_oral_sim.75bp_s1.fq.gz --file2 mock_oral_sim.75bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_75bp &
AdapterRemoval --file1 mock_oral_sim.100bp_s1.fq.gz --file2 mock_oral_sim.100bp_s2.fq.gz --trimns --trimqualities --minquality 20 --gzip --collapse --basename mock_oral_100bp
```

Concatenate wanted sequences

```bash
#keeping collapsed, singletons, and high qualtiy forward reads
ls *.collapsed.gz | sed 's/.collapsed.gz//' | while read line; do cat $line.collapsed.gz $line.collapsed.truncated.gz $line.pair1.gz $line.pair1.truncated.gz $line.singletons.gz $line.singletons.truncated.gz > $line.fq.gz; done 
ls *bp.fq.gz | parallel 'gzip -d {}'
```

## Dataset Spiking

Add sequence counts to headers

```bash
# TB reads
sed 's/:.*//' tb_30bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > tb_30bp.fq.fix
sed 's/:.*//' tb_50bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > tb_50bp.fq.fix
sed 's/:.*//' tb_75bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > tb_75bp.fq.fix
sed 's/:.*//' tb_100bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > tb_100bp.fq.fix
# mock oral
sed 's/:.*//' mock_oral_30bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > mock_oral_30bp.fq.fix &
sed 's/:.*//' mock_oral_50bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > mock_oral_50bp.fq.fix &
sed 's/:.*//' mock_oral_75bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > mock_oral_75bp.fq.fix &
sed 's/:.*//' mock_oral_100bp.fq | awk '/@/{print $0=$0"_"(++i)}!/@/' > mock_oral_100bp.fq.fix
```

Prepare TB for spike in

```bash
# 200000 1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 200000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.1per.fq; done
# 100000 0.5%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 100000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.05per.fq; done
# 20000 0.1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 20000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.01per.fq; done
grep "^@" *per*fq -c 
```


Prepare mock oral communities

```bash
cd /Volumes/histolytica/terry_tb_sim
# 19800000 1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19800000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.1per.fq; done
# 19900000 0.5%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19900000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.05per.fq; done
 # 19980000 0.1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 19980000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.01per.fq; done
grep "^@" *per*fq -c 
```

Finally, concatenate appropriate sets

```bash
cd ..
# 1 percent
cat adapterremoval/tb_30bp.1per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_30bp.1per.fa > /Volumes/histolytica/terry_tb_sim/spike/30bp.1per.fa &
cat adapterremoval/tb_50bp.1per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_50bp.1per.fa > /Volumes/histolytica/terry_tb_sim/spike/50bp.1per.fa &
cat adapterremoval/tb_75bp.1per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_75bp.1per.fa > /Volumes/histolytica/terry_tb_sim/spike/75bp.1per.fa &
cat adapterremoval/tb_100bp.1per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_100bp.1per.fa > /Volumes/histolytica/terry_tb_sim/spike/100bp.1per.fa &
# 0.5 percent
cat adapterremoval/tb_30bp.05per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_30bp.05per.fa > /Volumes/histolytica/terry_tb_sim/spike/30bp.05per.fa &
cat adapterremoval/tb_50bp.05per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_50bp.05per.fa > /Volumes/histolytica/terry_tb_sim/spike/50bp.05per.fa &
cat adapterremoval/tb_75bp.05per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_75bp.05per.fa > /Volumes/histolytica/terry_tb_sim/spike/75bp.05per.fa &
cat adapterremoval/tb_100bp.05per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_100bp.05per.fa > /Volumes/histolytica/terry_tb_sim/spike/100bp.05per.fa &
# 0.1 percent
cat adapterremoval/tb_30bp.01per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_30bp.01per.fa > /Volumes/histolytica/terry_tb_sim/spike/30bp.01per.fa &
cat adapterremoval/tb_50bp.01per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_50bp.01per.fa > /Volumes/histolytica/terry_tb_sim/spike/50bp.01per.fa &
cat adapterremoval/tb_75bp.01per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_75bp.01per.fa > /Volumes/histolytica/terry_tb_sim/spike/75bp.01per.fa &
cat adapterremoval/tb_100bp.01per.fa /Volumes/histolytica/terry_tb_sim/mock_oral_100bp.01per.fa > /Volumes/histolytica/terry_tb_sim/spike/100bp.01per.fa
```



## TO DO: cleanup
