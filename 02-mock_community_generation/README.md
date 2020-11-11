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
ls *fq | while read line; do sed 's/:.*//' $line | awk '/@/{print $0=$0"_"(++i)}!/@/' > $line.fix 
```

Prepare TB for spike in

```bash
# 100000 1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 100000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.1per.fq; done
# 50000 0.5%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 50000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.05per.fq; done
# 10000 0.1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 10000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.01per.fq; done
grep "^@" *per*fq -c 
```

Prepare mock oral communities

```bash
cd /Volumes/histolytica/terry_tb_sim
# 9900000 1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 9900000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.1per.fq; done
# 9950000 0.5%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 9950000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.05per.fq; done
 # 9990000 0.1%
ls *fix | sed 's/.fq.fix//' | while read line; do cat $line.fq.fix | awk '/^@/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n 9950000 | awk '{printf("%s\n%s\n",$1,$2)}' | sed 's/+/\n+\n/' > $line.01per.fq; done
grep "^@" *per*fq -c 
```

Finally, concatenate appropriate sets

```bash
# 1 percent
printf "30\n50\n75\n100\n" | while read line; do cat tb_$line\bp.1per.fq /Volumes/histolytica/terry_tb_sim/spike/mock_oral_$line\bp.1per.fq > /Volumes/histolytica/terry_tb_sim/spike/$line\bp.1per.fq ; done
# 0.5 percent
printf "30\n50\n75\n100\n" | while read line; do cat tb_$line\bp.05per.fq /Volumes/histolytica/terry_tb_sim/spike/mock_oral_$line\bp.05per.fq > /Volumes/histolytica/terry_tb_sim/spike/$line\bp.05per.fq ; done
# 0.1 percent
printf "30\n50\n75\n100\n" | while read line; do cat tb_$line\bp.01per.fq /Volumes/histolytica/terry_tb_sim/spike/mock_oral_$line\bp.01per.fq > /Volumes/histolytica/terry_tb_sim/spike/$line\bp.01per.fq ; done
```

## TO DO: cleanup
