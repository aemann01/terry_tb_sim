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

Prepare TB for spike in

```bash
# 100000 1%
ls *fq | parallel 'seqtk sample -s167 {} 100000 > {}.1per'
# 50000 0.5%
ls *fq | parallel 'seqtk sample -s116 {} 50000 > {}.05per'
# 10000 0.1%
ls *fq | parallel 'seqtk sample -s339 {} 10000 > {}.01per'
wc -l *per
```

Prepare mock oral communities # NOTE: running all of these in parallel at once crashes my computer (16G RAM, 8 Cores) so proceed with caution

```bash
cd /Volumes/histolytica/terry_tb_sim
# 9900000 1%
ls *fq | parallel 'seqtk sample -s320 {} 9900000 > {}.1per'
# 9950000 0.5%
ls *fq | parallel 'seqtk sample -s469 {} 9950000 > {}.05per'
 # 9990000 0.1%
ls *fq | parallel 'seqtk sample -67 {} 9990000 > {}.01per'
wc -l *per
```

Finally, concatenate appropriate sets

```bash
# 1 percent
printf "30\n50\n75\n100\n" | while read line; do cat tb_$line\bp.fq.1per /Volumes/histolytica/terry_tb_sim/mock_oral_$line\bp.fq.1per > /Volumes/histolytica/terry_tb_sim/spike/$line\bp.1per.fq ; done
# 0.5 percent
printf "30\n50\n75\n100\n" | while read line; do cat tb_$line\bp.fq.05per /Volumes/histolytica/terry_tb_sim/mock_oral_$line\bp.fq.05per > /Volumes/histolytica/terry_tb_sim/spike/$line\bp.05per.fq ; done
# 0.1 percent
printf "30\n50\n75\n100\n" | while read line; do cat tb_$line\bp.fq.01per /Volumes/histolytica/terry_tb_sim/mock_oral_$line\bp.fq.01per > /Volumes/histolytica/terry_tb_sim/spike/$line\bp.01per.fq ; done
cd /Volumes/histolytica/terry_tb_sim/spike
grep "^@" *fq -c
ls *fq | parallel 'gzip {}'
```

## TO DO: cleanup
