# Mock community generation

Allison E. Mann (2020)

## Environment Setup

```bash
mkdir 01-simulate_aDNA/sim && mkdir 01-simulate_aDNA/sim/endo && mkdir 01-simulate_aDNA/sim/cont && mkdir 01-simulate_aDNA/sim/bact
```

## Data Downloading and Preparation

We will download the genomes of the taxa will be putting into the synthetic
dataset.

```bash
cd 01-simulate_aDNA/mock_oral
wget -i bac.query -q
ls *gz | parallel 'gzip -d {}' 
cd ../tb_genomes
wget -i tb.query -q 
ls *gz | parallel 'gzip -d {}' 
# only do below if we are also pulling the sra data for reference
# cat tb.sra.query | xargs -I{} prefetch {}
cd ..
```

Now we can tag each contig by what genome it is derived from

```bash
tagSeq()
{
    set $file
    for i in $name; do
        sed "s/^>/>${i}_/" ${1} > ${1}_fix
        shift
    done
}
cd mock_oral 
name=$(awk '{print $2}' bac.rename)
file=$(awk '{print $1}' bac.rename)
tagSeq
rm *fna
cat *fix > mock_oral.fa
rm *fix
cd ../tb_genomes
name=$(awk '{print $2}' tb.rename)
file=$(awk '{print $1}' tb.rename)
tagSeq
cat *fix TBancestor.fasta > tb_genomes.fa
rm *fna
rm *fix
cd ..
```

## Simulation

Now we can use Gargammel to simulate aDNA fragmentation and metagenome creation.

Running Gargammel for generating aDNA from TB ancestral genome

Generate simulated data 30bp

```bash
cp mock_oral/mock_oral.fa sim/bact/
cp mock_oral/mock_oral.fa sim/cont/
cp tb_genomes/TBancestor.fasta sim/endo/
gargammel -n 1000000 --comp 0,0,1 -l 30 -damage 0.03,0.4,0.01,0.3 -o sim/tb_sim sim
mkdir samples
mv sim/tb_s*fq.gz samples
mv samples/tb_sim_s1.fq.gz samples/tb_sim_s1.30bp.fq.gz
mv samples/tb_sim_s2.fq.gz samples/tb_sim_s2.30bp.fq.gz
```

50bp

```bash
gargammel -n 1000000 --comp 0,0,1 -l 50 -damage 0.03,0.4,0.01,0.3 -o sim/tb_sim sim
mv sim/tb_s*fq.gz samples
mv samples/tb_sim_s1.fq.gz samples/tb_sim_s1.50bp.fq.gz
mv samples/tb_sim_s2.fq.gz samples/tb_sim_s2.50bp.fq.gz
```

75bp

```bash
gargammel -n 1000000 --comp 0,0,1 -l 75 -damage 0.03,0.4,0.01,0.3 -o sim/tb_sim sim
mv sim/tb_s*fq.gz samples
mv samples/tb_sim_s1.fq.gz samples/tb_sim_s1.75bp.fq.gz
mv samples/tb_sim_s2.fq.gz samples/tb_sim_s2.75bp.fq.gz
```

100bp

```bash
gargammel -n 1000000 --comp 0,0,1 -l 100 -damage 0.03,0.4,0.01,0.3 -o sim/tb_sim sim
mv sim/tb_s*fq.gz samples
mv samples/tb_sim_s1.fq.gz samples/tb_sim_s1.100bp.fq.gz
mv samples/tb_sim_s2.fq.gz samples/tb_sim_s2.100bp.fq.gz
```

Do the same with the bacterial genomes (@ different read lengths). These files are huge so save (-o) to a disk with enough space (here /Volumes/histolytica/terry_tb_sim).

```bash
rm sim/endo/TBancestor.fasta
cp mock_oral/mock_oral.fa sim/endo
gargammel -n 25000000 --comp 0,0,1 -l 30 -damage 0.03,0.4,0.01,0.3 -o /Volumes/histolytica/terry_tb_sim/mock_oral_sim.30bp sim &
gargammel -n 25000000 --comp 0,0,1 -l 50 -damage 0.03,0.4,0.01,0.3 -o /Volumes/histolytica/terry_tb_sim/mock_oral_sim.50bp sim &
gargammel -n 25000000 --comp 0,0,1 -l 75 -damage 0.03,0.4,0.01,0.3 -o /Volumes/histolytica/terry_tb_sim/mock_oral_sim.75bp sim &
gargammel -n 25000000 --comp 0,0,1 -l 100 -damage 0.03,0.4,0.01,0.3 -o /Volumes/histolytica/terry_tb_sim/mock_oral_sim.100bp sim &

```

Now we can move onto the next step of the analysis: 02-mock_community_generation
