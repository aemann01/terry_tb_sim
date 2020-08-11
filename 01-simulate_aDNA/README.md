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
rm *fna
cat *fix > tb_genomes.fa
rm *fix
cd ..
```

## Simulation

Now we can use Gargammel to simulate aDNA fragmentation and metagenome creation.

Running Gargammel for generating aDNA from TB genomes

```bash
cp mock_oral/mock_oral.fa sim/bact/
cp mock_oral/mock_oral.fa sim/cont/
cp tb_genomes/tb_genomes.fa sim/endo/
gargammel -n 1000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/tb_sim sim
```

And finally generating aDNA from Bacterial genomes

```bash
rm sim/endo/tb_genomes.fa
cp mock_oral/mock_oral.fa sim/endo
gargammel -n 6000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_oral_sim sim
```

Now we can move onto the next step of the analysis