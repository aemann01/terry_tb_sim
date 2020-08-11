# Mock community generation

By Allison E. Mann

## Data Downloading and Preparation

We will download the genomes of the taxa will be putting into the synthetic
dataset.

```bash
cd mock_oral
wget -i bac.query -q
ls *gz | parallel 'gzip -d {}' 
cd ../tb_genomes
wget -i tb.query -q 
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
name=$(awk '{print $2}' ../bac.rename)
file=$(awk '{print $1}' ../bac.rename)
tagSeq
rm *fna
cat *fix > mock_oral.fa
rm *fix
cd ../mock_euks
name=$(awk '{print $2}' ../euk.rename)
file=$(awk '{print $1}' ../euk.rename)
tagSeq

#wheat has to be processed separately (12GB file)
mv GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna_fix ../sim/endo/
rm *fna
cat *fix > mock_euks.fa
rm *fix
cd ..
```

## Simulation

Now we can use Gargammel to simulate aDNA fragmentation and metagenome creation.

Running Gargammel for generating aDNA from Eukaryotic genomes

First generate reads for wheat


```bash
cp mock_oral/mock_oral.fa sim/bact/
cp mock_oral/mock_oral.fa sim/cont/
mv sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna_fix sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna
gargammel -n 1000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/wheat_sim sim
mv sim/endo/GCA_900519105.1_iwgsc_refseqv1.0_genomic.fna* mock_euks/
cp mock_euks/mock_euks.fa sim/endo/
```

Now generate reads for all other Euk. genomes

```bash
gargammel -n 5005000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_euks_sim sim
```

And finally generating aDNA from Bacterial/Archaeal genomes

```bash
rm sim/endo/mock_euks.fa
cp mock_oral/mock_oral.fa sim/endo
gargammel -n 6000000 --comp 0,0,1 --minsize 25 --maxsize 125 -damage 0.03,0.4,0.01,0.3 --loc 4.106487474 --scale 0.358874723 -o sim/mock_oral_sim sim
```

Now we can move onto the next step of the analysis, which is under 04-synthetic_real_dataset_processing