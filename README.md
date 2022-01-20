# SIMONE (Structural varIants and Methylation using lONg rEads)

Snakemake pipeline for investigating structural variation and methylation using long-reads and tumor-control samples

## Description

SIMONE uses [minimap2](https://github.com/lh3/minimap2) to align long reads to the human reference genome. SV identification is performed using [sniffles](https://github.com/fritzsedlazeck/Sniffles) and [cuteSV](https://github.com/tjiangHIT/cuteSV). Methylation frequencies are calculated using [nanopolish](https://github.com/jts/nanopolish). Regions with differential methylation (tumor vs control) are identified by means of [pycoMeth](https://github.com/snajder-r/pycoMeth).  

## Get Data

To get the Oxford Nanopore PromethION ultra long reads: 

```bash
pip install gdown 
```
then run:

```bash
cd SIMONE
mkdir -p data
cd data
# get GRCh38 reference (hg38 with decoys)
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
#index the reference genome
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa
#get chromosomes that we need to exclude later
cut -f1 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | tail -n +25 > GRCh38_full_analysis_set_plus_decoy_hla.exclude.txt
#get chromosomes we need
cut -f1 GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | head -24 > classic.chrs.txt
# get GRCh38 gff
wget http://ftp.ensembl.org/pub/release-105/gff3/homo_sapiens/Homo_sapiens.GRCh38.105.gff3.gz
bgzip -d Homo_sapiens.GRCh38.105.gff3.gz
sed -i -e '/#/! s/^/chr/' Homo_sapiens.GRCh38.105.gff3
#get control fastq
gdown https://drive.google.com/uc?id=1fbdD4TZMSrx8LbeF0bsfwYr4WqDc8yX0
#get tumor fastq
gdown https://drive.google.com/uc?id=1tQfbXyTyxyc-XEuaKHaFKGry1RVQMcZa
#get control,tumor fast5
gdown https://drive.google.com/drive/folders/1kHrnQrGYfPFWFMpsOAOANdBUpVxeULDE --folder
cd -

```
The control and the tumor files are 1G.pass.fastq.gz (control) and 1S.pass.fastq.gz (tumor) with correponding fast5 files.

## RUN
To use nanopolish with conda env:

```bash
export HDF5_PLUGIN_PATH=/usr/local/hdf5/lib/plugin
```

then run:

### Align FASTQ to GRCh38

``` bash
snakemake align --cores 20 --use-conda 
```

### Calculate/plot coverage statistics

``` bash
snakemake depth --cores 20 --use-conda
```

### Call SVs with sniffles, cuteSV and make a consensus SVs callset selecting only tumor SVs

``` bash
snakemake svcall --cores 20 --use-conda
```

### Methylation analysis with nanopolish and pycoMeth analysis for comparing methylation status between normal and tumor

``` bash
snakemake meth --cores 20 --use-conda
```
