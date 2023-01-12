# RRNPP_detector

## 1. What does ```RRNPP_detector``` do?

```RRNPP_detector``` aims at identifying known and novel RRNPP-type quorum sensing systems in chromosomes, mobile genetic elements and bacteriophages of Firmicutes.

```RRNPP_detector``` defines candidate RRNPP-type quorum sensing systems as tandems of adjacent ORFs encoding a candidate receptor (250-500aa protein matching HMMs of peptide-binding tetraticopeptide repeats (TPRs)) and a candidate pro-peptide (10-100aa protein predicted to be excreted via the SEC-translocon or matching the amino-acid profile of SHPs) 

```RRNPP_detector``` outputs candidate RRNPP-type systems with 3 different detection strictness:  
* ```strict```: All receptor-propeptide pairs oriented in a divergent or co-directional context with the propeptide downstream from the receptor and matching the HMM profile of SHPs or predicted to undergo a SEC/SPI-dependent secretion according to ```SignalP```
* ```relaxed```: All remaining receptor-propeptide pairs in which the propeptide harbors any of the SP(Sec/SPI), TAT(Tat/SPI) or LIPO(Sec/SPII) signal sequence according to ```PrediSi``` or ```SignalP```    
* ```loose```: All remaining receptor-propeptide pairs in which the propeptide is preceded by a RBS motif with both a high score (according to ```Prodigal```) and a high usage across prokaryotes (according to Omotajo et al.)


## 2. How to install ```RRNPP_detector```?

### Create a dedicated and isolated conda environment for ```RRNPP_detector``` (recommended but not mandatory)

```bash
conda create --name rrnpp_detector
conda activate rrnpp_detector
````

### Install dependencies 

Using conda

```bash
conda install -c bioconda pandas orfipy prodigal hmmer blast openjdk

```


```orfipy``` needs ```sqlite3``` to be compiled. If the installation exits with the "sqlite3.h: No such file or directory" error message, please try the following:

```bash
sudo apt-get install libsqlite3-dev
conda install -c bioconda pandas orfipy prodigal hmmer blast openjdk
```

Alternative to installation via conda:

```bash
sudo apt-get install prodigal hmmer ncbi-blast+
pip3 install pandas orfipy
```

### Clone the repository

```bash
cd ~
git clone https://github.com/TeamAIRE/RRNPP_detector.git
cd RRNPP_detector/
````

### Uncompress SignalP and compile PrediSi


```bash
tar -xvzf signalp-5.0b.tar.gz

# optional: if you intend to use Predisi in addition to SignalP
cd predisi
javac JSPP.java
cd ../
```

```RRNPP_detector``` comes by default with the binary of ```SignalP version 5.0b Linux x86_64```. However, if you want to use another version of SignalP, the software is freely available for academic users at: https://services.healthtech.dtu.dk/software.php

## 3. How to use ```RRNPP_detector```?

### General usage

```RRNPP_detector``` must be run with python3 against one or multiple contig(s)/(meta)genome(s)
```
usage: rrnpp_detector.py [-h] [--version] [-o OUT_DIR] [--fna FNA] [--faa FAA] [--gff GFF] [--ft FEATURE_TBL] [--cpu CPU] [--preserve_ram] [--keep_working_dir] [--min_pl MIN_PROPEPTIDE_LEN]
                         [--max_pl MAX_PROPEPTIDE_LEN] [--min_rl MIN_RECEPTOR_LEN] [--max_rl MAX_RECEPTOR_LEN] [--min_igd MIN_INTERGEN_DIST] [--max_igd MAX_INTERGEN_DIST] [--expand_to_homologs]
                         [--tprpred] [--predisi]

RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in chromosomes, plasmids and bacteriophages of Firmicutes

Main arguments
 -h, --help               show this help message and exit
  --version               print version number and exit.
  -o OUT_DIR              path to output directory (default is current directory)
  --fna FNA               path to the fasta of the target genome(s) (will run Prodigal to detect CDSs if faa not provided)
  --faa FAA               path to fasta of the protein sequences of the target genome(s) (requires additional --gff or --ft option)
  --gff GFF               path to the annotations of the target genome(s) in gff
  --ft FEATURE_TBL        path to the annotations of the target genome(s) in the NCBI_assembly feature_table format
  --cpu CPU               number of cpu to use (default is 1)
  --chunk_size CHUNK_SIZE nb target genomes to be processed altogether to preserve RAM usage 
                          (e.g. if --chunk_size 100, then the program will process 100 genomes by 100 genomes instead of all together)
  --keep_working_dir      keep the directory of intermediate files

Specific optional parameters for searching RRNPP systems
  --min_pl  MIN_PROPEPTIDE_LEN  minimal propeptide length (default=10)
  --max_pl  MAX_PROPEPTIDE_LEN  maximal propeptide length (default=70)
  --min_rl  MIN_RECEPTOR_LEN    minimal receptor length (default=250)
  --max_rl  MAX_RECEPTOR_LEN    maximal receptor length (default=500)
  --min_igd MIN_INTERGEN_DIST   minimal intergenic distance (default=-20)
  --max_igd MAX_INTERGEN_DIST   maximal intergenic distance (default=400)
  --start_codons START_CODONS   comma-separated list of start codons to consider for ORF calling (default=ATG)
  --rbs_bins RBS_BINS           comma-separated list of Prodigal's RBS bins to consider for ORF calling (default=27,24,23,22,20,19,16,15,14,13,12,6), 
                                to by bypass the filter, use --rbs_bins 27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0
  --expand_to_homologs          use detected systems as seeds to detect putative homologous systems missed by RRNPP_detector
  --tprpred                     run tprpred in addition to hmmsearch for TPR motifs detection
  --predisi                     run predisi in addition to signalp for detection of propeptides with a signal sequence (warning: this increases the risk of false positives)

```

### Best option: Provide one or multiple genome(s) (```fna```) along with annotation(s) (```gff```) and proteome(s) (```faa```)

When ```rrnpp_detector``` is fed with one or multiple genome(s) (concatenated in one ```fna``` file) along with their corresponding annotated proteome(s) (concatenated in one ```faa``` file) and annotations (concatenated in one ```gff``` or one NCBI_assembly ```feature_table file```), both annotated proteins and unannotated proteins (searched in the genomic vicinity of annotated receptors) will be present in the output.

```bash
python rrnpp_detector.py --fna bacillus_subtilis/genome.fna --faa bacillus_subtilis/proteome.faa --gff bacillus_subtilis/annotations.gff
```

### Alternative 1: Provide one or multiple genome(s) (```fna```)

If you provide only a nucleotide fasta file (```fna```) to rrnpp_detector, the tool won't work with already annotated proteins. Instead, the proteins will be detected with ```prodigal```. Search for unannotated small ORFs in the vicinity of receptors will be performed.

```bash
python rrnpp_detector.py --fna bacillus_subtilis/genome.fna
```

### Alternative 2: Provide all proteins from one or multiple genome(s) (```faa```) + the annotations of the genome(s) (```gff```)

If you don't provide a nucleotide fasta file but only a proteome file (```faa```) and an annotation file, ```RRNPP_detector``` will work only with annotated proteins. The downside of this option is that it won't search for unannotated small ORFs in the vicinity of receptors

```bash
# with a gff
python rrnpp_detector.py --faa bacillus_subtilis/proteome.faa --gff bacillus_subtilis/annotations.gff
# with a feature table
python rrnpp_detector.py --faa bacillus_subtilis/proteome.faa --ft bacillus_subtilis/feature_table.txt
`````` 

## 4. Custom search

If you wish to design a custom search with very specific parameters, you can change the parameters hard coded in ```rrnpp_detector/preprocessing.py```

## 5. Practical example of analysis

In this example, we will propose to use ```RRNPP_detector``` against all genomes of Viruses available on the NCBI.

First, we will use the ```ncbi_datasets``` command line tool to fetch the genomes from the NCBI assembly database.
Here is a tutorial on how to install ```ncbi_datasets```: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

```bash
conda activate ncbi_datasets

DOWNLOAD_DIRECTORY="/home/viral_genomes";
mkdir "$DOWNLOAD_DIRECTORY";

# will donwload genomes (fna), protein sequences (faa) and annotations (gff3) of complete genomes of Viruses available on Genbank
datasets download genome taxon "Viruses" \
--annotated \
--assembly-level complete \
--assembly-source genbank \
--include genome,protein,gff3 \
--filename "$DOWNLOAD_DIRECTORY"/assemblies.zip

```

Next, we will extract the files from this archive and concatenate all genomes in one ```genomes.fna``` file, all protein fastas in one ```proteomes.faa``` file and all gffs in one ```annotations.gff``` file

```bash
cd "$DOWNLOAD_DIRECTORY"
unzip assemblies.zip
cd ncbi_dataset/data
for ASSEMBLY in GCA_*; do 
  cat "$ASSEMBLY"/*.fna >> ../../viral_genomes.fna; rm "$ASSEMBLY"/*.fna;
  cat "$ASSEMBLY"/*.faa >> ../../viral_proteomes.faa; rm "$ASSEMBLY"/*.faa;
  cat "$ASSEMBLY"/*.gff >> ../../viral_annotations.gff; rm "$ASSEMBLY"/*.gff;
done
```

All that remains to do is to execute ```RRNPP_detector``` against this dataset. For instance:

```bash
cd ~/Programs/rrnpp_detector
python rrnpp_detector.py --fna "$DOWNLOAD_DIRECTORY"/viral_genomes.fna --faa "$DOWNLOAD_DIRECTORY"/viral_proteomes.faa --gff "$DOWNLOAD_DIRECTORY"/viral_annotations.gff -o ~/rrnpp_detector_vs_viruses --cpu 20
```

## 6. Help improving RRNPP_detector

If you have suggestions to improve the tool or would like to report bugs, please post your message on the Issues section of this repository.

## 7. Historic of versions

### v1.1.0
* Receptor detection:
  - A filter of >65% coverage of the query HMM has been introduced to minimize false positives
  - The additional Com_TPR (PF18710) HMM has been included within the library of HMMs of TPRs

* RAM Usage:
  - the ```--chunk_size``` option has been introduced, which enables to divide the target dataset into chunks of N genomes in an effort to preserve RAM usage. 

* Search options:
  - the ```--rbs_bins``` and ```start_codons``` options have been introduced, which enable to explictily specify the RBS motifs and start codons to consider for the detection of small ORFs encoded in the vicinity of candidate receptors


### v1.0.0
* Receptor detection: 
  - The library of HMMs for TPRs has been updated using the HMMs present in ```interproscan-5.56-89.0```
  - 20 additional HMMs of TPRs from Pfam have been included: TPR_1 (PF00515.30), TPR_2 (PF07719.19), TPR_3 (PF07720.14), TPR_4 (PF07721.16), TPR_5 (PF12688.9), TPR_6 (PF13174.8), TPR_7 (PF13176.8), TPR_9 (PF13371.8), TPR_10 (PF13374.8), TPR_11 (PF13414.8), TPR_14 (PF13428.8), TPR_15 (PF13429.8), TPR_16 (PF13432.8), TPR_17 (PF13431.8), TPR_18 (PF13512.8), TPR_19 (PF14559.8), TPR_20 (PF14561.8), TPR_21 (PF09976.11), TPR_22 (PF18833.3) and TPR_MalT (PF17874.3)
  - The ```tprpred``` software has been integrated and can be called in complement of ```hmmsearch``` to increase the sensitivity of the tool
  - An HMM of the AimR family has been built and is now used to identify these receptors with more sensitivity

* Propeptide detection:
  - The main improvement of ```v1.0.0``` lies in a new algorithm to detect small peptides encoded in the vicinity of receptors, as such small peptides are typically absent from annotation files. We called this method ```SPRAT``` for *Small Peptides with RBS Annotation Tool*. This method identifies peptides preceded by a Shine-Dalgarno RBS in the flanking regions of each receptor, using the 27 hierarchical regular expressions introduced by Prodigal to detect SD-RBS motifs. This is justified by the fact that 90% of the canonical genes encoded by Firmicutes have an SD-RBS upstream. Optionally, the user can submit a list of possible start codons to consider for the detection of putative small pepite-coding ORFs (by default, only ATG is considered). 
  - The ```PrediSi``` software has been integrated an can be called in complement of ```SignalP``` to increase the sensitivity of the tool
  - An HMM of SHP propeptides has been built and is now used to identify SHP propeptides since SHPs are not exported via the SEC-translocon and are therefore not returned by ```SignalP``` or ```PrediSi```

* Iterative Search:
  - When the target database is large, it may be relevant to use the detected systems as baits to fish homologous systems that did not pass the conservative thresholds of ```RRNPP_detector```. If this option increases the sensitivity of the tool, it also increases the risk of false positives


### v0.0.1

Initial push of ```RRNPP_detector```
