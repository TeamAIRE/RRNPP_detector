# RRNPP_detector

## 1. What does ```RRNPP_detector``` do?

```RRNPP_detector``` aims at identifying candidate RRNPP-type quorum sensing systems in chromosomes or in mobile genetic elements (plasmids, phages) of gram-positive bacteria (e.g. Firmicutes).

```RRNPP_detector``` defines candidate RRNPP-type quorum sensing systems as tandems of adjacent ORFs encoding a candidate receptor (250-500aa protein matching HMMs of peptide-binding tetraticopeptide repeats (TPRs)) and a candidate pro-peptide (10-70aa protein predicted to be excreted via the SEC-translocon) 

```RRNPP_detector``` outputs three types of files:  
\- ```summary.tsv```: a tabular summary of the candidate RRNPP-type quorum sensing systems detected in the target genome(s)  
\- ```candidate_receptors.faa```: a fasta of the protein sequences of the candidate receptors  
\- ```candidate_propeptides.faa```: a fasta of the protein sequences of the candidate propeptides

As a few experimentally-validated RRNPP-type propeptides are known to be excreted via other secretion systems than the SEC-translocon (e.g. PrgQ and Shp), ```RRNPP_detector``` was designed to specify two kinds of prefixes for each output files:  
\- ```conservative_```: satisfying the condition that the propeptides must be predicted to be secreted via the SEC-translocon  
\- ```exploratory_```: without the above condition

## 2. How to install ```RRNPP_detector```?

### Clone the repository

```bash
cd ~
git clone https://gitlab.com/charles-bernard/rrnpp_detector.git
```

### Manual installation of the dependencies

* ```pandas```: a python library to manipulate dataframes

```bash
pip3 install pandas
```

* ```prodigal```: a tool to detect ORFs in prokaryotic genomes

```bash
sudo apt-get install prodigal
```

* ```HMMER3```: a toolbox to manipulate HMMs

```bash
sudo apt-get install hmmer
```

* ```BLAST+```: a toolbox to run sequence homology searches

```bash
sudo apt-get install ncbi-blast+
```

* ```Signalp_5.0```: a tool to predict the presence of signal peptides in proteins

```RRNPP_detector``` comes by default with the binary of ```SignalP version 5.0b Linux x86_64```. However, if you want to use another version of SignalP, the software is freely available for academic users at: https://services.healthtech.dtu.dk/software.php

## 3. How to use ```RRNPP_detector```?

### General usage

```RRNPP_detector``` must be run with python3 against one or multiple (meta)genome(s)
```
usage: rrnpp_detector.py [-h] [-o OUT_DIR] [--fna FNA] [--faa FAA] [--gff GFF]
                         [--ft FEATURE_TBL] [--cpu CPU] [--preserve_ram]

RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-
positive bacteria and bacteriophages genomes

optional arguments:
  -h, --help        show this help message and exit
  -o OUT_DIR        path to the output directory (default is current
                    directory)
  --fna FNA         path to the fasta of the target genome(s) (will run
                    Prodigal to detect ORFs)
  --faa FAA         path to the fasta of the protein sequences of the target
                    genome(s) (requires additional --gff or --ft option)
  --gff GFF         path to the annotations of the target genome(s) in gff
  --ft FEATURE_TBL  path to the annotations of the target genome(s) in the
                    NCBI_assembly feature_table format
  --cpu CPU         number of cpu to use (default is 1)
  --preserve_ram    minimize RAM usage at the expense of speed (will process
                    target genomes one by one instead of all together)
```

### Option 1: Provide one or multiple genome(s)

Ìf you provide one or multiple genome(s) (stored in one ```fna```), ```RRNPP_detector``` will run ```prodigal``` to predict the ORFs of each genome and produce an annotation file (```gff3```) and a fasta of the predicted protein sequences (```faa```)

```bash
python rrnpp_detector.py --fna ~/bacillus_subtilis.fna
```

### Option 2: Provide all proteins from one or multiple genome(s) + the annotations of the genome(s)

If you provide all proteins predicted from one or multiple genome(s) (stored in one ```faa```) and the annotations of the genomes (concatenated in one ```gff3``` or one NCBI_assembly ```feature_table```), ```RRNPP_detector``` will integrate the pre-computed annotations. Option 2 has the advantage to use the reference names/ids of the annotated proteins of the target genome(s) along the analysis. 

```bash
# with a gff
python rrnpp_detector.py --faa ~/bacillus_subtilis.faa --gff ~/bacillus_subtilis.gff
# with a feature table
python rrnpp_detector.py --faa ~/bacillus_subtilis.faa --ft ~/bacillus_subtilis_feature_table.txt
`````` 

## 4. Practical example of analysis

In this example, we will propose to use ```RRNPP_detector``` against all genomes of Viruses available on the NCBI.

To this end, we will download the protein sequences and the annotated features of these genomes from the NCBI assembly database: https://www.ncbi.nlm.nih.gov/assembly/. Once on the web page, enter the following query search:

```
Viruses[ORGN] AND "Latest GenBank"[Filter]
```

Then, we will click the button "Download Assemblies" and download the two following files:

* "Protein FASTA (.faa)" with "Genbank" selected as source database
* "Feature table (.txt)" with "Genbank" selected as source database

Now we will move these two archives in a dedicated directory

```bash
mkdir ~/genomes_of_viruses
mv ~/Downloads/genome_assemblies_prot_fasta.tar ~/genomes_of_viruses
mv ~/Download/genome_assemblies_features.tar ~/genomes_of_viruses
```

Next, we will extract the files from these two archives and concatenate all protein fastas in one ```faa``` and all feature tables in one ```feature_table```

```bash
cd ~/genomes_of_viruses
for FILE in *.tar; do tar -xvf "$FILE"; done
cd ncbi-genomes*
for FILE in *.gz; do gzip -d "$FILE"; done
for FILE in *_feature_table.txt; do cat "$FILE" >> ../viruses_feature_table.txt; done
for FILE in *_protein.faa; do cat "$FILE" >> ../viruses_protein.faa; done
```

All that remains to do is to execute ```RRNPP_detector``` against this dataset. For instance:

```bash
cd ~/Programs/rrnpp_detector
python rrnpp_detector.py --faa ~/genomes_of_viruses/viruses_protein.faa --ft ~/genomes_of_viruses/viruses_feature_table.txt -o ~/genomes_of_viruses --cpu 20
```
