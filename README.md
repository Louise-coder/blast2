**💻🧬 BLAST 2**\
_Basic Local Alignment Search Tool 2_
==============
<ins>Student</ins> : Louise LAM (21980795) M2BI

## Getting started
These instructions will help you get a functional version of my short-term Blast2 project running on your local machine. Please follow the instructions below for an optimal experience.

BLAST (Basic Local Alignment Search Tool) is a fundamental tool in bioinformatics used for performing local pairwise alignments to detect homologies between a query sequence and a database. BLAST2, also known as Gapped-BLAST, is an optimized version of this tool, offering improved performance: it is on average three times faster and has better sensitivity for detecting weak similarities. This version introduces new criteria for selecting hits to extend, aiming to reduce their number, as well as heuristics for generating alignments with gaps. 

The goal of this project was to re-implement BLAST2 for protein sequences in Python, to understand the challenges associated with reproducing programs from scientific publications and to become familiar with the practical aspects of this process.


## Installation (UNIX console)
1. Clone the project repository:
```bash
git clone https://github.com/Louise-coder/blast2.git
```

2. Navigate to the `blast2` directory:
```bash
cd blast2
```

3. Create the `conda` environment from the `environment.yml` file using the command:
```bash
conda env create -f environment.yml
```


## Usage
1. Activate the conda environment using the command:
```bash
conda activate blast2-env
```

2. Run Blast2 on the human P53 protein sequence:
```bash
# for results in the console
python src/main.py -d data/subset.fasta -q data/p53_human.fasta

# for results in a file results/out.txt
python src/main.py -d data/subset.fasta -q data/p53_human.fasta -o results/out.txt 

# to redefine the word size (k), substitution matrix (m), and e-value threshold (e)
# note: with certain parameters, execution time may be very long...
python src/main.py -d data/subset.fasta -q data/p53_human.fasta -k 4 -m pam250 -e 0.01 
```

## What about BLAST+ ?
It is also possible to perform the same search using the BLAST+ software. Since it is large and takes a long time to download, I have not included it in the `.yml` file. However, if you wish to use BLAST+ yourself, you can download it and:

1. Convert your database to binary files:
```bash
makeblastdb -in data/subset.fasta -dbtype prot -out data/db
```

2. Run `blastp`:
```bash
# le choix des paramètres a été expliqué dans le rapport
blastp -query data/p53_human.fasta -db data/db -out results/blastplus_out.txt -word_size 3 -threshold 11 -window_size 40 -evalue 0.05
```
