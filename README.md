# BioE_week_4 - Gene Finder Tool

This tool finds Open Reading Frames (ORFs) in FASTA files.

## Requirements
- Python 3.6+
- BioPython

## Installation
1. Clone this repository
2. Install dependencies: `pip install biopython`

## Usage
```console
python genefinder.py input_file.fasta
```

# Steps followed to create this repo

## Initialize Directory for Git

```bash
mkdir BioE_week_4
cd BioE_week_4
git init
touch genefinder.py README.md
```
## Implementing genefinder

```bash
nano genefinder.py
git add genefinder.py README.md
git commit -m "added genefinder.py"

# usage
python genefinder.py /home/ashhadm/genomes/e.coli.fna > output1.txt
```
## Implementing gene finder with reverse complements

```bash
touch genefinder_reverse.py
nano genefinder_reverse.py
git add genefinder_reverse.py 
git commit -m "added genefinder_reverse.py"

# usage
python genefinder_reverse.py /home/ashhadm/genomes/e.coli.fna > output2.txt
```
## Applying code to all 14 downloaded genomes

```bash
find /home/ashhadm/in_class/genomes -type f -name "*GCF*.fna" | while read genome; do python genefinder_reverse.py "$genome"; done > all_orfs.txt
```
## Implementing gene finder with length filter

```bash
touch genefinder_filtered.py
nano genefinder_filtered.py
git add genefinder_reverse.py 
git commit -m "added genefinder_filtered.py"

# usage
python genefinder_filtered.py /home/ashhadm/genomes/e.coli.fna -l 100
```
## Implementing gene finder with length, rbs site and rbs type filter

```bash
touch genefinder_rbs.py
nano genefinder_rbs.py
git add genefinder_rbs.py 
git commit -m "added genefinder_rbs.py"

# usage
python genefinder_filtered.py /home/ashhadm/genomes/e.coli.fna -l 150 -u 25 -r AGGAGG
```
## Push repo to github

```bash
git remote add origin https://github.com/ashhadm/BioE_week_4.git
git branch -M main
git push -u origin main
```
