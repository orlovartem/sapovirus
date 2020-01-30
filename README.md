# Sapovirus (term project).

This repository includes files of the term project.

## *parser_orf.py*

### Description:

Converts file with sequences in GenBank-format to fasta-format with sequence names in the following format: GenbankAccessionNumber.

The script creates 5 files with different open reading frames or untranslated regions. Each file contains one ORF (1, 2 or 1+2) or UTR (5'or 3').

Saves output-files in the directory of input-file.

### Usage:

```
parser_gb.py [-h] -input INPUT_FILE

  -input <str>          Path to input-file in GenBank-format

  -r                    Remove sequences from exceptions.csv
```
## *resolve_ambiguous.py*

### Description:

Resolves ambiguous nucleotides in nucleotide sequences according to consensus in most related sequences to the region with ambiguous nt.

### Usage:

```
resolve_ambiguous.py [-h] -input INPUT_FILE -pout DIR_OUT [-w WIND_SIZE]

  -input <str>          Path to file with alignment of nt sequences in fasta-format

  -pout  <str>          Output directory. If not defined the output files will be saved in 'years' folder in the directory of input file

  -w     <str>          Window size
```
## *get_orf.ps1*

### Description:

Runs scripts (parcer_orf.py, resolve_ambiguous.py) and MAFFT.

Requires GenBank file (sapovirus_genomes.gb).

### Usage:

```
./get_orf.ps1
```
## *genotyping.py*

### Description:

The script genotypes fasta records in alignments using colored tree.

### Usage:

```
genotyping.py -in_rep REPOSITORY -in_tree TREE -in_csv CSV_TABLE

  -in_rep  <str>        Repository containing fasta alignments

  -in_tree <str>        Colored phegenetic tree in nwc format

  -in_csv  <str>        CSV table with colors and genogroups
```

## Requirements

* Python 3
* MAFFT
*
