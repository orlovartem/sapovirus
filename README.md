# Sapovirus (term project).

This repository includes files of the term project. 

## *parser_orf.py*

### Description:

Converts file with sequences in GenBank-format to fasta-format with sequence names in the following format: GenbankAccessionNumber.

The script creates 5 files with different open reading frames or untranslated regions. Each file contains one ORF (1, 2 or 1+2) or UTR (5'or 3').

Saves output-files in the directory of input-file.

### Usage:

```
parser_gb.py [-h] -input INPUT_FILE -min MIN_LENGTH -max MAX_LENGTH

  -input <str>          Path to input-file in GenBank-format
                        
  -min   <int>          Minimal length of sequence. Sequences shorter than min
                        length will not be included in the final dataset
                        (doesn't  work yet)
                        
  -max   <int>          Maximal length of sequence. Sequences longer than max
                        length will not be included in the final dataset
                        (doesn't  work yet)
```
## Requirements

* Python 3
