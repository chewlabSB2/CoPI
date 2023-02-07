# Counting Peptide Inserts (CoPI)

Counting Peptide Inserts (CoPI) is a simple text-based search and analysis program which does quality filtering and counts the occurrences of peptide inserts in a NGS dataset for any protein engineering (eg. Directed Evolution). CoPI is written in Python 3.0 and executable on Windows or Linux platforms.

Installation 
------------

Installation:
```bash 
git clone https://github.com/chewlabSB2/CoPI.git
cd CoPI 
python setup.py install
```

Simple Usage
------------

```bash
CoPI --parameters sample.config -d
```

Config File
-----------

A sample of the config file could be found [here](https://github.com/chewlabSB2/CoPI/blob/main/sample.config).

Brief Explaination of each parameters:

| Key        | Values          |
|-----------:|:---------------:|
| reference | Nucleotide Sequence |
| anchor1 | Flanking Sequences (9bp) |
| anchor2 | Flanking Sequences (9bp) |
| read1 | Paired-Ends Read1 (Fastq) |
| read2 | Paired-Ends Read2 (Fastq) |
| unpaired | Unpaired Reads (List of Fastq) |
| paired | Read1 and Read2 are paired ends (Default: True) |
| merge | Option to merge paired ends (Default: True) |
| prefix | Prefix of Output |
| threads | Number of CPU Cores |
| threshold | Lower bound of counts not to collapse (Default: 50) |
| collapse | Collpase low count reads (Default: True) |
| plot | Plot insertion and in-frame Amino Acid Distribution (Default: True) |
| label1 | Label for in-frame Amino Acid Distribution plot (Default: Label1) |
| label2 | Label for in-frame Amino Acid Distribution plot (Default: Label2) |

When choosing flanking anchors, it is recommended to leave 6bp (2AAs) between the anchors. 

If you have two unpaired set of reads, create a new line with unpaired=<>

Expected Output
---------------
1. Summary of Analysis
2. Pre-Collapse Counts
3. Pre-Collapse Counts (Inclusive of Stop Codon)
4. Post-Collapse Counts 
5. Post-Collapse Counts (Inclusive of Stop Codon)
6. Plot of the Insertion Distribution
7. Plot of the in-frame Amino Acid Distribtuion
8. Fastp output (If Fastp is available)

Issues Faced
------------
Keyword Error 48
Overloaded RAM
