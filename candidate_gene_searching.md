# Searching Candidate Gene
> ***Note*** : I am using local HPCC server's NCBI database to run BLASTn.
## Overview

This workflow identifies potential candidate genes within a genomic region of interest by:

1. Define your genomic interval/ know your region of interest.
2. Extracting exon coordinates from a genome annotation file.
3. Retrieving exon sequences from the reference genome.
4. Running BLAST searches against the NCBI nucleotide database.
5. Obtaining annotations for BLAST hits.
6. Prioritizing candidate genes based on biological relevance.

---

# Example Region

Suppose we found a GWAS-consistent region on chromosome Chr2A at position 13,160,000 bp and we want to identify the genes located within 3kb of this region. 

| Chromosome | Position (bp) | Position - 3kb | Position + 3kb |
|------------|---------------|------------------|------------------|
|Chr2A | 13,160,000 | 13,157,000 | 13,163,000 |

| Chromosome | Start | End | Size |
|------------|--------|------|------|
| Chr2A | 13,157,000 | 13,163,000 | 6 kb |

> ***Note***: Replace these coordinates with the genomic interval associated with your QTL, GWAS consistent region, marker interval, or region of interest.

---

# Required Input Files

## 1. Reference Genome

Genome sequence in FASTA format.

```text
genome.fa
```

---

## 2. Genome Annotation File

Genome annotation in GTF/GFF format.

```text
genome.gtf
```

---

## 3. Software Requirements

```bash
bedtools
BLAST+ (v2.16.0 or later)
```

Load modules:

```bash
module load bedtools
module load blast+/2.16.0
```

## 4. BLAST Database

NCBI nucleotide database. HPCC server I used has a database for blasting nucloetide sequence (BLASTn) called **core_nt**. 

```text
core_nt
```

Verify the NCBI database of HPCC server in your linux environment:
>***Note*** : On HPCC server I used, NCBI database is host-integrated with blast+ version 2.16.0 or above. 

```bash
module load blast+/2.16.0
blastdbcmd -db core_nt -info
```

---

---

# Workflow

## Step 1. Define the Genomic Region

Set the chromosome and coordinates of the target interval.

```bash
CHR="Chr2A"
START=13158500
END=13161500
```

---

## Step 2. Extract Exons from the Region

Extract exon coordinates located within the target interval.

```bash
awk -v OFS="\t" \
-v chr="$CHR" \
-v start="$START" \
-v end="$END" \
'$1==chr && $3=="exon" && $4>=start && $5<=end \
{print $1,$4,$5}' \
genome.gtf > region_exons.bed
```

Output:

```text
region_exons.bed
```

Example:

```text
Chr2A    13158620    13158950
Chr2A    13160175    13160490

```
### Understanding the `awk` Command ###

#### `-v OFS="\t"`

Sets the **Output Field Separator (OFS)** to a tab character (`\t`). This ensures that the output is written as **tab-separated values**, which is the standard format required for BED files.

---

#### `-v chr="$CHR" -v start="$START" -v end="$END"`

Passes the shell variables (`$CHR`, `$START`, and `$END`) into `awk` as internal variables (`chr`, `start`, and `end`), allowing them to be used within the filtering conditions.

---

#### `$1==chr && $3=="exon" && $4>=start && $5<=end`

Applies the following filters to each line of the GTF file:

| Condition | Description |
|------------|------------|
| `$1 == chr` | Chromosome name (Column 1) matches the target chromosome. |
| `$3 == "exon"` | Feature type (Column 3) is an exon. |
| `$4 >= start` | Exon start coordinate (Column 4) is greater than or equal to the region start position. |
| `$5 <= end` | Exon end coordinate (Column 5) is less than or equal to the region end position. |

Only exons that satisfy **all four conditions** are retained.

---

#### `{print $1,$4,$5}`

Prints the selected records in a **3-column BED format**:

| Output Column | Description |
|---------------|-------------|
| `$1` | Chromosome name |
| `$4` | Exon start position |
| `$5` | Exon end position |

The resulting file can be directly used as input for tools such as **BEDTools**.
---

## Step 3. Extract Exon Sequences

Retrieve nucleotide sequences corresponding to the exon coordinates.

```bash
bedtools getfasta \
-fi genome.fa \
-bed region_exons.bed \
-fo region_exons.fa
```

Output:

```text
region_exons.fa
```

Example:

```fasta
>Chr2A:1000150-1000450
ATGCGATCGATCGATCGATCGATCG...
```

---

## Step 4. Run BLAST Search

Search the extracted exon sequences against the NCBI nucleotide database.

```bash
blastn \
-query region_exons.fa \
-db core_nt \
-outfmt 6 \
-num_threads 32 \
-out exon_blast_results.tsv
```

Output:

```text
exon_blast_results.tsv
```

---

## Step 5. Extract Unique Accession IDs

Retrieve unique subject accession IDs from the BLAST output.

```bash
cut -f2 exon_blast_results.tsv \
| sort -u \
> accession_ids.tsv
```

Output:

```text
accession_ids.tsv
```

---

## Step 6. Retrieve Hit Annotations

Obtain accession titles and descriptions for the BLAST hits.

```bash
blastdbcmd \
-db core_nt \
-entry_batch accession_ids.tsv \
-outfmt "%a\t%t" \
> accession_titles.tsv
```

Output:

```text
accession_titles.tsv
```

Example:

```text
XM_123456    Auxin-responsive protein IAA17
XM_987654    ABC transporter family protein
XM_555555    Heat shock protein 70
```

---

## Step 7. Identify Candidate Genes

Review the annotations and prioritize genes based on their known biological functions.

Examples of relevant candidate gene categories:

| Trait of Interest | Candidate Gene Types |
|------------------|---------------------|
| Root growth | Auxin signaling genes |
| Shoot growth | Cytokinin signaling genes |
| Flowering time | FT, CO, FLC genes |
| Stress tolerance | HSP, WRKY, NAC genes |
| Disease resistance | NBS-LRR genes |
| Plant architecture | Gibberellin pathway genes |

---

# Complete SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name="Candidate_Gene_Search"
#SBATCH -p bigmem
#SBATCH -t 168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

module load bedtools
module load blast+/2.16.0

GFF_FILE="genome.gtf"
REF_GENOME="genome.fa"

CHR="Chr2A"
START=1000000
END=1002000

echo "Extracting exons..."

awk -v OFS="\t" \
-v chr="$CHR" \
-v start="$START" \
-v end="$END" \
'$1==chr && $3=="exon" && $4>=start && $5<=end \
{print $1,$4,$5}' \
$GFF_FILE > region_exons.bed

echo "Extracting sequences..."

bedtools getfasta \
-fi $REF_GENOME \
-bed region_exons.bed \
-fo region_exons.fa

echo "Running BLAST..."

blastn \
-query region_exons.fa \
-db core_nt \
-outfmt 6 \
-num_threads 32 \
-out exon_blast_results.tsv

echo "Retrieving accession IDs..."

cut -f2 exon_blast_results.tsv \
| sort -u \
> accession_ids.tsv

echo "Retrieving annotations..."

blastdbcmd \
-db core_nt \
-entry_batch accession_ids.tsv \
-outfmt "%a\t%t" \
> accession_titles.tsv

echo "Analysis completed."
```

---

# Input and Output Files

## Input Files

```text
genome.fa
genome.gtf
core_nt
```

## Output Files

```text
region_exons.bed
region_exons.fa
exon_blast_results.tsv
accession_ids.tsv
accession_titles.tsv
```

---

# Workflow Summary

```text
Region of Interest
        │
        ▼
Extract Exons from GTF
        │
        ▼
Retrieve Sequences
        │
        ▼
BLAST Search
        │
        ▼
Extract Accessions
        │
        ▼
Retrieve Annotations
        │
        ▼
Identify Candidate Genes
```
