# Searching for Candidate Gene in a Genomic Region
### NOTE: I am using local HPCC server's NCBI database to run BLASTn.
## Overview

This workflow identifies potential candidate genes within a genomic region of interest by:

1. Defining a genomic interval.
2. Extracting exon coordinates from a genome annotation file.
3. Retrieving exon sequences from the reference genome.
4. Running BLAST searches against the NCBI nucleotide database.
5. Obtaining annotations for BLAST hits.
6. Prioritizing candidate genes based on biological relevance.

---

# Example Region

For this example, we will investigate a **2 kb region on chromosome 2A**.

| Chromosome | Start | End | Size |
|------------|--------|------|------|
| Chr2A | 1,000,000 | 1,002,000 | 2 kb |

Replace these coordinates with the genomic interval associated with your QTL, GWAS consistent region, marker interval, or region of interest.

---

# Required Input Files

## 1. Reference Genome

Genome sequence in FASTA format.

```text
genome.fa
```

---

## 2. Genome Annotation File

Genome annotation in GTF format.

```text
genome.gtf
```

---

## 3. BLAST Database

NCBI nucleotide database (or any custom nucleotide database).

```text
core_nt
```

Verify the database:

```bash
blastdbcmd -db core_nt -info
```

---

# Software Requirements

```bash
bedtools
BLAST+ (v2.16.0 or later)
awk
```

Load modules:

```bash
module load bedtools
module load blast+/2.16.0
```

---

# Workflow

## Step 1. Define the Genomic Region

Set the chromosome and coordinates of the target interval.

```bash
CHR="Chr2A"
START=1000000
END=1002000
```

Calculate interval size:

```bash
echo $((END-START))
```

Expected output:

```text
2000
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
Chr2A    1000150    1000450
Chr2A    1000800    1001100
```

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
