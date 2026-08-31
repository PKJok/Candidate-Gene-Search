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

Suppose we found a GWAS-consistent regions on chromosome Chr2A, Chr3A, and Chr6B at position 13,160,000 bp, 2,516,000 bp, 27,209,000 bp respectively and we want to identify the genes located within 20kb of this regions. 

| Chromosome | Position (bp) | Position - 10kb (Start) | Position + 10kb (End) |
|------------|---------------|------------------|------------------|
|Chr2A | 13,160,000 | 13,150,000 | 13,170,000 |
| Chr3A | 2,516,000 | 2,506,000 | 2,526,000 |
| Chr6B | 27,209,000 | 27,199,000 | 27,219,000 |

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

We need NCBI nucleotide sequence database to blast against (blastn). HPCC server I used has a database for blasting nucloetide sequence (BLASTn) called **core_nt**. 

```text
core_nt
```

Verify the NCBI database of HPCC server in your linux environment:
>***Note*** : On HPCC server I used, NCBI database is host-integrated with blast+ version 2.16.0 or above. So, when you load the blast+/2.16.0 NCBi database is also loaded.

```bash
module load blast+/2.16.0
blastdbcmd -db core_nt -info
```

---

---

# Workflow

## Step 1. Define the Genomic Region

We formed a bash array contatining 3 string elements. Each element holds space-separated genomic coordinates. 

```bash
regions=(
  "Chr2A 13150000 13170000"
  "Chr3A 2506000 2526000 "
  "Chr6B 27199000 27219000" )
```

---
## Step 2. Running the loop to extract the sequences from all the regions
We also formed a ```for``` loop that iterates over each element in the ```regions``` array.

```bash
for region in "${regions[@]}" ; do
    read -r CHROM START END <<< "$region"
    PREFIX="${CHROM}_${START}"

    echo "$PREFIX"
    echo "$CHROM"
    echo "$START"
    echo "$END"
    ..............
    .............
 done

```
### Understanding the ```loop``` code
``` bash
read -r CHROM START END <<< "$region"
```
#### This ```read``` command with a Here-String (```<<<```) splits the ```$region``` string by spaces into three distinct variables: ```$CHROM```,```$START```,and ```$END```.

``` bash
PREFIX="${CHROM}_${START}"
````
#### Combines the chromosome name and start coordinate with an underscore (e.g.,```Chr2A_13150000```) to create unique prefix for naming file.
---
## Step 3. Extract Exons from the Region

Extract exon coordinates located within the target interval. ```$3==="exon"``` ensures we are only extracting the exons from GTF/GFF file.

```bash
awk -v OFS="\t" \
-v chrom="$CHROM" \
-v start="$START" \
-v end="$END" \
'$1==chrom && $3=="exon" && $4>=start && $5<=end {print $1,$4,$5}' \
genome.gtf > "${PREFIX}.bed"
```

Output:

```text
Chr2A_13150000.bed
Chr3A_2506000.bed
Chr6B_27199000.bed
```

Example:

```text
$more Chr2A_13150000.bed
```
```bash
Chr2A 13154242 13156423
Chr2A 13167915 13168522
Chr2A 13168927 13169155
Chr2A 13169256 13169504
Chr2A 13169597 13169647
.....
```
### Understanding the `awk` Command ###

#### `-v OFS="\t"`

Sets the **Output Field Separator (OFS)** to a tab character (`\t`). This ensures that the output is written as **tab-separated values**, which is the standard format required for BED files.

---

#### `-v chrom="$CHROM" -v start="$START" -v end="$END"`

Passes the shell variables (`$CHROM`, `$START`, and `$END`) into `awk` as internal variables (`chrom`, `start`, and `end`), allowing them to be used within the filtering conditions.

---

#### `$1==chrom && $3=="exon" && $4>=start && $5<=end`

Applies the following filters to each line of the GTF/GFF file:

| Condition | Description |
|------------|------------|
| `$1 == chrom` | Chromosome name (Column 1) matches the target chromosome. |
| `$3 == "exon"` | Feature type (Column 3) is an exon. |
| `$4 >= start` | Exon start coordinate (Column 4) is greater than or equal to the region start position. |
| `$5 <= end` | Exon end coordinate (Column 5) is less than or equal to the region end position. |

Only exons that satisfy **all four conditions** are retained.

---

#### `{print $1,$4,$5}`

Prints the selected records og GTF/GFF file in a **3-column BED format**:

| Output Column | Description |
|---------------|-------------|
| `$1` | Chromosome name |
| `$4` | Exon start position |
| `$5` | Exon end position |

> The resulting file can be directly used as input for tools such as **BEDTools**.
---

## Step 3. Extract Exon Sequences

Retrieve nucleotide sequences corresponding to the exon coordinates.

```bash
bedtools getfasta \
 -fi genome.fa \
 -bed "${PREFIX}.bed" \
 -fo "${PREFIX}_exons.seqs"

```

Output:

```text
Chr2A_13150000_exons.seqs
Chr3A_2506000_exons.seqs
Chr6B_ 27199000_exons.seqs
```

Example:

```text
$more Chr2A_13150000_exons.seqs
```
```fasta
>Chr2A:13154242-13156423
GCCAAGTAGAGAGGTCCAACAGATAGTGTCAACGAAATGGAAATTAACAAAGCCGATCAA
CAGTTATCAAAATTAGGCCACTCCGTGCATGCTGGATGAACT..................
>Chr2A:13167915-13168522
TTGGAATACCGTTTTTGTTCCTGAGCCGGAAGATCCCATACGCACTTTTCCAGCGACAAA
ACTCCAAAAAGAACACAGATATGCTCATCCGACGGCGTAAAC..................
......
```

---

## Step 4. Run BLAST Search

Search the extracted exon sequences against the NCBI nucleotide database.

```bash
blastn \
-query "${PREFIX}_exons.seqs"  \
-db core_nt \
-outfmt 6 \
-num_threads 32 \
-out "${PREFIX}_exon_blast_results.tsv"

```

Output:

```text
Chr2A_13150000_exon_blast_results.tsv
Chr3A_2506000_exon_blast_results.tsv
Chr6B_ 27199000_exon_blast_results.tsv
```
Example:

```text
$more Chr2A_13150000_exon_blast_results.tsv
```
```.tsv file
Chr2A:13154242-13156423 CP137581.1      86.735  294     24      3       3       281     10615510        10615217        1.19e-80        313
Chr2A:13167915-13168522 CP137580.1      91.150  113     10      0       7       119     13473675        13473563        2.64e-33        154
.......
```
---

## Step 5. Extract Unique Accession IDs

Since ```blastn``` outputs only NCBI accession IDs, we must retrieve the descriptions for each ID. We will first take out the all unique NCBI accesion IDs from the output.

```bash
cut -f 2 "${PREFIX}_exon_blast_results.tsv" > "${PREFIX}_IDs.tsv"
```

Output:

```text
Chr2A_13150000_IDs.tsv
Chr3A_2506000_IDs.tsv
Chr6B_27199000_IDs.tsv
```
Example:
```text
$more Chr2A_13150000_IDs.tsv
```

```_ids.tsv
NM_001365546.1
XM_987654
PP763331.1
NM_001361952.1
CP137586.1
....
```
---

## Step 6. Retrieve Hit Annotations

Obtain accession titles and descriptions for the BLAST hits.

```bash
 blastdbcmd -db core_nt \
 -entry_batch "${PREFIX}_IDs.tsv" \
 -outfmt "%a\t%t" > "${PREFIX}_titles.tsv"

```

Output:

```text
Chr2A_13150000_titles.tsv
Chr3A_2506000_titles.tsv
Chr6B_27199000_titles.tsv
```

Example:
```text
$more Chr2A_13150000_titles.tsv
```

```text
NM_001365546.1        Zea mays Auxin response factor 3 (LOC100502387), mRNA
XM_987654        ABC transporter family protein
PP763331.1        Cymbopogon flexuosus clone 1907 auxin response factor ARF12 mRNA, complete cds
NM_001361952.1        Zea mays auxin response factor 11 (ARF11) gene, complete cds
CP137586.1        Eragrostis tef cultivar Dabbi chromosome 3A
.....
```

---

---

## Step 7. Merge BLAST Results, Annotations, and Exon Information in R

After obtaining the BLAST results and accession descriptions, I used R to combine:

1. BLAST alignment results (`exon_blast_results.tsv`)
2. Accession descriptions (`accession_titles.tsv`)
3. Genome annotation information (`genome.gtf`)

This step creates a comprehensive table containing:
- Exon IDs, BLAST hits, Functional descriptions, Alignment statistics (e-value, percent identity, bitscore, etc.)

### Required Files
> Just download the files from linux environment and load in R. 
```text
exon_blast_results.tsv
accession_titles.tsv
genome.gtf
```

### Load Required Package

```r
install.packages("tidyverse")
library(tidyverse)
```

### Read Input Files

```r
# load BLASTn output in R environemt
exon_blast_result <- read.table("exon_blast_results.tsv")

# load accession descriptions in R
accession_titles <- read_tsv("accession_titles.tsv", col_names = FALSE) %>%
  separate(X1, into = c("id", "description"),sep = "\\\\t|\\t")

# load genome annotation file in R
gtf_file <- read.table("genome.gtf", sep = "\t", header = FALSE)

# setting the column names of gtf_file
colnames(gtf_file) <- c(
  "seqid", "source", "type",
  "start", "end", "score",
  "strand", "phase", "attribute")
```

### Prepare the GTF Annotation

Create a genomic location identifier that can be matched to the exon coordinates extracted earlier.

```r
gtf_file <- gtf_file %>%
  mutate(location = paste0(start, "-", end)) %>%
  select(location, attribute, type)
```

### Assign BLAST Column Names

The BLAST output was generated using `-outfmt 6`, which contains 12 standard columns.

```r
outfmt6_colnames <- c(
  "qseqid", "sseqid", "pident",
  "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart",
  "send", "evalue", "bitscore")

colnames(exon_blast_result) <- outfmt6_colnames
```

### Merge All Information

```r
final_df <- exon_blast_result %>%
 left_join(
    accession_titles,
    by = c("sseqid" = "id"),
    relationship = "many-to-many") %>%
  mutate(
    location = str_split_i(qseqid, ":", 2)) %>%
  left_join(
    gtf_file, by = "location",relationship = "many-to-many") %>%
  mutate(
    gene_id = sub(
      ".*gene_id ([^;]+);.*",
      "\\1",
      attribute)) %>%
  select(
    !c(location, sstart, send, attribute)) %>%
  filter(
    type == "exon",
    !str_detect(description, "PREDICTED"),
    !str_detect(description, "hypothetical"),
    !str_detect(description, "chromosome"),
    !str_detect(description, "predicted protein")) %>%
  filter(
    length > 30,
    evalue < 1e-10) %>%
  distinct(
    qseqid,
    description,
    .keep_all = TRUE) %>%
  group_by(qseqid) %>%
  slice_head(n = 10) %>%
  ungroup()
```

### Explanation of the R Workflow

| Step | Purpose |
|--------|----------|
| Join BLAST hits with accession titles | Adds functional descriptions to each BLAST hit. |
| Extract exon coordinates from `qseqid` | Creates a common key for merging. |
| Join with GTF annotations | Links each exon to its genomic annotation. |
| Extract `gene_id` | Retrieves the corresponding gene identifier from the GTF attributes column. |
| Remove low-quality annotations | Excludes predicted, hypothetical, and chromosome-level descriptions. |
| Apply alignment filters | Retains hits with alignment length > 30 bp and e-value < 1e-10. |
| Remove duplicate descriptions | Keeps unique functional annotations per exon. |
| Retain top hits | Limits the output to the first 10 hits per exon. |

### Output

```r
head(final_df)
```
Example:

```text
qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	description	exon
Chr3B:8825429-8826763	CP137587.1	80.846	851	105	23	499	1303	24628072	24627234	4.28E-171	616	Eragrostis tef cultivar Dabbi chromosome 3B	Chr3B_Cda01350.1-2
Chr3B:8825429-8826763	CP137586.1	81.356	708	76	26	652	1312	26067452	26066754	7.41E-144	525	Eragrostis tef cultivar Dabbi chromosome 3A	Chr3B_Cda01350.1-2
Chr3B:8825429-8826873	CP137587.1	79.49	980	140	26	499	1431	24628072	24627107	2.75E-178	640	Eragrostis tef cultivar Dabbi chromosome 3B	Chr3B_Cda01350.2-2
Chr3B:8825429-8826873	CP137586.1	81.356	708	76	26	652	1312	26067452	26066754	8.04E-144	525	Eragrostis tef cultivar Dabbi chromosome 3A	Chr3B_Cda01350.2-2
Chr3B:8833583-8833738	CP137586.1	84.314	153	12	7	1	153	26059337	26059197	1.05E-28	139	Eragrostis tef cultivar Dabbi chromosome 3A	Chr3B_Cda01352.1-1
....
```


This table represents the final candidate-gene annotation dataset used for downstream biological interpretation.

---

# Complete Linux Script

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
START=13150000
END=13170000

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


        ▼
Retrieve Annotations
  
```
