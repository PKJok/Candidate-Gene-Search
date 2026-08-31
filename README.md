# Searching Candidate Gene
### I have provided full bash and R scripts below this README.md file.
> ***Note*** : I am using local HPCC server's NCBI database to run BLASTn. This contains bash scripts and R codes to achieve the objective.
##📖 Overview

This workflow identifies potential candidate genes within a genomic region of interest by:

1. Define your genomic interval/ know your region of interest.
2. Extracting exon coordinates from a genome annotation file.
3. Retrieving exon sequences from the reference genome.
4. Running BLAST searches against the NCBI nucleotide database.
5. Obtaining annotations for BLAST hits.
6. Prioritizing candidate genes based on biological relevance.

---

# 🧬 Example Region

Suppose we found GWAS-consistent regions on chromosome ```Chr2A```, ```Chr3A```, and ```Chr6B``` at position ```13,160,000 bp```, ```2,516,000 bp```, ```27,209,000 bp``` respectively and we want to identify the genes located within ```20 kb``` of these regions. 

| Chromosome | Position (bp) | Position - 10kb (Start) | Position + 10kb (End) |
|------------|---------------|------------------|------------------|
|Chr2A | 13,160,000 | 13,150,000 | 13,170,000 |
| Chr3A | 2,516,000 | 2,506,000 | 2,526,000 |
| Chr6B | 27,209,000 | 27,199,000 | 27,219,000 |

> ***Note***: Replace these coordinates with the genomic interval associated with your QTL, GWAS consistent region, marker interval, or region of interest.

---

# 📝 Required Input Files

## 1. Reference Genome

Genome sequence in FASTA format.

```text
genome.fa
```


## 2. Genome Annotation File

Genome annotation in GTF/GFF format.

```text
genome.gtf
```


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
>***Note*** : On HPCC server I used, NCBI database is host-integrated with ```blast+ version 2.16.0 or above```. So, when you load the ```blast+/2.16.0``` NCBI database is also loaded.

```bash
module load blast+/2.16.0
blastdbcmd -db core_nt -info
```

---
---
# 💻 Workflow for Bash Script
>***Note***: Whole script is provided below this ```README.md``` file.
## Step 1. Define the Genomic Region

We formed a bash array contatining 3 string elements. Each element holds space-separated genomic coordinates. 

```bash
regions=(
  "Chr2A 13150000 13170000"
  "Chr3A 2506000 2526000 "
  "Chr6B 27199000 27219000" )
```

## Step 2. Running the loop to extract the sequences from all the regions
We also formed a ```for``` ..... ```done``` loop that iterates over each element in the ```regions``` array.

```bash
for region in "${regions[@]}" ; do
    read -r CHROM START END <<< "$region"
    PREFIX="${CHROM}_${START}"

    echo "$PREFIX"
    echo "$CHROM"
    echo "$START"
    echo "$END"
    ....
 done

```
### ➡️Understanding the ```loop``` code
``` bash
read -r CHROM START END <<< "$region"
```
#### This ```read``` command with a Here-String (```<<<```) splits the ```$region``` string by spaces into three distinct variables: ```$CHROM```,```$START```,and ```$END```.
---
``` bash
PREFIX="${CHROM}_${START}"
````
#### Combines the chromosome name and start coordinate with an underscore (e.g.,```Chr2A_13150000```) to create a unique prefix for naming files.
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
### ➡️ Understanding the `awk` Command ###

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
| `$1 == chrom` | Chromosome name ```(Column 1)``` matches the target chromosome. |
| `$3 == "exon"` | Feature type ```(Column 3)``` is an exon. |
| `$4 >= start` | Exon start coordinate ```(Column 4)``` is greater than or equal to the region start position. |
| `$5 <= end` | Exon end coordinate ```(Column 5)``` is less than or equal to the region end position. |

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

## Step 4. Extract Exon Sequences

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

## Step 4.⚡Run BLAST Search

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
Chr6B_27199000_exon_blast_results.tsv
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

## Step 5. 📝 Extract Unique Accession IDs

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

## Step 6. 📝 Retrieve Hit Annotations

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


# ⚡Automating Candidate Gene Annotation in R
>***Note***: Full R script is provided below this ```README``` file
## Overview

This script automatically:

1. Detects all BLAST result files in the working directory.
2. Loads the corresponding accession descriptions.
3. Merges BLAST results with genome annotations.
4. Filters low-quality and uninformative hits.
5. Saves a candidate gene table for each genomic region.

---
> **Note:** Before running the R workflow, it is recommended to create a dedicated R project to keep all files organized and ensure reproducibility.

### 🚨 IMPORTANT: Creating a New R Project

1. In **RStudio**, click **File** → **New Project...**
2. Select **New Directory**.
3. Click **New Project**.
4. Enter a name for your project.
5. Choose the location where you would like the project folder to be created.
6. Click **Create Project**.


### 📌Organizing Project Files

After the project is created:

1. Open the newly created project.
2. Copy or upload all required input files into the project directory, including:
   - BLAST output files (`*_exon_blast_results.tsv`)
   - Accession title files (`*_titles.tsv`)
   - Genome annotation file (`genome.gtf`)
     
### Example:
```text
Chr2A_13150000_exon_blast_results.tsv
Chr3A_2506000_exon_blast_results.tsv
Chr6B_27199000_exon_blast_results.tsv
Chr2A_13150000_titles.tsv
Chr3A_2506000_titles.tsv
Chr6B_27199000_titles.tsv
genome.gtf
```

## Step 1. Load Required Package

```r
library(tidyverse)
```

---

## Step 2. Identify File Prefixes

```r
prefixes <- list.files(pattern = "_exon_blast_results.tsv") %>%
  str_remove("_exon_blast_results.tsv")
```

This extracts prefixes from BLAST result files.

### Example

Input files:

```text
Chr2A_13150000_exon_blast_results.tsv
Chr3A_2506000_exon_blast_results.tsv
Chr6B_27199000_exon_blast_results.tsv
```

Extracted prefixes:

```output
Chr2A_13150000
Chr3A_2506000
Chr6B_27199000
```

---

## Step 3. Load and Prepare the GTF File

```r
gtf_file <- read.table("genome.gtf", sep = "\t", header = FALSE)

colnames(gtf_file) <- c(
  "seqid","source","type","start","end",
  "score","strand","phase","attribute"
)

gtf_file <- gtf_file %>%
  mutate(location = paste0(start, "-", end)) %>%
  select(location, attribute, type)
```

The `location` column is used to match exon coordinates between the BLAST results and genome annotation.

---

## Step 4. BLAST Column Names

```r
outfmt6_colnames <- c(
  "qseqid", "sseqid", "pident", "length",
  "mismatch", "gapopen", "qstart", "qend",
  "sstart", "send", "evalue", "bitscore"
)
```

These are the standard column names for BLAST output generated using `-outfmt 6`.

---

## Step 5.Candidate Gene Processing Function

```r
process_candidate_genes <- function(prefix){
  
  exon_blast_result <- read.table(
    paste0(prefix, "_exon_blast_results.tsv"))
    
  colnames(exon_blast_result)<-outfmt6_colnames
  
  accession_titles <- read_tsv(paste0(prefix, "_titles.tsv"),   # getting the titles loaded 
                               col_names = FALSE) %>%
    separate(X1,
             into = c("id", "description"),
             sep = "\\\\t|\\t")


    final_df<- exon_blast_result %>%
    left_join(accession_titles,
            by = c("sseqid" = "id"),
            relationship = "many-to-many") %>%
    mutate(location = str_split_i(qseqid, ":", 2)) %>%
    left_join(gtf_file,
            by = "location",
            relationship = "many-to-many") %>%
    select(!c(location, sstart, send, attribute)) %>%
    filter(
    type == "exon",
    !str_detect(description, regex("PREDICTED|hypothetical|chromosome|predicted protein",
                                   ignore_case = TRUE))) %>%
    filter(length > 30,evalue < 1e-10) %>%
    distinct(qseqid,description,.keep_all = TRUE) %>%
    group_by(qseqid) %>%
    slice_head(n = 5) %>%
    ungroup()
    
    write_tsv(final_df,paste0(prefix, "_candidate_genes.tsv") )
    
    return(final_df)
}
```

The function performs the following steps:

| Step | Purpose |
|--------|----------|
| Load BLAST results | Reads `*_exon_blast_results.tsv` |
| Load accession descriptions | Reads `*_titles.tsv` |
| Merge data | Combines BLAST hits with descriptions |
| Match exon coordinates | Links BLAST hits to GTF annotations |
| Filter results | Removes predicted, hypothetical, and chromosome-level matches |
| Quality filtering | Keeps hits with `length > 30` and `evalue < 1e-10` |
| Remove duplicates | Keeps unique descriptions per exon |
| Retain top hits | Keeps the first 5 hits per exon |
| Save output | Writes `*_candidate_genes.tsv` |

---

## Step 6. Save the Final Results

```r
write_tsv(
  final_df,
  paste0(prefix, "_candidate_genes.tsv")
)
```

### Example Outputs

```text
Chr2A_13150000_candidate_genes.tsv
Chr3A_2506000_candidate_genes.tsv
Chr6B_27199000_candidate_genes.tsv
```

---

## Run the Analysis for All Regions

```r
lapply(prefixes, process_candidate_genes)
```

This automatically processes every BLAST result file in the working directory and generates a corresponding candidate gene table.

### Final Output

Each output file contains:

- Exon IDs
- BLAST accession IDs
- Functional descriptions
- Alignment statistics
- Candidate gene annotations

These files serve as the final candidate gene datasets for downstream biological interpretation.

### Example:

```text
qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	evalue	bitscore	description	type
LG03:52986100-52986214	NM_001409131.1	87.5	88	8	3	23	110	1.2e-16	99	Oryza sativa Japonica Group cyclin-B1-1-like (LOC4327547), transcript variant 1, mRNA	exon
LG03:52986100-52986214	NM_001409132.1	87.5	88	8	3	23	110	1.2e-16	99	Oryza sativa Japonica Group cyclin-B1-1-like (LOC4327547), transcript variant 2, mRNA	exon
LG03:52988801-52988950	NM_001143365.1	90.604	149	14	0	1	149	1.65e-46	198	Zea mays uncharacterized LOC100216987 (LOC100216987), mRNA	exon
....
```


This table represents the final candidate-gene annotation dataset used for downstream biological interpretation.

---

# Complete Linux Script
>***Note***: This contains script upto the blastn and extraction of acccession titles from accession Ids. 

```bash
#!/bin/bash
#SBATCH --job-name="Ultimate blast searching"
#SBATCH -p bigmem
#SBATCH -t 168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=prajjwal.koirala10@okstate.edu
#SBATCH --mail-type=end

# module load
#echo "Loading modules .........................................."
module load bedtools
module load blast+/2.16.0
which blastn
blastn -version

# Important Files
echo "lodaing files .........................................."
GFF_FILE="/scratch/pkjok/Candidate_gene/genome.gtf"
REF_GENOME="/scratch/pkjok/Candidate_gene/genome.fa"


# Input your Chromosome and region of interest
regions=(
  "Chr2A 13150000 13170000"
  "Chr3A 2506000 2526000 "
  "Chr6B 27199000 27219000" )

for region in "${regions[@]}" ; do
    read -r CHROM START END <<< "$region"
    PREFIX="${CHROM}_${START}"

    echo "$PREFIX"
    echo "$CHROM"
    echo "$START"
    echo "$END"

    # extracting the exons coordinates with in our region of interest and saving as bed file
    echo "extracting the bed file................"
    awk -v OFS="\t" -v chrom="$CHROM" -v start="$START" -v end="$END" \
    '$1==chrom && $3=="exon" &&  $4>=start && $5<=end {print $1, $4, $5}' "$GFF_FILE" > "${PREFIX}.bed"

    echo "extracting the sequences from bed file........"
    bedtools getfasta -fi "$REF_GENOME" -bed "${PREFIX}.bed" -fo "${PREFIX}_exons.seqs"

    # running the blastn
    echo "running blastn ..............................."
     blastn \
     -query "${PREFIX}_exons.seqs" \
     -db core_nt \
     -outfmt 6 \
     -num_threads 32 \
     -out "${PREFIX}_exon_blast_results.tsv"

    # getting IDS for getting the accession title
    cut -f 2 "${PREFIX}_exon_blast_results.tsv" > "${PREFIX}_IDs.tsv"

    # get accession titles
    echo "Getting accession titles ..................."
    blastdbcmd -db core_nt -entry_batch "${PREFIX}_IDs.tsv" -outfmt "%a\t%t" > "${PREFIX}_titles.tsv"

done

echo "Analysis completed."
```

---
# Complete R Script
```
# load the library
library(tidyverse)

# list out our blast result files and extract the prefixes from them
prefixes<-list.files(pattern = "_exon_blast_results.tsv")%>%
  str_remove("_exon_blast_results.tsv")


# load the gtf file

gtf_file<-read.table("genome.gtf", sep = "\t", header = F)
colnames(gtf_file)<-c("seqid","source","type","start","end","score","strand",
                      "phase","attribute")

# load the column names of outfmt6
outfmt6_colnames<-c("qseqid", "sseqid", "pident", "length", "mismatch", 
                    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")


gtf_file<-gtf_file%>%
  mutate(location=paste0(start,"-",end))%>%
  select(c(location,attribute,type))

process_candidate_genes <- function(prefix){
  
  exon_blast_result <- read.table(
    paste0(prefix, "_exon_blast_results.tsv"))
    
  colnames(exon_blast_result)<-outfmt6_colnames
  
  accession_titles <- read_tsv(paste0(prefix, "_titles.tsv"),
                               col_names = FALSE) %>%
    separate(X1,
             into = c("id", "description"),
             sep = "\\\\t|\\t")


    final_df<- exon_blast_result %>%
    left_join(accession_titles,
            by = c("sseqid" = "id"),
            relationship = "many-to-many") %>%
    mutate(location = str_split_i(qseqid, ":", 2)) %>%
    left_join(gtf_file,
            by = "location",
            relationship = "many-to-many") %>%
    select(!c(location, sstart, send, attribute)) %>%
    filter(
    type == "exon",
    !str_detect(description, regex("PREDICTED|hypothetical|chromosome|predicted protein",
                                   ignore_case = TRUE))) %>%
    filter(length > 30,evalue < 1e-10) %>%
    distinct(qseqid,description,.keep_all = TRUE) %>%
    group_by(qseqid) %>%
    slice_head(n = 5) %>%
    ungroup()
    
    write_tsv(final_df,paste0(prefix, "_candidate_genes.tsv") )
    
    return(final_df)
}

lapply(prefixes, process_candidate_genes) 
```
