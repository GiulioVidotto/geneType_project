# geneTypes Package

The `geneTypes` package provides a comprehensive framework for representing various gene types in the human genome using S4 classes. It supports modeling both protein-coding and non-coding genes, along with their products (e.g., proteins, long non-coding RNAs, microRNAs, rRNAs, and tRNAs).

In the human genome, nearly 100,000 genes have been identified. Among these, approximately 20,000-25,000 encode proteins, while the remaining genes give rise to various types of RNA molecules. The package models this diversity by providing specific classes for each gene type:

- `proteinCodingGene`: Represents genes that encode proteins.
- `longNonCodingRNAGene`: Represents genes that encode long non-coding RNAs.
- `microRNAGene`: Represents genes that encode microRNAs.
- `rRNAgene`: Represents genes that encode ribosomal RNAs.
- `tRNAgene`: Represents genes that encode transfer RNAs.

## Installation

### Installation via Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("geneTypes")
```

### Installation from GitHub

```r
devtools::install_github("GiulioVidotto/geneTypes")
```

## Loading the Package

Before using the package, load the necessary libraries:

```r
library(GenomicRanges)
library(Biostrings)
library(geneTypes)
```

## Package Structure and Usage

### Virtual Classes

- **gene**: Captures the common attributes of all genes, such as Ensembl ID, HUGO symbol, full gene name, description, and gene structure.
    - **Methods:**
        - `getGeneId(object)`
        - `getSymbol(object)` / `setSymbol(object, new_gene_symbol)`
        - `getFullName(object)` / `setFullName(object, new_gene_full_name)`
        - `getDescription(object)` / `setDescription(object, new_description)`
        - `getStructure(object)` / `setStructure(object, new_structure)`

- **product**: A virtual class for gene products that includes Ensembl ID as a common attribute.

### Specific Classes and Examples

#### 1. proteinCodingGene

Represents protein-coding genes, including a list of `protein` objects.

```r
gene1 <- proteinCodingGene(
  proteins = list(
    protein(id = "ENSP01234567891",
            protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
            description = "It has a specific function"),
    protein(id = "ENSP19876543210",
            protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
            description = "It has a specific function")
  ),
  geneID = "ENSG01234567891",
  gene_symbol = "GENE1",
  full_gene_name = "GENE ONE",
  description = "This gene has a specific function",
  structure = GRanges("chr1:1-1000")
)
```

#### 2. longNonCodingRNAGene

Represents genes that encode long non-coding RNAs.

```r
gene2 <- longNonCodingRNAGene(
  LongNonCodingRNAs = list(
    longNonCodingRNA(id = "ENST01234567891",
                     long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                     regulatory_mechanism = "Description of the regulatory mechanism"),
    longNonCodingRNA(id = "ENST19876543210",
                     long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                     regulatory_mechanism = "Description of the regulatory mechanism")
  ),
  geneID = "ENSG01234567891",
  gene_symbol = "LINC01018",
  full_gene_name = "LONG INTERGENIC NON-PROTEIN CODING RNA 01018",
  description = "This gene has a specific function",
  structure = GRanges("chr1:1-1000")
)
```

#### 3. microRNAGene

Represents genes that encode microRNAs.

```r
gene3 <- microRNAGene(
  microRNAs = list(
    microRNA(id = "ENST01234567891",
             microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
             silencing_mechanism = "The microRNA silences specific targets"),
    microRNA(id = "ENST19876543210",
             microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
             silencing_mechanism = "The microRNA silences specific targets")
  ),
  geneID = "ENSG01234567891",
  gene_symbol = "MIR 12",
  full_gene_name = "MICRORNA 23",
  description = "This gene has a specific function",
  structure = GRanges("chr1:1-1000")
)
```

## Contact

For any questions or further information, please contact:

**Giulio Vidotto**  
Email: [vidottogiulio@libero.it]
