# geneTypes Package

Il pacchetto `geneTypes` fornisce un framework completo per rappresentare i vari tipi di geni nel genoma umano utilizzando le classi S4. Supporta la modellazione sia dei geni codificanti proteine che dei geni non codificanti, insieme ai loro prodotti (ad esempio, proteine, RNA lunghi non codificanti, microRNA, rRNA e tRNA).

Nel genoma umano sono stati identificati quasi 100.000 geni. Tra questi, circa 20.000-25.000 codificano per proteine, mentre il resto produce vari tipi di molecole di RNA. Il pacchetto modella questa diversità fornendo classi specifiche per ciascun tipo di gene:

- `proteinCodingGene`: Rappresenta i geni che codificano per proteine.
- `longNonCodingRNAGene`: Rappresenta i geni che codificano per RNA lunghi non codificanti.
- `microRNAGene`: Rappresenta i geni che codificano per microRNA.
- `rRNAgene`: Rappresenta i geni che codificano per rRNA.
- `tRNAgene`: Rappresenta i geni che codificano per tRNA.

## Installation

### Installazione tramite Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("geneTypes")
```

### Installazione da GitHub

```r
devtools::install_github("GiulioVidotto/geneTypes")
```

## Loading the Package

Prima di utilizzare il pacchetto, carica le librerie necessarie:

```r
library(GenomicRanges)
library(Biostrings)
library(geneTypes)
```

## Package Structure and Usage

### Classi Virtuali

- **gene**: Cattura gli attributi comuni a tutti i geni, come l'Ensembl ID, il simbolo HUGO, il nome completo del gene, la descrizione e la struttura del gene.
    - **Metodi:**
        - `getGeneId(object)`
        - `getSymbol(object)` / `setSymbol(object, new_gene_symbol)`
        - `getFullName(object)` / `setFullName(object, new_gene_full_name)`
        - `getDescription(object)` / `setDescription(object, new_description)`
        - `getStructure(object)` / `setStructure(object, new_structure)`

- **product**: Una classe virtuale per i prodotti genici che include l'Ensembl ID come attributo comune.

### Classi Specifiche ed Esempi

#### 1. proteinCodingGene

Rappresenta i geni codificanti per proteine, includendo una lista di oggetti di classe `protein`.

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

Rappresenta i geni che codificano per RNA lunghi non codificanti.

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

Rappresenta i geni che codificano per microRNA.

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

## Contributing

Contributi e feedback sono ben accetti! Per favore apri una issue o invia una pull request sul repository GitHub.

## License

[Inserisci qui le informazioni sulla licenza]

## Contact

Per qualsiasi domanda o ulteriore informazione, contatta:

**Giulio Vidotto**  
Email: [vidottogiulio@libero.it](mailto:vidottogiulio@libero.it)
