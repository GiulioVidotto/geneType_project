# GeneTypes package


**Description of the package**

  As of today, in the human genome, there have been identified around 100.000 genes. One of the main differences, among all
  of the genes, is regarding the functions of the product or products transcribed and, in some cases, translated.
  For example, it is estimated that about 20.000 to 25,000 genes codify for proteins, and for this reason, these genes
  are classified as "Protein Coding Genes". On the other end, there are other genes classified as "Non-Protein
  Coding Gene" whose products are RNA molecules that have different functions such as, for example, regulatory function.
  Among the products of these genes there are the so-called: microRNA (miRNA), long non-coding RNA (lncRNA), ribosomal 
  RNA (rRNA) and transfer RNA (tRNA).
  The 'r Rpackage("geneTypes")' package provides a set of S4 classes for each different gene type described before. 
  Because there is some information that is common to all the genes, it was created a virtual class called "gene"
  that represents the genes in general. For the sake of a clear and organized structure of the products of each gene,
  a class for each one of them has been created and inherited from a virtual class called "product".


**gene class**


**Description of the class**

  This class represents the genes in a general way. Indeed, it stores all the information that is common to all the genes.
  In particular, it defines the Ensmbl ID, the HUGO symbol, the name of the gene, a short description, and all the information
  about the structure of the gene (e.g. the chromosome, start and stop, etc.). It is important to remember that it is a virtual
  class, so this means that all the attributes and methods of the class will be inherited by all the specific gene classes. 
  Because, as said before, this class represents the genes in a general way, the only additional constraint is regarding the
  Ensembl ID of the gene. It must start with "ENSG" and be followed by 11 digits otherwise, the following error will rise:
  "The gene enemblID is not correct".


**Methods of the class**

### - getGeneId

  This is a method of the "gene" class, called **"getGeneId"**, that can be used to access information about the Ensembl ID gene. 
  The get method is the only method specified for the Ensembl ID gene because usually the ID of a gene should not change, even
  if the gene is updated.



### - getSymbol

  This is a method of the "gene" class called **"getSymbol"**. It can be used to access information about the the HUGO gene symbol.
  This method has only one parameter, called "object", and it is the gene for which the user wants to know the HUGO gene symbol.



### - setSymbol

  This is a method of the "gene" class, called **"setSymbol"**, that can be used to change the information about the HUGO gene symbol.
  This method requires to specify two paramters. the first one, called "object" which rapresents a gene and the second one, called
  "new_gene_symbol" which is the new HUGO symbol for the gene.



### - getFullName

  This is a method of the "gene" class, called **"getFullName"**, that allows the user to extract the information about the full name
  of the gene. This method has only one parameter, called "object", which is a gene.



### - setFullName

  This is a method of the "gene" class, called **"setFullName"**, that allows the user to change the information about the full name of
  the gene. This method has two parameters, the first one is called "object" which represents a gene, and the second parameter is called
  "new_gene_full_name" which, as the name suggests, is the new full gene name.



### - getDescription

  This is a method of the "gene" class, called **"getDescription"**, and it allows the user to extract the description of a specific 
  gene. This method has only one parameter, which is called "object", and it is the gene for which the user wants to inspect the
  description.



### - setDescription

  This is a method of the "gene" class, called **"setDescription"**, and it allows the user to change the description of a specific gene. It
  takes as parameters, a gene ("object") and a new description for the gene ("new_description"). The new description must be of class "character",
  if not the following error will rise: "The description is not correct".



### - getStructure

  This is a method of the "gene" class, called **"getStructure"**, and it allows the user to extract from a specific object all the information on
  the structure of the gene. It reqeuires only one parameter, which is called "object" and it is the gene for which the user wants to know the
  structure.



### - setStructure

  This is a method of the "gene" class, called **"setStructure"**, that is useful to extract from a specific gene all the information about the structure.
  It requires two paramters: the "object" parameter, which is the gene, and the "new_structure" parameter which is the new structure of the gene. This 
  second parameter must be of class GRanges, otherwise the method will rise the followin "The new defined structure of the gene is not correct".


# - product class



**Description of the class**

  Because for most of the genes, the products are lists of specific elements, to create and specify a clear structure of these lists it was decided to define a new class called "product".
  This class is a virtual class that represents the products of the different genes in a general way. Indeed, in this virtual class, the only attribute specified is the Ensembl ID of
  the product of the genes. This is because it is the only attribute common to all the products (proteins, microRNA, lncRNA, rRNA, tRNA).


  From this virtual class, other classes represent the specific products for the different types of genes:

### - Protein

  This class is called "Protein" and it inherits from the virtual class "product". It inherits the attribute called Ensembl ID from the virtual class "product" and it has also other two 
  attributes which are the sequence of the protein ("protein_sequence") and the description of the protein ("description"). To further specify the object of the class, two additional
  constraints have been added:
    * The Ensembl ID of the protein must start with "ENSP" followed by 11 digits, otherwise, it will raise the following: "The protein enemblID is not correct".
    * The aminoacid of the protein sequence must be valid, otherwise, it will raise the following: "In the protein sequence there are one or more invalid aminoacids".


### - LongNonCodingRNA

  This class is called "LongNonCodingRNA" and it inherits from the virtual class "product". It inherits the attribute called Ensembl ID from the virtual class "product" and it has also other two 
  attributes which are the sequence of the LongNonCodingRNA ("long_non_coding_RNA_sequence") and the regulatory mechanism of the LongNonCodingRNA ("regulatory_mechanism"). To further specify 
  the object of the class, one more additional constraint has been added:
    * The Ensembl ID of the Long Non-Coding RNA must start with "ENST" followed by 11 digits, otherwise, it will raise the following: "The long noncoding RNA enemblID is not correct".


### - microRNA

  This class is called "MicroRNA" and it inherits from the virtual class "product". It inherits the attribute called Ensembl ID from the virtual class "product" and it has also other two 
  attributes which are the sequence of the MicroRNA ("microRNA_sequence") and the description of the silencing mechanism of the MicroRNA ("silencing_mechanism"). To further specify the
  object of the class, one additional constraint has been added:
    * The Ensembl ID of the microRNA must start with "ENST" followed by 11 digits, otherwise, it will raise the following: "The microRNA enemblID is not correct".
  

### - tRNA

  This class is called "tRNA" and it inherits from the virtual class "product". It inherits the attribute called Ensembl ID from the virtual class "product" and it has also other three 
  attributes which are the sequence of the tRNA ("tRNA_sequence"), the aminoacid carried by the tRNA ("aminoacid"), and the anticodon present on the tRNA ("anticodon").
  To furthermore specify the object of the class, other three additional constraints have been added:
    * The Ensembl ID of the tRNA must start with "ENST" followed by 11 digits, otherwise, it will raise the following: "The tRNA enemblID is not correct".
    * The aminoacid must be valid, otherwise, it will raise the following: "Invalid aminoacid".
    * The anticodon must have a length of 3 base pairs, otherwise, it will raise the following: "Invalid anticodon".
  

### - rRNA

  This class is called "rRNA" and it inherits from the virtual class "product". It inherits the attribute called Ensembl ID from the virtual class "product" and it has also other two 
  attributes which are the sequence of the rRNA ("rRNA_sequence") and the ribosomal subunit ("ribosomal_subunit").
  To furthermore specify the object of the class, other three additional constraints have been added:
    * The Ensembl ID of the rRNA must start with "ENST" followed by 11 digits, otherwise, it will raise the following: "The rRNA enemblID is not correct".
    * The subunit must be valid, otherwise, it will raise the following: "The ribosomal subunit is not correct".
    * Only one subunit, otherwise, it will raise the following: "The ribosomal subunit is not correct".


# proteinCodingGene class


**Description of the class**

  This class rapresents the genes whose transcript(s) codify for proteins. This class inherits all the attributes of the virtual class "gene" and it furthermore
  specifies the type of the gene by adding information about the protein(s), codified by this type of gene. The "proteins" attribute is a list of elements of 
  the class "Protein". This class inherits from the virtual class "product" and it has different parameters such as the "id" of the protein, the sequence of the
  protein ("protein_sequence") and lastly the description of the protein ("description").
  This class inherits the attributes but also the constrains of the "gene" class. Other than that, the objects of this class must safisfy other constrains such as:
    * The HUGO gene symbol must be uppercase
    * The full gene name must be uppercase
    * The elemnts of the list of proteins must be of class "protein"
    * No replicates in the list of proteins

  An example of how to create an object of class "proteinCodingGene":
```{r}
gene1 <- proteinCodingGene(
        proteins = list(
            protein(id = "ENSP01234567891", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function"),
            protein(id = "ENSP19876543210", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
```

**Methods of the class**

### - getGeneId

```{r}
getGeneId(gene1)
```

### - getSymbol

```{r}
getSymbol(gene1)
```

### - setSymbol

  This is a method inherited from the "gene" class, For the class "ProteinCodingGene" there have been added additional constrains regarding the new HUGO
  gene symbol. Indeed it must be of class "character" and also must be uppercase, otherwise, the following will raise: "Invalid gene symbol".

```{r}
gene1 <- setSymbol(gene1, "GENE2")
gene1
```

### - getFullName

```{r}
getFullName(gene1)
```

### - setFullName

  This is a method inherited from the "gene" class, For the class "ProteinCodingGene" there have been added additional constrains regarding the new full
  gene name. Indeed it must be of class "character" and also must be uppercase, otherwise, the following will raise: "Invalid full gene name".

```{r}
gene1 <- setFullName(gene1, "GENE TWO")
gene1
```

### - getDescription

```{r}
getDescription(gene1)
```

### - setDescription

```{r}
gene1 <- setDescription(gene1, "This gene has no specific function")
gene1
```

### - getStructure

```{r}
getStructure(gene1)
```

### - setStructure

```{r}
gene1 <- setStructure(gene1, GRanges("chr1:100-200"))
gene1
```

### - getProtein

  This is a method of the class "ProteinCodingGene" class. It allows the user to acccess the information about the protein(s) codified by the gene.

```{r}
getProteins(gene1)
```

### - setProtein

  This is a method of the class "ProteinCodingGene" class. It allows the user to modify the information about a specific protein codified by the gene.
  The new protein must be an object of class "Protein", otherwise, the following will raise: "Invalid protein".

```{r}
gene1 <- setProtein(gene1, "ENSP01234567891", protein(id = "ENSP00000000001", protein_sequence = AAString("ARNDCQ"), description = "It has no specific function"))
gene1
```


# longNonCodingRNAGene class


**Description of the class**

  This class represents the genes whose transcript(s) are classified as long non coding RNA. This class inherits all the attributes
  of the virtual class "gene" and it furthermore specifies the type of the gene by adding information about the long non coding RNA(s).
  The "LongNonCodingRNAs" attribute is a list of elements of the class "LongNonCodingRNA". This class inherits from the virtual class
  "product" and it has different parameters such as the "id" of the long non coding RNA, the sequence of the long non coding RNA 
  ("long_non_coding_RNA_sequence") and lastly its regulator mechanism ("regulatory_mechanism").
  This class inherits the attributes but also the constrains of the "gene" class. Other than that, the objects of this class must safisfy
  other constrains about the HUGO gene symbol and the full name gene depending on the type of long non coding RNA gene:
    * long intergenic non-protein coding RNA
    * LncRNA genes that are antisense to the genomic span of a protein coding gene
    * LncRNA genes that are divergent to (share a bidirectional promoter with) a protein coding gene
    * LncRNA genes that are contained within an intron of a protein coding gene
    * LncRNA genes that overlap a protein coding gene on the same strand

  An example of how to create an object of class "proteinCodingGene":
```{r}
gene2 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(id = "ENST01234567891", long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(id = "ENST19876543210", long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "LINC01018",
    full_gene_name = "long intergenic non-protein coding RNA 01018",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))
```

**Methods of the class**

### - getGeneId

```{r}
getGeneId(gene2)
```

### - getSymbol

```{r}
getSymbol(gene2)
```

### - setSymbol

  This is a method inherited from the "gene" class, For the class "longNonCodingRNAGene" there have been added additional constrains regarding the new HUGO
  gene symbol. Indeed it must be of class "character" and also must contain specific patters of words (depending on the type of lncRNA gene), otherwise, the
  following will raise: "Invalid gene symbol".

```{r}
# The example is made with a long intergenic non-protein coding RNA gene
gene2 <- setSymbol(gene2, "LINC11111")
gene2
```

### - getFullName

```{r}
getFullName(gene2)
```

### - setFullName

  This is a method inherited from the "gene" class, For the class "longNonCodingRNAGene" there have been added additional constrains regarding the new
  full gene name. Indeed it must be of class "character" and also must contain specific patters of words (depending on the type of lncRNA gene), otherwise, the
  following will raise: "Invalid full gene name".

```{r}
# The example is made with a long intergenic non-protein coding RNA gene
gene2 <- setFullName(gene2, "long intergenic non-protein coding RNA 01234")
gene2
```

### - getDescription

```{r}
getDescription(gene2)
```

### - setDescription

```{r}
gene2 <- setDescription(gene2, "This gene has no specific function")
gene2
```

### - getStructure

```{r}
getStructure(gene2)
```

### - setStructure

```{r}
gene2 <- setStructure(gene2, GRanges("chr1:100-200"))
gene2
```

### - getLongNonCodingRNA

  This method of class "longNonCodingRNAGene" which allows the users to extract and inspect the information about the long Non Coding RNAs.

```{r}
getLongNonCodingRNA(gene2)
```

### - setLongNonCodingRNA

  This method of class "longNonCodingRNAGene" which allows the users to change the information about the long Non Coding RNAs.

```{r}
gene2 <- setlongNonCodingRNA(gene2, "ENST01234567891", longNonCodingRNA(id = "ENST00000000001", long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), regulatory_mechanism = "Description of the specific regulatory mechanism"))
gene2
```

### lengthProductLNCRG

  This function helps in calculating the length of the RNA sequence of all the long non coding RNA that are codified by the gene

```{r}
lengthProductLNCRG(gene2)
```


# microRNA class


**Description of the class**

  This class represents the genes whose transcript(s) are classified as microRNA. This class inherits all the attributes
  of the virtual class "gene" and it furthermore specifies the type of the gene by adding information about the 
  microRNA(s). The "microRNAs" attribute is a list of elements of the class "microRNA". This class inherits from the virtual
  class "product" and it has different parameters such as the "id" of the microRNA, the sequence of the microRNA ("microRNA_
  sequence") and lastly its silencing mechanism ("silencing_mechanism").
  This class inherits the attributes but also the constrains of the "gene" class. Other than that, the objects of this class must safisfy
  other constrains:
    * The HUGO gene symbol must start with "MIR" followed by a space and then a number
    * The full name of the gene must start with "microRNA", followed by a space and then a number
    * The elements of the list "microRNAs" must be of class "microRNA"
    * No microRNA replicates

  An example of how to create an object of class "microRNAGene":
```{r}
gene3 <- microRNAGene(
    microRNAs = list(
        microRNA(id = "ENST01234567891", microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), silencing_mechanism = "The microRNA has .."),
        microRNA(id = "ENST19876543210", microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), silencing_mechanism = "The microRNA has ..")),
    geneID = "ENSG01234567891",
    gene_symbol = "MIR 12",
    full_gene_name = "microRNA 23",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000")
    )
```

**Methods of the class**

### - getGeneId

```{r}
getGeneId(gene3)
```

### - getSymbol

```{r}
getSymbol(gene3)
```

### - setSymbol

  This is a method inherited from the "gene" class, For the class "MicroRNAGene" the new gene symbol has to respect
  different conditions, for example it has to be a character object and also start with "MIR " followed by a number.
  If this is not the case, then the following will raise: "Invalid gene symbol".

```{r}
gene3 <- setSymbol(gene3, "MIR 123")
gene3
```

### - getFullName

```{r}
getFullName(gene3)
```

### - setFullName

  This is a method inherited from the "gene" class, For the class "MicroRNAGene" the new full name has to respect
  different conditions such as it has to be a character and also starts with "microRNA" followed by a space and then
  a number.

```{r}
gene3 <- setFullName(gene3, "microRNA 1234")
gene3
```

### - getDescription

```{r}
getDescription(gene3)
```

### - setDescription

```{r}
gene3 <- setDescription(gene3, "This gene has no specific function")
gene3
```

### - getStructure

```{r}
getStructure(gene3)
```

### - setStructure

```{r}
gene3 <- setStructure(gene3, GRanges("chr1:100-200"))
gene3
```

### - getMicroRNA

  This is a method of the class "microRNAgene" class. It allows the users to extract and inspect the information
  about the microRNAs. This object will contain the information about the ID, sequence and silencing mechanism
  of each microRNA defined for a specific gene.

```{r}
getMicroRNA(gene3)
```

### - setMicroRNA

  This method is called "setMicroRNA" which allows the users to change the information about a specific
  microRNA. The user can change the information about the ID, sequence and silencing mechanism of the
  microRNA

```{r}
gene3 <- setMicroRNA(gene3, "ENST01234567891", microRNA(id = "ENST00000000001", microRNA_sequence = RNAString("ACG"), silencing_mechanism = "Its silencing mechanism consists of"))
gene3
```

### - lengthProductMRG

  This function helps in calculating the length of the RNA sequence of all the microRNA that are codified by the gene

```{r}
lengthProductMRG(gene3)
```


# rRNA class


**Description of the class**

  This class rapresents the genes whose transcript(s) are classified as rRNA (ribosomal RNA).. This class inherits all the attributes
  of the virtual class "gene" and it furthermore specifies the type of the gene by adding information about the rRNA. The "rRNA"
  attribute is a single object of class "rRNA". This class inherits from the virtual class "product" and it has different parameters
  such as the "id" of the rRNA, the sequence of the rRNA ("rRNA_sequence") and lastly its ribosomal subunit ("ribosomal_subunit").
  This class inherits the attributes but also the constrains of the "gene" class. Other than that, the objects of this class must safisfy
  other constrains:
    * The HUGO gene symbol must start with "RNA" followed by a the ribosomal subunit.
    * The full name of the gene must start with "RNA," followed by the number of the ribosomal subunit and then the word "ribosomal".

  An example of how to create an object of class "microRNAGene":
```{r}  
gene4 <- rRNAgene(
    rRNA = rRNA(
        id = "ENST01234567891", rRNA_sequence = RNAString("AUUUGA"), ribosomal_subunit = "28S"),
    geneID = "ENSG01234567891",
    gene_symbol = "RNA28S",
    full_gene_name = "RNA, 28S ribosomal",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))
```

**Methods of the class**

### - getGeneId

```{r}
getGeneId(gene4)
```

### - getSymbol

```{r}
getSymbol(gene4)
```

### - setSymbol

  This is a method inherited from the "gene" class, for the class "rRNAGene" the new full name has to respect
  different conditions such as it has to be a character and also starts with "RNA" followed by a valid ribosomal
  subunits

```{r}
gene4 <- setSymbol(gene4, "RNA45S")
gene4
```

### - getFullName

```{r}
getFullName(gene4)
```

### - setFullName

  This is a method inherited from the "gene" class, for the class "rRNAGene" the new full name has to respect
  different conditions such as it has to be a character and also starts with "RNA" followed by a valid ribosomal
  subunits and then "ribosomal".

```{r}
gene4 <- setFullName(gene4, "RNA, 45S ribosomal")
gene4
```

### - setRibosomalSub

  This is a method for the class "rRNAgene" that allows the user to change the ribosomal subunit. It has to respect
  different conditions, for example it has to be a character object and also be a valid ribosomal subunits. If this is
  not the case, then the following will raise: "Invalid Ribosomal Subunit".

```{r}
gene4 <- setRibosomalSub(gene4, "45S")
gene4
```

### - getDescription

```{r}
getDescription(gene4)
```

### - setDescription

```{r}
gene4 <- setDescription(gene4, "This gene has no specific function")
gene4
```

### - getStructure

```{r}
getStructure(gene4)
```

### - setStructure

```{r}
gene4 <- setStructure(gene4, GRanges("chr1:100-200"))
gene4
```

### - getRRNA

  This is a method of the class "rRNAgene" class. It allows the users to extract and inspect the information
  about the rRNA. This object will contain the information about the ID and ribosomal subunit of the rRNA
  for a specific gene.

```{r}
getRRNA(gene4)
```

### - setRRNA

  This method is called "setRRNA" which allows the users to change the information about the rRNA.
  The user can change the information about the ID ad the ribosomal subunit of the rRNA.

```{r}
gene4 <- setRRNA(gene4, rRNA(id = "ENST00000000001", rRNA_sequence = RNAString("AUUUGA"), ribosomal_subunit = "45S"))
gene4
```

### - lengthProductRRG

  This function helps in calculating the length of the RNA sequence of the rRNA that is codified by the gene.

```{r}
lengthProductRRG(gene4)
```


# tRNA class


**Description of the class**

  This class rapresents the genes whose transcript(s) are classified as tRNA (transfer RNA). This class inherits all the attributes
  of the virtual class "gene" and it furthermore specifies the type of the gene by adding information about the tRNA. The "tRNA"
  attribute is a single object of class "tRNA". This class inherits from the virtual class "product" and it has different parameters
  such as the "id" of the tRNA, the sequence of the tRNA ("tRNA_sequence"), the aminoacid carried by the tRNA ("aminoacid") and the
  anticodon sequence ("anticodon").

  An example of how to create an object of class "tRNA":
```{r}  
gene5 <- tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), aminoacid = AAString("A"), anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000"))
```

**Methods of the class**

### - getGeneId

```{r}
getGeneId(gene5)
```

### - getSymbol

```{r}
getSymbol(gene5)
```

### - setSymbol

  This is a method inherited from the "gene" class, for the class "tRNAGene" the new HUGO gene symbol has to respect
  different conditions such as it has to be a character and also starts with "TR" followed by a valid aminoacid symbol,
  a valid anticodon and additional digits.

```{r}
gene5 <- setSymbol(gene5, "TRM-AUU1-1")
gene5
```

### - setAminoacid

  This is a method of the class "tRNAgene". The new aminoacid has to respect different conditions, for example
  it has to be of class AAString and also has to be a valid aminoacid. If this is not the case, then the following
  will raise: "Invalid aminoacid".

```{r}
gene5 <- setAminoacid(gene5, AAString("M"))
gene5
```

### - setAnticodon

  This is a method of the class "tRNAgene". The new anticodon has to respect different conditions, for example
  it has to be of class RNAString and also its length has to be of 3 base. If this is not the case, then the
  following will raise: "Invalid aminoacid".

```{r}
gene5 <- setAnticodon(gene5, RNAString("AUU"))
gene5
```

### - getFullName

```{r}
getFullName(gene5)
```

### - setFullName

  This is a method inherited from the "gene" class, for the class "tRNAGene" the new full name has to respect
  different conditions such as it has to be a character and also starts with "tRNA-" followed by a valid aminoacid,
  a valid anticodon and additional digits.

```{r}
gene5 <- setFullName(gene5, "tRNA-METAUU 1-1")
gene5
```

### - getDescription

```{r}
getDescription(gene1)
```

### - setDescription

```{r}
gene5 <- setDescription(gene5, "This gene has no specific function")
gene5
```

### - getStructure

```{r}
getStructure(gene5)
```

### - setStructure

```{r}
gene5 <- setStructure(gene5, GRanges("chr1:100-200"))
gene5
```

### - getRRNA

  This is a method of the class "tRNAgene" class. It allows the user to acccess the information about the tRNA codified by the gene.
  This object will contain the information about the ID, sequence of the tRNA and the aminoacid and anticodon carried by it.

```{r}
getTRNA(gene5)
```

###  - setRRNA

  This is a method of the class "tRNAgene" class. It allows the user to acccess the information about the tRNA codified by the gene.
  This object will contain the information about the ID, sequence of the tRNA and the aminoacid and anticodon carried by it.

```{r}
gene5 <- setTRNA(gene5, tRNA(id = "ENST00000000001", tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), aminoacid = AAString("M"), anticodon = RNAString("AGC")))
gene5
```

### - lengthProductTRG

  This function helps in calculating the length of the RNA sequence of the tRNA that is codified by the gene.

```{r}
lengthProductTRG(gene5)
```

## Contact

For any questions or further information, please contact:

**Giulio Vidotto**  
Email: [vidottogiulio@libero.it]
