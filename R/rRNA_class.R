#' rRNA gene class
#'
#' This class rapresents the genes whose transcript(s) are classified as rRNA
#' (ribosomal RNA). This class inherits all the attributes of the virtual class
#' "gene" and it furthermore specifies the type of the gene by adding
#' information about the rRNA. The "rRNA" attribute is a single object of class
#' "rRNA". This class inherits from the virtual class "product" and it has
#' different parameters such as the "id" of the rRNA, the sequence of the rRNA
#' ("rRNA_sequence") and lastly its ribosomal subunit ("ribosomal_subunit").
#'
#' @param geneID Ensembl gene ID
#' @param gene_symbol HUGO gene symbol
#' @param full_gene_name The full name of the gene
#' @param description A brief description of the gene
#' @param structure chromosome, start, stop, etc.
#' @param rRNA All the information about the rRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' @importFrom Biostrings RNAString
#' @importFrom GenomicRanges GRanges
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods callNextMethod
#' @export rRNAgene

rRNAgene <- methods::setClass("rRNAgene", contains = "gene", slots = list(rRNA = "rRNA"))

methods::setValidity("rRNAgene", function(object) {

    # All the letters in the HUGO gene symbol must start with RNA and it muts contain the specific ribosomal subunit
    if (!grepl(paste0("^RNA", object@rRNA@ribosomal_subunit, "$"), object@gene_symbol)) {
        stop("Invalid gene symbol")
    }

    # The full gene name should start with "RNA," followed by the number of the ribosomal subunit and then the word "ribosomal".
    if (!grepl(paste("^RNA,", object@rRNA@ribosomal_subunit, "ribosomal$"), object@full_gene_name)) {
        stop("Invalid full gene name symbol")
    }
 
    return(TRUE)

})

#' setSymbol
#'
#' This is a method inherited from the "gene" class. For the class
#' "rRNAgene" there have been added additional constraints
#' regarding the new HUGO gene symbol. Indeed it must be of class
#' "character" and also must have a specific structure.
#'
#' @param object An rRNA gene
#' @param new_gene_symbol A new symbol for a specific gene
#' @return The gene with the new HUGO gene symbol
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' setRibosomalSub(gene1, "45S")
#' gene1 <- setSymbol(gene1, "RNA45S")
#' @rdname setSymbol
#' @export setSymbol

methods::setMethod(
    f = "setSymbol",
    signature = "rRNAgene",
    definition = function(object, new_gene_symbol) {

        all_ribosomal_subunits <- c("45S", "28S", "18S", "5S", "5-8S")
        if ((!is.character(new_gene_symbol)) | !(sub("^RNA", "", new_gene_symbol) %in% all_ribosomal_subunits)) {
            stop("Invalid gene symbol")
        }
        else {
            new_ribosomal_sub <- sub("^RNA", "", new_gene_symbol)
            object@rRNA@ribosomal_subunit <- new_ribosomal_sub
            object@full_gene_name <- paste("RNA,", new_ribosomal_sub, "ribosomal")
            return(methods::callNextMethod())
        }
    }
)

#' setFullName
#'
#' This is a method inherited from the "gene" class. For the class
#' "rRNAgene" there have been added additional constraints
#' regarding the new full gene name. Indeed it must be of class
#' "character" and also must have a specific structure.
#' 
#' @param object An rRNA gene
#' @param new_gene_full_name The new full name of the gene
#' @return The gene with the new full name
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setFullName(gene1, "RNA, 45S ribosomal")
#' @rdname setFullName
#' @export setFullName

methods::setMethod(
    f = "setFullName",
    signature = "rRNAgene",
    definition = function(object, new_gene_full_name) {
        all_ribosomal_subunits <- c("45S", "28S", "18S", "5S", "5-8S")
        if ((!is.character(new_gene_full_name)) | !(sub("^RNA, (\\d+)S ribosomal$", "\\1S", new_gene_full_name) %in% all_ribosomal_subunits)) {
            stop("Invalid full gene name")
        } 
        else {
            new_ribosomal_sub <- sub("^RNA, (\\d+)S ribosomal$", "\\1S", new_gene_full_name)
            object@rRNA@ribosomal_subunit <- new_ribosomal_sub
            object@full_gene_name <- paste("RNA,", new_ribosomal_sub, "ribosomal")
            return(methods::callNextMethod())
        }
    }
)

#' setRibosomalSub
#'
#' This is a method inherited from the "gene" class. For the class "rRNAgene"
#' the new ribosomal subunit has to respect different conditions, for example
#' it has to be a character object and also be a valid ribosomal subunits.
#' If this is not the case, then the following will raise: "Invalid Ribosomal
#' Subunit".
#'
#' @param object A rRNA gene gene
#' @param new_RibosomalSub A new ribosomal subunit for the rRNA gene
#' @return the gene with the new informations
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setRibosomalSub(gene1, "45S")
#' @rdname setRibosomalSub
#' @export setRibosomalSub

methods::setGeneric("setRibosomalSub", function(object, new_RibosomalSub) {
    standardGeneric("setRibosomalSub")
})

#' @rdname setRibosomalSub
methods::setMethod(
    f = "setRibosomalSub",
    signature = "rRNAgene",
    definition = function(object, new_RibosomalSub) {
        all_ribosomal_subunits <- c("45S", "28S", "18S", "5S", "5-8S")
        if ((!is.character(new_RibosomalSub) | !(new_RibosomalSub %in% all_ribosomal_subunits))) {
            stop("Invalid Ribosomal Subunit")
            }
        else {
            object@rRNA@ribosomal_subunit <- new_RibosomalSub 
            object@gene_symbol <- paste0("RNA", new_RibosomalSub)
            object@full_gene_name <- paste("RNA,", new_RibosomalSub, "ribosomal")
            return(object)
        }
    }
)

#' getRRNA
#'
#' This is a method of the class "rRNAgene" class. It allows the users to
#' extract and inspect the information about the rRNA. This object will
#' contain the information about the ID and ribosomal subunit of the rRNA
#' for a specific gene.
#'
#' @param object A gene
#' @return The information about the rRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' getRRNA(gene1)
#' @rdname getRRNA
#' @export getRRNA
 
methods::setGeneric("getRRNA", function(object) {
    standardGeneric("getRRNA")
})

#' @rdname getRRNA
methods::setMethod(
    f = "getRRNA",
    signature = "rRNAgene",
    definition = function(object) {
        return(object@rRNA)
    }
)

#' setRRNA
#'
#' This method is called "setRRNA" which allows the users to change the
#' information about a specific rRNA. The user can change the information
#' about the ID, sequence and silencing mechanism of the rRNA.
#'
#' @param object A gene
#' @param new_rRNA The new rRNA
#' @return A gene with a new rRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' setRRNA(gene1, rRNA(id = "ENST00000000001",
#'      rRNA_sequence = RNAString("AUUUGA"),
#'      ribosomal_subunit = "45S"))
#' @rdname setRRNA
#' @export setRRNA
 
methods::setGeneric("setRRNA", function(object, new_rRNA) {
    standardGeneric("setRRNA")
})

#' @rdname setRRNA
methods::setMethod(
    f = "setRRNA",
    signature = "rRNAgene",
    definition = function(object, new_rRNA) {

        if(!is(new_rRNA, "rRNA")) {
            stop("Invalid rRNA")
        }
        else {
            object@rRNA <- new_rRNA 
            object@gene_symbol <- paste0("RNA", new_rRNA@ribosomal_subunit)
            object@full_gene_name <- paste("RNA,", new_rRNA@ribosomal_subunit, "ribosomal")
            return(object)
        }
    }
)

#' lengthProductRRG
#'
#' This function helps in calculating the length of the RNA sequence of the
#' rRNA that is codified by the gene.
#'
#' @param gene A gene
#' @return The length of the rRNA sequence
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- rRNAgene(
#'    rRNA = rRNA(
#'          id = "ENST01234567891",
#'          rRNA_sequence = RNAString("AUUUGA"),
#'          ribosomal_subunit = "28S"),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "RNA28S",
#'    full_gene_name = "RNA, 28S ribosomal",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' lengthProductRRG(gene1)
#' @export lengthProductRRG

lengthProductRRG <- function(gene) {

                        if (!is(gene, "rRNAgene")) {
                            stop("Invalid function object. It should be of class 'rRNAgene'")}

                        else {
                            length_product <- length(gene@rRNA@rRNA_sequence)
                            names(length_product) <- gene@rRNA@id
                            return (length_product)
                        }
                    }