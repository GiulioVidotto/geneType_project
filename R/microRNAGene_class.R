#' microRNA gene class
#'
#' This class represents the genes whose transcript(s) are classified as
#' microRNA. This class inherits all the attributes of the virtual class
#' "gene" and it furthermore specifies the type of the gene by adding
#' information about the microRNA(s). The "microRNAs" attribute is a list
#' of elements of the class "microRNA". This class inherits from the virtual
#' class "product" and it has different parameters such as the "id" of the
#' microRNA, the sequence of the microRNA ("microRNA_sequence") and lastly
#' its silencing mechanism ("silencing_mechanism").
#'
#' @param geneID Ensembl gene ID
#' @param gene_symbol HUGO gene symbol
#' @param full_gene_name The full name of the gene
#' @param description A brief description of the gene
#' @param structure chromosome, start, stop, ect.
#' @param microRNAs microRNAs of a specific gene
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' microRNAGene(
#'    microRNAs = list(
#'        microRNA(id = "ENST01234567891",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has .."),
#'        microRNA(id = "ENST19876543210",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has ..")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "MIR 12",
#'    full_gene_name = "microRNA 23",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' @importFrom Biostrings RNAString
#' @importFrom GenomicRanges GRanges
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods callNextMethod
#' @export microRNAGene

microRNAGene <- methods::setClass(
    "MicroRNAGene",
    contains = "gene",
    slots = list(microRNAs = "list")
)

methods::setValidity("MicroRNAGene", function(object) {

    # Letters in the HUGO symbol must start with MIR followed by a number
    if (!grepl("^MIR\\s\\d*$", object@gene_symbol)) {
        stop("Invalid gene symbol")
    }

    # The full gene name must start with microRNA followed by a space and a number
    if (!grepl("^microRNA\\s\\d*$", object@full_gene_name)) {
        stop("Invalid full gene name symbol")
    }

    # Elements of the list microRNAs must be of class microRNA
    if (!all(vapply(object@microRNAs, function(x) {is(x, "microRNA")}, logical(1)))) {
        stop("The elements of the microRNAs list must be of class 'microRNA'")
    }
    
    # No replicates
    if (length(unique(vapply(object@microRNAs, function(x) {x@id}, character(1)))) != length(vapply(object@microRNAs, function(x) {x@id}, character(1)))) {
        stop("No replicates in microRNAs")
    }
 
    return(TRUE) 

})

#' setSymbol
#'
#' This is a method inherited from the "gene" class. For the class
#' "MicroRNAGene" the new gene symbol has to respect different
#' conditions, for example it has to be a character object and also
#' start with "MIR " followed by a number. If this is not the case,
#' then the following will raise: "Invalid gene symbol".
#'
#' @param object A microRNA gene
#' @param new_gene_symbol A new HUGO gene symbol for the microRNA gene
#' @return The microRNA gene with the new HUGO gene symbol
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- microRNAGene(
#'    microRNAs = list(
#'        microRNA(id = "ENST01234567891",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has .."),
#'        microRNA(id = "ENST19876543210",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has ..")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "MIR 12",
#'    full_gene_name = "microRNA 23",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setSymbol(gene1, "MIR 123")
#' @rdname setSymbol
#' @export setSymbol

methods::setMethod(f = "setSymbol",
        signature = "MicroRNAGene",
        definition = function(object, new_gene_symbol) {
            if ((!is.character(new_gene_symbol)) | (!grepl("^MIR\\s\\d*$", new_gene_symbol))) {
                stop("Invalid gene symbol")
            }
            else {
                return(methods::callNextMethod())
        }
    }
)


#' setFullName
#'
#' This is a method inherited from the "gene" class. For the class
#' "MicroRNAGene" the new full name has to respect different conditions
#' such as it has to be a character and also start with "microRNA"
#' followed by a space and then a number.
#'
#' @param object A microRNA gene
#' @param new_gene_full_name The new full name of the microRNA gene
#' @return The microRNA gene with the new full gene name
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- microRNAGene(
#'    microRNAs = list(
#'        microRNA(id = "ENST01234567891",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has .."),
#'        microRNA(id = "ENST19876543210",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has ..")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "MIR 12",
#'    full_gene_name = "microRNA 23",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setFullName(gene1, "microRNA 123")
#' @rdname setFullName
#' @export setFullName


methods::setMethod(f = "setFullName",
        signature = "MicroRNAGene",
        definition = function(object, new_gene_full_name) {
            if ((!is.character(new_gene_full_name)) | (!grepl("^microRNA\\s\\d*$", new_gene_full_name))) {
                stop("Invalid full gene name")
            }

            else {
                return(methods::callNextMethod())
        }
    }
)

#' getMicroRNA
#'
#' This is a method of the class "MicroRNAGene" class. It allows the users to
#' extract and inspect the information about the microRNAs. This object will
#' contain the information about the ID, sequence and silencing mechanism
#' of each microRNA defined for a specific gene.
#'
#' @param object A gene
#' @return The silencing mechanisms of all the possible microRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- microRNAGene(
#'    microRNAs = list(
#'        microRNA(id = "ENST01234567891",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has .."),
#'        microRNA(id = "ENST19876543210",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has ..")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "MIR 12",
#'    full_gene_name = "microRNA 23",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' getMicroRNA(gene1)
#' @rdname getMicroRNA
#' @export getMicroRNA
 
methods::setGeneric("getMicroRNA", function(object) {
    standardGeneric("getMicroRNA")
})

#' @rdname getMicroRNA
methods::setMethod(f = "getMicroRNA",
        signature = "MicroRNAGene",
        definition = function(object) {
            return(object@microRNAs)
        })

#' setMicroRNA
#'
#' This method is called "setMicroRNA" which allows the users to change the
#' information about a specific microRNA. The user can change the information
#' about the ID, sequence and silencing mechanism of the microRNA.
#'
#' @param object A gene
#' @param miRNA_EnsemblID The EnsemblID of the microRNA for which the user wants to chaange the information
#' @param new_microRNA The new microRNA
#' @return The gene with the new microRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- microRNAGene(
#'    microRNAs = list(
#'        microRNA(id = "ENST01234567891",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has .."),
#'        microRNA(id = "ENST19876543210",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has ..")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "MIR 12",
#'    full_gene_name = "microRNA 23",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setMicroRNA(
#'      gene1, "ENST01234567891",
#'      microRNA(
#'          id = "ENST00000000001",
#'          microRNA_sequence = RNAString("ACG"),
#'          silencing_mechanism = "Its silencing mechanism consists of")
#' )
#' @rdname setMicroRNA
#' @export setMicroRNA
 
methods::setGeneric("setMicroRNA", function(object, miRNA_EnsemblID, new_microRNA) {
    standardGeneric("setMicroRNA")
})

#' @rdname setMicroRNA
methods::setMethod(
    f = "setMicroRNA",
    signature = "MicroRNAGene",
    definition = function(object, miRNA_EnsemblID, new_microRNA) {

        if(!is(new_microRNA, "microRNA")) {
            stop("Invalid microRNA")
        }

        miRNA_index <- grep(miRNA_EnsemblID, object@microRNAs)
        object@microRNAs[[miRNA_index]] <- new_microRNA

        return(object)
    }
)

#' lengthProductMRG
#'
#' This function helps in calculating the length of the RNA sequence
#' of all the microRNA that are codified by the gene
#'
#' @param object A gene
#' @return The length of the microRNAs
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- microRNAGene(
#'    microRNAs = list(
#'        microRNA(id = "ENST01234567891",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has .."),
#'        microRNA(id = "ENST19876543210",
#'                  microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'                  silencing_mechanism = "The microRNA has ..")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "MIR 12",
#'    full_gene_name = "microRNA 23",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' lengthProductMRG(gene1)
#' @export lengthProductMRG

lengthProductMRG <- function(object) {

                        if (!is(object, "MicroRNAGene")) {
                            stop("Invalid function object. It should be of class 'MicroRNAGene'")}

                        else {
                            length_product <- vapply(object@microRNAs, function(x) {length(x@microRNA_sequence)}, integer(1))
                            names(length_product) <- vapply(object@microRNAs, function(x) {x@id}, character(1))
                            return (length_product)
                        }
                    }