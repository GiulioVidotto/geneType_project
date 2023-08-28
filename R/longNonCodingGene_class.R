#' Long non coding RNA gene class
#'
#' This class represents the genes whose transcript(s) are classified as long
#' non coding RNA. This class inherits all the attributes of the virtual class
#' "gene" and it furthermore specifies the type of the gene by adding
#' information about the long non coding RNA(s). The "LongNonCodingRNAs"
#' attribute is a list of elements of the class "LongNonCodingRNA". This
#' class inherits from the virtual class "product" and it has different
#' parameters such as the "id" of the long non coding RNA, the sequence of the
#' long non coding RNA ("long_non_coding_RNA_sequence") and lastly its
#' regulatory mechanism ("regulatory_mechanism").
#'
#' @param geneID Ensembl gene ID
#' @param gene_symbol HUGO gene symbol
#' @param full_gene_name The full name of the gene
#' @param description A brief description of the gene
#' @param structure chromosome, start, stop, ect.
#' @param LongNonCodingRNAs Information about the long non coding RNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' longNonCodingRNAGene(
#'    LongNonCodingRNAs = list(
#'        longNonCodingRNA(
#'              id = "ENST01234567891",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"),
#'        longNonCodingRNA(
#'              id = "ENST19876543210",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "LINC01018",
#'    full_gene_name = "long intergenic non-protein coding RNA 01018",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' @importFrom Biostrings RNAString
#' @importFrom GenomicRanges GRanges
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods callNextMethod
#' @export longNonCodingRNAGene


longNonCodingRNAGene <- methods::setClass(
    "longNonCodingRNAGene",
    contains = "gene",
    slots = list(LongNonCodingRNAs = "list")
)

methods::setValidity("longNonCodingRNAGene", function(object) {

    # Checks on the HUGO gene symbol, depending on the type of long non coding RNA
    if (!grepl("^LINC\\d{5}$|^[A-Z]*-AS\\d*|^[A-Z]*\\d*-DT$|^[A-Z]*-IT\\d*$|^[A-Z]*\\d*-OT\\d*$", object@gene_symbol)) {
        stop("Invalid gene symbol")
    }
    # Checks for the full gene name for long intergenic non-protein coding RNA
    if(grepl("^LINC\\d{5}$", object@gene_symbol)) {
        if(!grepl(paste("^long intergenic non-protein coding RNA", sub("LINC", "", object@gene_symbol)), object@full_gene_name)) {
            stop("The full gene name of LncRNAgenes that are intergenic with respect to protein coding genes should start with 'long intergenic non-protein coding RNA'")
            }
        }
    # Checks for the full gene name for LncRNAs that are antisense to the genomic span of a protein coding gene
    if(grepl("^[A-Z]*-AS", object@gene_symbol)) {
        if(!grepl(paste(sub("-AS\\d*", "", object@gene_symbol), "antisense RNA", sub("[A-Z]*-AS", "", object@gene_symbol)), object@full_gene_name)) {
            stop("The full gene name of LncRNAs that are antisense to the genomic span of a protein coding gene should contain [first part of the gene symbol]'antisense RNA'[possible digits]")
            }
        }
    # Checks for the full gene name for LncRNAs that are divergent to (share a bidirectional promoter with) a protein coding gene
    if(grepl("^[A-Z]*\\d*-DT$", object@gene_symbol)) {
        if(!grepl(paste(sub("-DT", "", object@gene_symbol), "divergent transcript"), object@full_gene_name)) {
            stop("The full gene name of LncRNAs that are divergent to (share a bidirectional promoter with) a protein coding gene should contain [first part of the gene symbol]'divergent transcript'")
            }
        }
    # Checks for the full gene name for LncRNAs that are contained within an intron of a protein coding gene on the same strand
    if(grepl("^[A-Z]*-IT\\d*$", object@gene_symbol)) {
        if(!grepl(paste(sub("-IT\\d*", "", object@gene_symbol), "intronic transcript", sub("[A-Z]*-IT", "", object@gene_symbol)), object@full_gene_name)) {
            stop("The full gene name of LncRNAs that are contained within an intron of a protein coding gene on the same strand should contain [first part of the gene symbol]'intronic transcript'[possible digits]")
            }
        }
    # Checks for the full gene name for LncRNAs that overlap a protein coding gene on the same strand
    if(grepl("^[A-Z]*\\d*-OT\\d*$", object@gene_symbol)) {
        if(!grepl(paste(sub("-OT\\d*", "", object@gene_symbol), "3' UTR overlapping transcript", sub("[A-Z]*\\d*-OT", "", object@gene_symbol)), object@full_gene_name)) {
            stop("The full gene name of LLncRNAs that overlap a protein coding gene on the same strand should contain [first part of the gene symbol]'3' UTR overlapping transcript'[possible digits]")
            }
        }
    # Elements of the list LongNonCodingRNAs must be of class LongNonCodingRNA
    if (!all(vapply(object@LongNonCodingRNAs, function(x) {is(x, "LongNonCodingRNA")}, logical(1)))) {
        stop("The elements of the LongNonCodingRNAs list must be of class 'LongNonCodingRNA'")
    }
    # No replicates
    if (length(unique(vapply(object@LongNonCodingRNAs, function(x) {x@id}, character(1)))) != length(vapply(object@LongNonCodingRNAs, function(x) {x@id}, character(1)))) {
        stop("No replicates in LongNonCodingRNAs")}

    return(TRUE) 

})


#' setSymbol
#'
#' This method of class "longNonCodingRNAGene" which allows
#' the users to change the information about the
#' HUGO gene symbol.
#'
#' @param object A long Non Coding RNA gene
#' @param new_gene_symbol A new symbol for the long Non Coding RNA gene
#' @return The gene with the new HUGO gene symbol
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- longNonCodingRNAGene(
#'    LongNonCodingRNAs = list(
#'        longNonCodingRNA(
#'              id = "ENST01234567891",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"),
#'        longNonCodingRNA(
#'              id = "ENST19876543210",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "LINC01018",
#'    full_gene_name = "long intergenic non-protein coding RNA 01018",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setSymbol(gene1, "LINC01010")
#' @rdname setSymbol
#' @export setSymbol

methods::setMethod(
    f = "setSymbol",
    signature = "longNonCodingRNAGene",
    definition = function(object, new_gene_symbol) {
        if ((!is.character(new_gene_symbol)) | (!grepl("^LINC\\d{5}$|^[A-Z]*-AS\\d*|^[A-Z]*\\d*-DT$|^[A-Z]*-IT\\d*$|^[A-Z]*\\d*-OT\\d*$", new_gene_symbol))) {
            stop("Invalid gene symbol")
        }

        else{
            # Check for the symbol of long intergenic non-protein coding RNA  genes
            if(grepl("^LINC\\d{5}$", new_gene_symbol)){
                object@full_gene_name <- paste("long intergenic non-protein coding RNA", sub("LINC", "", new_gene_symbol))
                return(return(methods::callNextMethod()))
            }
            # Check for the symbol of LncRNAs that are antisense to the genomic span of a protein coding gene
            if(grepl("^[A-Z]*-AS", new_gene_symbol)){
                object@full_gene_name <- paste(sub("-AS\\d*", "", new_gene_symbol), "antisense RNA", sub("[A-Z]*-AS", "", new_gene_symbol))
                return(return(methods::callNextMethod()))
            }
            # Check for the symbol of LncRNAs that are divergent to (share a bidirectional promoter with) a protein coding gene
            if(grepl("^[A-Z]*\\d*-DT$", new_gene_symbol)){
                object@full_gene_name <- paste(sub("-DT", "", new_gene_symbol), "divergent transcript")
                return(return(methods::callNextMethod()))
            }
            # Check for the symbol of LncRNAs that are contained within an intron of a protein coding gene on the same strand
            if(grepl("^[A-Z]*-IT\\d*$", new_gene_symbol)){
                object@full_gene_name <- paste(sub("-IT\\d*", "", new_gene_symbol), "intronic transcript", sub("[A-Z]*-IT", "", new_gene_symbol))
                return(return(methods::callNextMethod()))
            }
            # Check for the symbol of LncRNAs that overlap a protein coding gene on the same strand
            if(grepl("^[A-Z]*\\d*-OT\\d*$", new_gene_symbol)){
                object@full_gene_name <- paste(sub("-OT\\d*", "", new_gene_symbol), "3' UTR overlapping transcript", sub("[A-Z]*\\d*-OT", "", new_gene_symbol))
                return(return(methods::callNextMethod()))
            }
        }
    }
)


#' setFullName
#'
#' This method of class "longNonCodingRNAGene" which allows
#' the users to change the information about the
#' full name of the gene.
#' 
#' @param object A long Non Coding RNA gene
#' @param new_gene_full_name A new full name for the long Non Coding RNA gene
#' @return The gene with the new full gene name
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- longNonCodingRNAGene(
#'    LongNonCodingRNAs = list(
#'        longNonCodingRNA(
#'              id = "ENST01234567891",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"),
#'        longNonCodingRNA(
#'              id = "ENST19876543210",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "LINC01018",
#'    full_gene_name = "long intergenic non-protein coding RNA 01018",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' gene1 <- setFullName(gene1, "long intergenic non-protein coding RNA 01234")
#' @rdname setFullName
#' @export setFullName

methods::setMethod(
    f = "setFullName",
    signature = "longNonCodingRNAGene",
    definition = function(object, new_gene_full_name) {

        # Check for the full gene name of long intergenic non-protein coding RNA genes
        if(grepl("^long intergenic non-protein coding RNA \\d{5}$", new_gene_full_name)){
            addtional_part <- sub(".+ (\\d{5})", "\\1", new_gene_full_name)
            object@gene_symbol <- paste0("LINC", addtional_part)
            return(methods::callNextMethod())
        }
        # Check for the full gene name of LncRNAs that are antisense to the genomic span of a protein coding gene
        else if(grepl("^.+ antisense RNA \\d*$", new_gene_full_name)){
            addtional_part_1 <- sub("(.+) antisense RNA \\d*", "\\1", new_gene_full_name)
            addtional_part_2 <- sub(".+ antisense RNA (\\d*)", "\\1", new_gene_full_name)
            object@gene_symbol <- paste0(addtional_part_1,"-AS", addtional_part_2)
            return(methods::callNextMethod())
        }
        # Check for the full gene name of LncRNAs that are divergent to (share a bidirectional promoter with) a protein coding gene
        else if (grepl("^.+ divergent transcript$", new_gene_full_name)){
            addtional_part <- sub(" divergent transcript", "", new_gene_full_name)
            object@gene_symbol <- paste0(addtional_part, "-DT")
            return(methods::callNextMethod())
        }
        # Check for the full gene name of LncRNAs that are contained within an intron of a protein coding gene on the same strand
        else if (grepl("^.+ intronic transcript \\d*$", new_gene_full_name)){
            addtional_part_1 <- sub("(.*) intronic transcript (\\d*)", "\\1", new_gene_full_name)
            addtional_part_2 <- sub("(.*) intronic transcript (\\d*)", "\\2", new_gene_full_name)
            object@gene_symbol <- paste0(addtional_part_1, "-IT", addtional_part_2)
            return(methods::callNextMethod())
        }
        # Check for the full gene name of LncRNAs that overlap a protein coding gene on the same strand
        else if (grepl("^.+ 3' UTR overlapping transcript \\d*$", new_gene_full_name)){
            addtional_part_1 <- sub("(.*) 3' UTR overlapping transcript (\\d*)", "\\1", new_gene_full_name)
            addtional_part_2 <- sub("(.*) 3' UTR overlapping transcript (\\d*)", "\\2", new_gene_full_name)
            object@gene_symbol <- paste0(addtional_part_1, "-OT", addtional_part_2)
            return(methods::callNextMethod())
        }

        else{
            stop("Invalid full gene name")
        }

        }
)


#' getLongNonCodingRNA
#
#' This method of class "longNonCodingRNAGene" which allows the users to
#' extract and inspect the information about the long Non Coding RNAs.
#' This object will contain the information about the ID, sequence and
#' regulatory mechanism of each long non coding RNA defined for a specific gene.
#'
#' @param object A gene
#' @return The information about of all the possible long non coding RNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- longNonCodingRNAGene(
#'    LongNonCodingRNAs = list(
#'        longNonCodingRNA(
#'              id = "ENST01234567891",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"),
#'        longNonCodingRNA(
#'              id = "ENST19876543210",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "LINC01018",
#'    full_gene_name = "long intergenic non-protein coding RNA 01018",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' getLongNonCodingRNA(gene1)
#' @rdname getLongNonCodingRNA
#' @export getLongNonCodingRNA
 
methods::setGeneric("getLongNonCodingRNA", function(object) {
    standardGeneric("getLongNonCodingRNA")
})

#' @rdname getLongNonCodingRNA
methods::setMethod(f = "getLongNonCodingRNA",
        signature = "longNonCodingRNAGene",
        definition = function(object) {

            return(object@LongNonCodingRNAs)
        })

#' setlongNonCodingRNA
#'
#' This method of class "longNonCodingRNAGene" which allows the users
#' to change the information about the long Non Coding RNAs.
#'
#' @param object A gene
#' @param lncRNA_EnsemblID The EnsemblID of the long non coding RNA that the user wants to change
#' @param new_lncRNA The new long non coidng RNA
#' @return The gene with the new long non coding RNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- longNonCodingRNAGene(
#'    LongNonCodingRNAs = list(
#'        longNonCodingRNA(
#'              id = "ENST01234567891",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"),
#'        longNonCodingRNA(
#'              id = "ENST19876543210",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "LINC01018",
#'    full_gene_name = "long intergenic non-protein coding RNA 01018",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' setlongNonCodingRNA(
#'          gene1, "ENST01234567891",
#'          longNonCodingRNA(
#'              id = "ENST00000000001",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"))
#' @rdname setlongNonCodingRNA
#' @export setlongNonCodingRNA
 
methods::setGeneric("setlongNonCodingRNA", function(object, lncRNA_EnsemblID, new_lncRNA) {
    standardGeneric("setlongNonCodingRNA")
})

#' @rdname setlongNonCodingRNA
methods::setMethod(
    f = "setlongNonCodingRNA",
    signature = "longNonCodingRNAGene",
    definition = function(object, lncRNA_EnsemblID, new_lncRNA) {

        if(!is(new_lncRNA, "LongNonCodingRNA")) {
            stop("Invalid longNonCodingRNA")
        }

        lncRNA_index <- grep(lncRNA_EnsemblID, object@LongNonCodingRNAs)
        object@LongNonCodingRNAs[[lncRNA_index]] <- new_lncRNA

        return(object)
    }
)

#' lengthProductLNCRG
#'
#' This function helps in calculating the length of the RNA sequence
#' of all the long non coding RNA that are codified by the gene
#'
#' @param object A gene
#' @return The length of the long non coding RNA sequence
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- longNonCodingRNAGene(
#'    LongNonCodingRNAs = list(
#'        longNonCodingRNA(
#'              id = "ENST01234567891",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism"),
#'        longNonCodingRNA(
#'              id = "ENST19876543210",
#'              long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'              regulatory_mechanism = "Description of the specific regulatory mechanism")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "LINC01018",
#'    full_gene_name = "long intergenic non-protein coding RNA 01018",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000"))
#' lengthProductLNCRG(gene1)
#' @export lengthProductLNCRG

lengthProductLNCRG <- function(object) {

                        if (!is(object, "longNonCodingRNAGene")) {
                            stop("Invalid function object. It should be of class 'longNonCodingRNAGene'")}

                        else {
                            length_product <- vapply(object@LongNonCodingRNAs, function(x) {length(x@long_non_coding_RNA_sequence)}, integer(1))
                            names(length_product) <- vapply(object@LongNonCodingRNAs, function(x) {x@id}, character(1))
                            return (length_product)
                        }
                    }