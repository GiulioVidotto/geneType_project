#' tRNA gene class
#'
#' This class rapresents the genes whose transcript(s) are classified as tRNA
#' (transfer RNA). This class inherits all the attributes of the virtual class
#' "gene" and it furthermore specifies the type of the gene by adding
#' information about the tRNA. The "tRNA" attribute is a single object of class
#' "tRNA". This class inherits from the virtual class "product" and it has
#' different parameters such as the "id" of the tRNA, the sequence of the tRNA
#' ("tRNA_sequence"), the aminoacid carried by the tRNA ("aminoacid") and the
#' anticodon sequence ("anticodon").
#'
#' @param geneID Ensembl gene ID
#' @param gene_symbol HUGO gene symbol
#' @param full_gene_name The full name of the gene
#' @param description A brief description of the gene
#' @param structure chromosome, start, stop, etc.
#' @param tRNA All the information about the tRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"),
#'            anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' @importFrom Biostrings RNAString
#' @importFrom Biostrings AAString
#' @importFrom GenomicRanges GRanges
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods callNextMethod
#' @export tRNAgene

tRNAgene <- methods::setClass("tRNAgene", contains = "gene", slots = list(tRNA = "tRNA"))

methods::setValidity("tRNAgene", function(object) {

    # List of all the aminoacid
    aa <- list(A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE",
                        G = "GLY", H = "HIS", I = "ILE", K = "LYS", L = "LEU",
                        M = "MET", N = "ASN", P = "PRO", Q = "GLN", R = "ARG",
                        S = "SER", T = "THR", V = "VAL", W = "TRP", Y = "TYR")

    # All the letters in the HUGO gene symbol must start with TR[one letter amino acid code]‐[anticodon][gtrnadb gene identifier]
    if (!grepl(paste0("^TR", as.character(object@tRNA@aminoacid), "-", as.character(object@tRNA@anticodon)), object@gene_symbol)) {
        stop("Invalid gene symbol")
    }

    # The full gene name should start with "tRNA-" followed by the aminocid symbol, the anticodon and then additional digits
    if (!grepl(paste0("tRNA-",aa[as.character(object@tRNA@aminoacid)],as.character(object@tRNA@anticodon)), object@full_gene_name)) {
        stop("Invalid full gene name symbol")
    }
 
    return(TRUE)

})


#' setSymbol
#'
#' This is a method inherited from the "gene" class. For the class
#' "tRNAgene" there have been added additional constrains
#' regarding the new HUGO gene symbol. Indeed it must be of class
#' "character" and also must have a specif structure.
#'
#' @param object An tRNA gene
#' @param new_gene_symbol A new symbol for a specific gene
#' @return The gene with the new HUGO gene symbol
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"), anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' gene1 <- setSymbol(gene1, "TRM-AUU1-1")
#' @rdname setSymbol
#' @importFrom Biostrings AAString
#' @importFrom Biostrings RNAString
#' @export setSymbol

methods::setMethod(
    f = "setSymbol",
    signature = "tRNAgene",
    definition = function(object, new_gene_symbol) {

        # List of all the aminoacid

        aa <- list(A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE",
                        G = "GLY", H = "HIS", I = "ILE", K = "LYS", L = "LEU",
                        M = "MET", N = "ASN", P = "PRO", Q = "GLN", R = "ARG",
                        S = "SER", T = "THR", V = "VAL", W = "TRP", Y = "TYR")
        
        # Checking the conditions for the new gene symbol, it must have valid aminoacid and valid anticodon

        if ((!is.character(new_gene_symbol)) |
             !(sub("^TR(\\w)-.*", "\\1", new_gene_symbol) %in% Biostrings::AA_ALPHABET) |
             (nchar(sub("^TR.-(.*?)(\\d*-\\d*)", "\\1", "TRM-AUU1-1")) != 3)) {
            stop("Invalid gene symbol")
        }

        else {
            # Additional part of the gene symbol (Ex. "TRA-AGC1-1", the addition part is "1-1")
            additional_part_gene_symbol <- sub("^TR.-.*?(\\d*-\\d*)", "\\1", new_gene_symbol)

            # Taking the aminoacid from the gene symbol (Ex. "TRA-AGC1-1", the addition part is "A")
            new_aminoacid <- sub("^TR(\\w)-.*", "\\1", new_gene_symbol)

            # Taking the anticodon from the gene symbol (Ex. "TRA-AGC1-1", the addition part is "AGC")
            new_anticodon <- sub("^TR.-(.*?)(\\d*-\\d*)", "\\1", new_gene_symbol)

            # Changing the attribute of the gene object
            object@tRNA@aminoacid <- Biostrings::AAString(new_aminoacid)
            object@tRNA@anticodon <- Biostrings::RNAString(new_anticodon)
            object@full_gene_name <- paste0("tRNA-", aa[as.character(new_aminoacid)], as.character(object@tRNA@anticodon), " ", additional_part_gene_symbol)
            return(methods::callNextMethod())
        }
    }
)

#' setFullName
#'
#' This is a method inherited from the "gene" class. For the class
#' "tRNAgene" there have been added additional constrains
#' regarding the new full gene name. Indeed it must be of class
#' "character" and also must have a specific structure.
#' 
#' @param object An tRNA gene
#' @param new_gene_full_name The new full name of the gene
#' @return The gene with the new full name
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"), anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' gene1 <- setFullName(gene1, "tRNA-METAUU 1-1")
#' @rdname setFullName
#' @importFrom Biostrings AAString
#' @importFrom Biostrings RNAString
#' @export setFullName

methods::setMethod(
    f = "setFullName",
    signature = "tRNAgene",
    definition = function(object, new_gene_full_name) {

        # List of all the aminoacid

        aa <- list(A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE",
                        G = "GLY", H = "HIS", I = "ILE", K = "LYS", L = "LEU",
                        M = "MET", N = "ASN", P = "PRO", Q = "GLN", R = "ARG",
                        S = "SER", T = "THR", V = "VAL", W = "TRP", Y = "TYR")

        # Checking the conditions for the new gene symbol, it must have valid aminoacid and valid anticodon

        if ((!is.character(new_gene_full_name)) |
             !(sub("^tRNA-(\\w{3})(\\w{3}) (\\d*-\\d*)", "\\1", new_gene_full_name) %in% aa) |
             (nchar(sub("^tRNA-(\\w{3})(\\w{3}) (\\d*-\\d*)", "\\2", new_gene_full_name)) != 3)) {
            stop("Invalid full gene name")
        } 
        else {
            # Additional part of the gene full name (Ex. "tRNA-METAUU 1-1", the addition part is "1-1")
            additional_part_gene_symbol <- sub("^tRNA-.* (\\d*-\\d*)", "\\1", new_gene_full_name)

            # Taking the aminoacid from the gene full name (Ex. "tRNA-METAUU 1-1", the addition part is "MET")
            extended_aminoacid <- sub("^tRNA-(\\w{3})(\\w{3}) (\\d*-\\d*)", "\\1", new_gene_full_name)
            # Searching which element of the list aa is equal to the extended_aminoacid (Ex. "MET")
            index_aa <- vapply(aa, function(x) {x == extended_aminoacid}, logical(1))
            # Finding the new aminoacid which corrispond to the name of the element of the list aa with its content equals to extended_aminoacid
            new_aminoacid <- names(which(index_aa))
            
            # Taking the anticodon from the gene symbol (Ex. "TRA-AGC1-1", the addition part is "AGC")
            new_anticodon <- sub("^tRNA-(\\w{3})(\\w{3}) (\\d*-\\d*)", "\\2", new_gene_full_name)

            # Changing the attribute of the gene object
            object@tRNA@aminoacid <- Biostrings::AAString(new_aminoacid)
            object@tRNA@anticodon <- Biostrings::RNAString(new_anticodon)
            object@gene_symbol <- paste0("TR", new_aminoacid, "-", new_anticodon, additional_part_gene_symbol)
            return(methods::callNextMethod())
        }
    }
)


#' setAminoacid
#'
#' This is a method of the class "tRNAgene". The new aminoacid has to respect
#' different conditions, for example it has to be of class AAString and also
#' has to be a valid aminoacid. If this is not the case, then the following
#' will raise: "Invalid aminoacid".
#'
#' @param object A tRNA gene gene
#' @param new_aminoacid A aminoacid for the tRNA gene
#' @return the gene with the tRNA with a new aminoacid
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"), anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' gene1 <- setAminoacid(gene1, AAString("M"))
#' @rdname setAminoacid
#' @export setAminoacid

methods::setGeneric("setAminoacid", function(object, new_aminoacid) {
    standardGeneric("setAminoacid")
})

#' @rdname setAminoacid
methods::setMethod(
    f = "setAminoacid",
    signature = "tRNAgene",
    definition = function(object, new_aminoacid) {
        if ((!is(new_aminoacid, "AAString"))) {
            stop("Invalid aminoacid")}
        else {
            aa <- list(A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE",
                        G = "GLY", H = "HIS", I = "ILE", K = "LYS", L = "LEU",
                        M = "MET", N = "ASN", P = "PRO", Q = "GLN", R = "ARG",
                        S = "SER", T = "THR", V = "VAL", W = "TRP", Y = "TYR")
            # Extracting the additional part of the gene symbol (usually digits)
            additional_part_gene_symbol <- sub(paste0("^.*-", object@tRNA@anticodon), "", object@gene_symbol)
            # new aminoacid
            object@tRNA@aminoacid <- new_aminoacid
            # Creating the new gene symbol, with the new aminoacid
            object@gene_symbol <- paste0("TR", as.character(new_aminoacid), "-", as.character(object@tRNA@anticodon), additional_part_gene_symbol)
            # Creating the new full gene name, with the new aminoacid
            object@full_gene_name <- paste0("tRNA-", aa[as.character(new_aminoacid)], as.character(object@tRNA@anticodon), " ", additional_part_gene_symbol)
            return(object)
        }
    }
)


#' setAnticodon
#'
#' This is a method of the class "tRNAgene". The new anticodon has to respect
#' different conditions, for example it has to be of class RNAString and also
#' its length has to be of 3 base. If this is not the case, then the following
#' will raise: "Invalid aminoacid".
#'
#' @param object A tRNA gene gene
#' @param new_anticodon An anticodon for the tRNA gene
#' @return the gene with the tRNA with a new anticodon
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"),
#'            anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' gene1 <- setAnticodon(gene1, RNAString("AUU"))
#' @rdname setAnticodon
#' @export setAnticodon

methods::setGeneric("setAnticodon", function(object, new_anticodon) {
    standardGeneric("setAnticodon")
})

#' @rdname setAnticodon
methods::setMethod(
    f = "setAnticodon",
    signature = "tRNAgene",
    definition = function(object, new_anticodon) {
        if (!is(new_anticodon, "RNAString") | length(unlist(new_anticodon)) != 3) {
                stop("Invalid anticodon")}
        else {
            aa <- list(A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE",
                        G = "GLY", H = "HIS", I = "ILE", K = "LYS", L = "LEU",
                        M = "MET", N = "ASN", P = "PRO", Q = "GLN", R = "ARG",
                        S = "SER", T = "THR", V = "VAL", W = "TRP", Y = "TYR")
            # Extracting the additional part of the gene symbol (usually digits)
            additional_part_gene_symbol <- sub(paste0("^.*-", as.character(object@tRNA@anticodon)), "", object@gene_symbol) 
            # new anticodon
            object@tRNA@anticodon <- new_anticodon
            # Creating the new gene symbol, with the new anticodon
            object@gene_symbol <- paste0("TR", as.character(object@tRNA@aminoacid), "-", as.character(new_anticodon), additional_part_gene_symbol)
            # Creating the new full gene name, with the new anticodon
            object@full_gene_name <- paste0("tRNA-", aa[as.character(object@tRNA@aminoacid)], as.character(new_anticodon), " ", additional_part_gene_symbol)
            return(object)
        }

})

#' getTRNA
#'
#' This is a method of the class "tRNAgene" class. It allows the user to acccess
#' the information about the tRNA codified by the gene. This object will contain
#' the information about the ID, sequence of the tRNA and the aminoacid and
#' anticodon carried by it.
#'
#' @param object A gene
#' @return The information about the tRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"),
#'            anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' getTRNA(gene1)
#' @rdname getTRNA
#' @export getTRNA
 
methods::setGeneric("getTRNA", function(object) {
    standardGeneric("getTRNA")
})

#' @rdname getTRNA
methods::setMethod(
    f = "getTRNA",
    signature = "tRNAgene",
    definition = function(object) {
        return(object@tRNA)
    }
)

#' setTRNA
#'
#' This is a method of the class "tRNAgene" class. It allows the user to change
#' the information about the tRNA. The user can change the information about the
#' ID, sequence, aminoacid and anticodon of the tRNA.
#'
#' @param object A gene
#' @param new_tRNA The new tRNA
#' @return The gene with the new tRNA
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"),
#'            anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' gene1 <- setTRNA(gene1, tRNA(id = "ENST00000000001",
#'      tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'      aminoacid = AAString("M"),
#'      anticodon = RNAString("AGC")))
#' @rdname setTRNA
#' @export setTRNA
 
methods::setGeneric("setTRNA", function(object, new_tRNA) {
    standardGeneric("setTRNA")
})

#' @rdname setTRNA
methods::setMethod(
    f = "setTRNA",
    signature = "tRNAgene",
    definition = function(object, new_tRNA) {

        if(!is(new_tRNA, "tRNA")) {
            stop("Invalid tRNA")
        }
        else {
            object@tRNA <- new_tRNA
            aa <- list(A = "ALA", C = "CYS", D = "ASP", E = "GLU", F = "PHE",
                        G = "GLY", H = "HIS", I = "ILE", K = "LYS", L = "LEU",
                        M = "MET", N = "ASN", P = "PRO", Q = "GLN", R = "ARG",
                        S = "SER", T = "THR", V = "VAL", W = "TRP", Y = "TYR")
            # Extracting the additional part of the gene symbol (usually digits)
            additional_part_gene_symbol <- sub(paste0("^.*-", object@tRNA@anticodon), "", object@gene_symbol)  
            # Creating the new gene symbol
            object@gene_symbol <- paste0("TR", as.character(object@tRNA@aminoacid), "-", as.character(object@tRNA@anticodon), additional_part_gene_symbol)
            # Creating the new fuull gene name
            object@full_gene_name <- paste0("tRNA-", aa[as.character(object@tRNA@aminoacid)], as.character(object@tRNA@anticodon), " ", additional_part_gene_symbol)
            return(object)
        }
    }
)

#' lengthProductTRG
#'
#' This function helps in calculating the length of the RNA sequence
#' of all the tRNA that are codified by the gene
#'
#' @param object A gene
#' @return The length of the microRNA sequence as the product(s) of the gene
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- tRNAgene(
#'        tRNA = tRNA(
#'            id = "ENST01987654321",
#'            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'            aminoacid = AAString("A"),
#'            anticodon = RNAString("AGC")),
#'        geneID = "ENSG01234567891",
#'        gene_symbol = "TRA-AGC1-1",
#'        full_gene_name = "tRNA-ALAAGC 1-1",
#'        description = "This gene has a specific function",
#'        structure = GRanges("chr1:1-1000"))
#' lengthProductTRG(gene1)
#' @export lengthProductTRG

lengthProductTRG <- function(object) {

                        if (!is(object, "tRNAgene")) {
                            stop("Invalid function object. It should be of class 'tRNAgene'")}

                        else {
                            length_product <- length(object@tRNA@tRNA_sequence)
                            names(length_product) <- object@tRNA@id
                            return (length_product)
                        }
                    }