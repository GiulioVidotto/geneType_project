#' product class
#'
#' This class is a virtual class for all the products of the genes specified
#' in each class. In this virtual class the only attribute specified is the
#' EnsemblID of the product of the genes. This is because it is the only
#' attribute common to all the products (proteis, microRNAs, lncRNAs, rRNAs,
#' tRNAs).
#'
#' @param id Ensembl ID of the product of a gene
#' @export

methods::setClass(
    "Product",
    methods::representation(id = "character"),
    "VIRTUAL"
)


#' protein Product Class
#'
#' This class is called "Protein" and it inherits from the virtual class
#' "product". It represents the product of the protein-coding genes.
#'
#' @param id Ensembl ID of the protein
#' @param protein_sequence Sequence of the protein
#' @param protein_function Description of the function of the protein
#' @examples
#' library(Biostrings)
#' protein(
#'      id = "ENSP01234567891",
#'      protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'      description = "It has a specific function")
#' @importFrom Biostrings AAString
#' @export protein

protein <- methods::setClass(
    "Protein",
    contains = "Product",
    slots = list(
        protein_sequence = "AAString",
        description = "character"
    )
)

methods::setValidity("Protein", function(object) {

    # The EnsemblID of the products of a protein coding gene must start with
    # ENS, followed by P (protein) and then 11 digits
    if (!grepl("^ENSP\\d{11}$", object@id)) {
        stop("The protein enemblID is not correct")
    }

    # Check if the aminoacids taken as input are valid
    if (!all(unlist(strsplit(as.character(object@protein_sequence), "")) %in% AA_ALPHABET)) {
        stop("In the protein sequence there are one or more invalid aminoacids")
    }

    return(TRUE)

})

#' longNonCodingRNA class
#'
#' This class is called "LongNonCodingRNA" and it inherits from the virtual
#' class "product". It represents the product of the Long Non Coding RNA genes.
#'
#' @param id Ensembl ID of the long non coding RNA
#' @param long_non_coding_RNA_sequence Sequence of the long non coding RNA
#' @param regulatory_mechanism Description of the regulatory mechanism of the long non coding RNA
#' @examples
#' library(Biostrings)
#' longNonCodingRNA(
#'      id = "ENST01234567891",
#'      long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'      regulatory_mechanism = "Description of the regulatory mechanism")
#' @importFrom Biostrings RNAString
#' @export longNonCodingRNA

longNonCodingRNA <- methods::setClass(
    "LongNonCodingRNA",
    contains = "Product",
    slots = list(
        long_non_coding_RNA_sequence = "RNAString",
        regulatory_mechanism = "character"
    )
)

methods::setValidity("LongNonCodingRNA", function(object) {

    # The EnsemblID of the long non coding RNA must start with ENS, followed
    # by T (transcript) and then 11 digits
    if (!grepl("^ENST\\d{11}$", object@id)) {
        stop("The long non codingRNA enemblID is not correct")
    }

    return(TRUE)

})

#' microRNA class
#'
#' This class is called "MicroRNA" and it inherits from the virtual class
#' "product". It represents the product of the microRNA genes.
#'
#' @param id Ensembl ID of the MicroRNA
#' @param microRNA_sequence Sequence of the MicroRNA
#' @param silencing_mechanism Description of the silencing mechanism of the MicroRNA
#' @examples
#' library(Biostrings)
#' microRNA(
#'      id = "ENST01234567891",
#'      microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
#'      silencing_mechanism = "The microRNA has ..")
#' @importFrom Biostrings RNAString
#' @export microRNA

microRNA <- setClass(
    "microRNA",
    contains = "Product",
    slots = list(
        microRNA_sequence = "RNAString",
        silencing_mechanism = "character"
    )
)

methods::setValidity("microRNA", function(object) {

    # The EnsemblID of the microRNA must start with ENS, followed by T
    # (transcript) and then 11 digits
    if (!grepl("^ENST\\d{11}$", object@id)) {
        stop("The microRNA enemblID is not correct")
    }

    return(TRUE)

})

#' tRNA class
#' 
#' This class is called "tRNA" and it inherits from the virtual class
#' "product". It represents the product of the tRNA genes.
#' 
#' @param id Ensembl ID of the tRNA
#' @param tRNA_sequence Sequence of the tRNA
#' @param aminoacid The aminoacid carried by the tRNA
#' @param anticodon The anticodon on the tRNA
#' @examples
#' library(Biostrings)
#' tRNA(
#'      id = "ENST00000000001",
#'      tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
#'      aminoacid = AAString("M"),
#'      anticodon = RNAString("AGC"))
#' @importFrom Biostrings RNAString
#' @importFrom Biostrings AAString
#' @export tRNA

tRNA <- methods::setClass(
    "tRNA",
    contains = "Product",
    slots = list(
        tRNA_sequence = "RNAString",
        aminoacid = "AAString",
        anticodon = "RNAString"
    )
)

methods::setValidity("tRNA", function(object) {

    # The EnsemblID of the tRNA must start with ENS,
    # followed by T (transcript) and then 11 digits
    if (!grepl("^ENST\\d{11}$", object@id)) {
        stop("The tRNA enemblID is not correct")
    }
    # Constraint for the aminoacid
    if (!(any(grepl(unlist(object@aminoacid), Biostrings::AA_ALPHABET))) | (length(object@aminoacid) != 1)) {
        stop("Invalid aminoacid")
        }

    # The anticodon must be of length 3 bases
    if (length(unlist(object@anticodon)) != 3) {
        stop("Invalid anticodon")
        }

    return(TRUE)

})

#' rRNA class
#'
#' This class is called "rRNA" and it inherits from the virtual
#' class "product". It represents the product of the rRNA genes.
#'
#' @param id Ensembl ID of the rRNA
#' @param microRNA_sequence Sequence of the rRNA
#' @param ribosomial_subunit The ribosomal subunit
#' @examples
#' library(Biostrings)
#' rRNA(
#'      id = "ENST01234567891",
#'      rRNA_sequence = RNAString("AUUUGA"),
#'      ribosomal_subunit = "28S")
#' @importFrom Biostrings RNAString
#' @export rRNA

rRNA <- methods::setClass(
    "rRNA",
    contains = "Product",
    slots = list(
        rRNA_sequence = "RNAString",
        ribosomal_subunit = "character"
    )
)

methods::setValidity("rRNA", function(object) {

    all_ribosomal_subunits <- c("45S", "28S", "18S", "5S", "5-8S")

    # The EnsemblID of the rRNA must start with ENS, followed by T
    # (transcript) and then 11 digits
    if (!grepl("^ENST\\d{11}$", object@id)) {
        stop("The rRNA enemblID is not correct")
    }

    # The subunit has to be valid (45S, 28S, 18S, 5S, 5-8S)
    if (!(object@ribosomal_subunit %in% all_ribosomal_subunits)) {
        stop("The ribosomal subunit is not correct")
    }

    # Only one subunits per rRNA
    if (length(object@ribosomal_subunit) != 1) {
        stop("The ribosomal subunit is not correct")
    }

    return(TRUE)

})