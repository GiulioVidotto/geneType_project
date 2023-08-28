#' protein coding gene class
#'
#' This class represents the genes whose transcript(s) codify for proteins.
#' This class inherits all the attributes of the virtual class "gene" and
#' it furthermore specifies the type of the gene by adding information about
#' the protein(s), codified by this type of gene. The "proteins" attribute is
#' a list of elements of the class "Protein". This class inherits from the
#' virtual class "product" and it has different parameters such as the "id"
#' of the protein, the sequence of the protein ("protein_sequence") and lastly
#' the description of the protein ("description").
#'
#' @param geneID Ensembl gene ID
#' @param gene_symbol HUGO gene symbol
#' @param full_gene_name The full name of the gene
#' @param description A brief description of the gene
#' @param structure chromosome, start, stop, etc.
#' @param proteins list of all the protein codified by the gene
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function"),
#'        protein(id = "ENSP19876543210",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' @importFrom Biostrings AAString
#' @importFrom GenomicRanges GRanges
#' @importFrom methods new
#' @importFrom methods is
#' @importFrom methods callNextMethod
#' @export proteinCodingGene

proteinCodingGene <- methods::setClass(
    "ProteinCodingGene", 
    contains = "gene", 
    slots = list(proteins = "list")
)

methods::setValidity("ProteinCodingGene", function(object) {

    # All the letters in the HUGO gene symbol must be all uppercase
    if (object@gene_symbol != toupper(object@gene_symbol)) {
        stop("Invalid gene symbol")
    }

    # The full gene name should be all uppercaase (Human genes)
    else if (object@full_gene_name != toupper(object@full_gene_name)) {
        stop("Invalid full gene name symbol")
    }

    # Elements of the list Proteins must be of class protein
    if (!all(vapply(object@proteins, function(x) {is(x, "Protein")}, logical(1)))) {
        stop("The elements of proteins must be of class 'Protein'")
    }
    
    # No replicates
    if (length(unique(vapply(object@proteins, function(x) {x@id}, character(1)))) != length(vapply(object@proteins, function(x) {x@id}, character(1)))) {
        stop("No replicates in proteins")
    }

    return(TRUE) 

})


#' setSymbol
#'
#' This is a method inherited from the "gene" class. For the class
#' "ProteinCodingGene" there have been added additional constrains
#' regarding the new HUGO gene symbol. Indeed it must be of class
#' "character" and also must be uppercase.
#'
#' @param object A protein coding gene
#' @param new_gene_symbol A new symbol for a specific gene
#' @return The gene with the new HUGO gene symbol
#' @rdname setSymbol
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function"),
#'        protein(id = "ENSP19876543210",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' gene1 <- setSymbol(gene1, "GENE1")
#' @export setSymbol

methods::setMethod(
    f = "setSymbol",
    signature = "ProteinCodingGene",
    definition = function(object, new_gene_symbol) {

        if ((!is.character(new_gene_symbol)) | (new_gene_symbol != toupper(new_gene_symbol))) {
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
#' "ProteinCodingGene" there have been added additional constrains
#' regarding the new full gene name. Indeed it must be of class
#' "character" and also must be uppercase.
#' 
#' @param object A protein coding gene
#' @param new_gene_full_name The new full name of the protein coding gene
#' @return The gene with the new full name
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function"),
#'        protein(id = "ENSP19876543210",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' gene1 <- setFullName(gene1, "GENE TWO")
#' @rdname setFullName
#' @export setFullName

methods::setMethod(
    f = "setFullName",
    signature = "ProteinCodingGene",
    definition = function(object, new_gene_full_name) {
        if ((!is.character(new_gene_full_name)) | (new_gene_full_name != toupper(new_gene_full_name))) {
            stop("Invalid full gene name")
        } 
        else {
            return(methods::callNextMethod())
        }
    }
)

#' getProteins
#'
#' This is a method of the class "ProteinCodingGene" class. It allows
#' the user to acccess the information about the protein(s) codified
#' by the gene. This object will contain the information about the ID,
#' sequence and description of each protein defined for a specific gene.
#'
#' @param object A gene
#' @return The protein function of all the possible isoform of a gene
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function"),
#'        protein(id = "ENSP19876543210",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' getProteins(gene1)
#' @rdname getProteins
#' @export getProteins
 
methods::setGeneric("getProteins", function(object) {
    standardGeneric("getProteins")
})

#' @rdname getProteins
methods::setMethod(
    f = "getProteins",
    signature = "ProteinCodingGene",
    definition = function(object) {
        return(object@proteins)
    }
)


#' setProtein
#'
#' This is a method of the class "ProteinCodingGene" class. It allows
#' the user to modify the information about a specific protein codified
#' by the gene.
#'
#' @param object A gene
#' @param protein_EnsemblID EnsemblID of the protein for which the user wants to change the function
#' @param new_protein The new function for the protein
#' @return The protein gene with a new protein
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function"),
#'        protein(id = "ENSP19876543210",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' gene1 <- setProtein(gene1, "ENSP01234567891",
#'      protein(id = "ENSP00000000001",
#'      protein_sequence = AAString("ARNDCQ"),
#'      description = "It has no specific function"))
#' @rdname setProtein
#' @export setProtein
 
methods::setGeneric("setProtein", function(object, protein_EnsemblID, new_protein) {
    standardGeneric("setProtein")
})

#' @rdname setProtein
methods::setMethod(
    f = "setProtein",
    signature = "ProteinCodingGene",
    definition = function(object, protein_EnsemblID, new_protein) {

        if(!is(new_protein, "Protein")) {
            stop("Invalid protein")
        }

        protein_index <- grep(protein_EnsemblID, object@proteins)
        object@proteins[[protein_index]] <- new_protein

        return(object)
    }
)


#' lengthProductPCG
#'
#' This function helps in calculating the length of the aminoacid
#' sequence of all the proteins that are codified by the gene.
#'
#' @param object A gene
#' @return The length of the aminoacid sequence of the product(s) of the gene
#' @examples
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function"),
#'        protein(id = "ENSP19876543210",
#'                  protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'                  description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' lengthProductPCG(gene1)
#' @export lengthProductPCG

lengthProductPCG <- function(object) {

                        if (!is(object, "ProteinCodingGene")) {
                            stop("Invalid function object. It should be of class 'ProteinCodingGene'")}

                        else {
                            length_product <- vapply(object@proteins, function(x) {length(x@protein_sequence)}, integer(1))
                            names(length_product) <- vapply(object@proteins, function(x) {x@id}, character(1))
                            return (length_product)
                        }
                    }