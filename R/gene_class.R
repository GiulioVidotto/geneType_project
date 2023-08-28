#' gene class
#'
#' This class represents the genes in a general way. Indeed, it stores all the
#' information that is common to all the genes. In particular, it defines the
#' Ensmbl ID, the HUGO symbol, the name of the gene, a short description, and
#' all the information about the structure of the gene (e.g. the chromosome,
#' start and stop, etc.). It is important to remember that it is a virtual
#' class, so this means that all the attributes and methods of the class will
#' be inherited by all the specific gene classes. Because, as said before, this
#' class represents the genes in a general way, the only additional constraint
#' is regarding the Ensembl ID of the gene, which must start with "ENSG" and be
#' followed by 11 digits.
#'
#' @param geneID Ensembl ID of the gene
#' @param gene_symbol HUGO gene symbol
#' @param full_gene_name The full name of the gene
#' @param description A brief description of the gene
#' @param structure chromosome, start, stop, ect.
#' @importFrom GenomicRanges GRanges
#' @importFrom methods is
#' @export

methods::setClass("gene", representation(geneID = "character",
                                gene_symbol = "character",
                                full_gene_name = "character",
                                description = "character",
                                structure = "GRanges"),
                                "VIRTUAL")

# Additional constrain about the attribute of the class
methods::setValidity("gene", function(object) {

    # The EnsemblID gene must start with "ENSG" (because it is a human gene)
    # followed by 11 digits
    if (!grepl("^ENSG\\d{11}$", object@geneID)) {
        stop("The gene enemblID is not correct")
        }

    else {
       return(TRUE)
    }
})


#' getGeneId
#' 
#' This is a method of the "gene" class that can be used to access information
#' about the Ensembl ID gene. The get method is the only method specified for
#' the Ensembl ID gene because usually the ID of a gene should not change, even
#' if the gene is updated.
#' 
#' @param object A gene
#' @return The Ensembl ID of the gene
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used.
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function"),
#'        protein(id = "ENSP19876543210", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' # getGeneId(gene1)
#' @rdname getGeneId
#' @export getGeneId

methods::setGeneric("getGeneId", function(object) {
    standardGeneric("getGeneId")
})

#' @rdname getGeneId
methods::setMethod(f = "getGeneId",
        signature = "gene",
        definition = function(object) {
            return(object@geneID)
        })


#' getSymbol
#'
#' This is a method of the "gene" class called "getSymbol". It can be used to
#' access information about the the HUGO gene symbol.
#'
#' @param object A gene
#' @return The new gene object with the new HUGO gene symbol
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used.
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function"),
#'        protein(id = "ENSP19876543210", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' getSymbol(gene1)
#' @rdname getSymbol
#' @export getSymbol

methods::setGeneric("getSymbol", function(object) {
    standardGeneric("getSymbol")
})

#' @rdname getSymbol
methods::setMethod(
    f = "getSymbol",
    signature = "gene",
    definition = function(object) {
        return(object@gene_symbol)
    }
)


#' setSymbol
#'
#' This is a method of the "gene" class called "setSymbol" that can be used to
#' change the information about the HUGO gene symbol.
#'
#' @param object A gene
#' @param new_gene_symbol A new symbol for the gene
#' @return The gene with the new HUGO gene symbol
#' @rdname setSymbol
#' @export setSymbol

methods::setGeneric("setSymbol", function(object, new_gene_symbol) {
    standardGeneric("setSymbol")
})

#' @rdname setSymbol
methods::setMethod(
    f = "setSymbol",
    signature = "gene",
    definition = function(object, new_gene_symbol) {
        object@gene_symbol <- new_gene_symbol
        return(object)
    }
)


#' getFullName
#'
#' This is a method of the "gene" class, called "getFullName", that allows the
#' user to extract the information about the full name of the gene.
#'
#' @param object A gene
#' @return The full name of the gene
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used.
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(id = "ENSP01234567891", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function"),
#'        protein(id = "ENSP19876543210", protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' # getFullName(gene1)
#' @rdname getFullName
#' @export getFullName

methods::setGeneric("getFullName", function(object) {
    standardGeneric("getFullName")
})

#' @rdname getFullName
methods::setMethod(
    f = "getFullName",
    signature = "gene",
    definition = function(object) {
        return(object@full_gene_name)
    }
)

#' setFullName
#'
#' This is a method of the "gene" class, called "setFullName", that allows the
#' user to change the information about the full name of the gene.
#'
#' @param object A gene
#' @param new_gene_full_name The new full name of the gene
#' @return The gene with the new full name
#' @rdname setFullName
#' @export setFullName

methods::setGeneric("setFullName", function(object, new_gene_full_name) {
    standardGeneric("setFullName")
})

#' @rdname setFullName
methods::setMethod(
    f = "setFullName",
    signature = "gene",
    definition = function(object, new_gene_full_name) {
        object@full_gene_name <- new_gene_full_name
        return(object)
    }
)


#' getDescription
#'
#' This is a method of the "gene" class, called "getDescription", and it allows
#' the user to extract the description of a specific gene.
#'
#' @param object A gene
#' @return The description of the gene
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(
#'              id = "ENSP01234567891",
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'              description = "It has a specific function"),
#'        protein(
#'              id = "ENSP19876543210", 
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
#'              description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' getDescription(gene1)
#' @rdname getDescription
#' @export getDescription

methods::setGeneric("getDescription", function(object) {
    standardGeneric("getDescription")
})

#' @rdname getDescription
methods::setMethod(
    f = "getDescription",
    signature = "gene",
    definition = function(object) {
        return(object@description)
    }
)

#' setDescription
#'
#' This is a method of the "gene" class, called "setDescription", and it allows
#' the user to change the description of a specific gene.
#'
#' @param object A gene
#' @param new_description A new description of the gene
#' @return The gene with the new description
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(
#'              id = "ENSP01234567891",
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'              description = "It has a specific function"),
#'        protein(
#'              id = "ENSP19876543210",
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'              description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' setDescription(gene1, "new description")
#' @rdname setDescription
#' @export setDescription

methods::setGeneric("setDescription", function(object, new_description) {
    standardGeneric("setDescription")
})

#' @rdname setDescription
methods::setMethod(
    f = "setDescription",
    signature = "gene",
    definition = function(object, new_description) {
        if (!is.character(new_description)) {
            stop("The description is not correct")}
        else {
            object@description <- new_description
            return(object)
        }
    }
)


#' getStructure
#'
#' This is a method of the "gene" class, called "getStructure", and it allows
#' the user to extract from a specific object all the information on the
#' structure of the gene.
#'
#' @param object A gene
#' @return The information about the structure of the gene
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(
#'              id = "ENSP01234567891", 
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
#'              description = "It has a specific function"),
#'        protein(
#'              id = "ENSP19876543210", 
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
#'              description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' getStructure(gene1)
#' @rdname getStructure
#' @export getStructure

methods::setGeneric("getStructure", function(object) {
    standardGeneric("getStructure")
})

#' @rdname getStructure
methods::setMethod(
    f = "getStructure",
    signature = "gene",
    definition = function(object) {
        return(object@structure)
    }
)


#' setStructure
#'
#' This is a method of the "gene" class, called "setStructure", that is useful
#' to extract from a specific gene all the information about the structure.
#'
#' @param object A gene
#' @param new_structure The new structure of the gene
#' @return The gene with the new structure
#' @examples
#' # Because the gene class is a VIRTUAL class, for the example, an object of
#' # class "proteinCodingGene" has been used
#' library(Biostrings)
#' library(GenomicRanges)
#' gene1 <- proteinCodingGene(
#'    proteins = list(
#'        protein(
#'              id = "ENSP01234567891",
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'              description = "It has a specific function"),
#'        protein(
#'              id = "ENSP19876543210",
#'              protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
#'              description = "It has a specific function")),
#'    geneID = "ENSG01234567891",
#'    gene_symbol = "GENE1",
#'    full_gene_name = "GENE ONE",
#'    description = "This gene has a specific function",
#'    structure = GRanges("chr1:1-1000")
#' )
#' setStructure(gene1, GRanges("chr1:1-100"))
#' @rdname setStructure
#' @export setStructure

methods::setGeneric("setStructure", function(object, new_structure) {
    standardGeneric("setStructure")
})

#' @rdname setStructure
methods::setMethod(
    f = "setStructure",
    signature = "gene",
    definition = function(object, new_structure) {
        if (!is(new_structure, "GRanges")) {
            stop("The new defined structure of the gene is not correct")}
        else {
            object@structure <- new_structure
            return(object)
        }
    }
)
