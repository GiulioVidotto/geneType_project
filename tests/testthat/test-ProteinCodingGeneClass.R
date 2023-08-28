# Testing if a ProteinCodingGene gene object is define in a correct way
test_that("correct ProteinCodingGene object", {

    gene1<- proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891",
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    )

    expect_no_error(gene1)

})

# Tests regarding the validity of object of class ProteinCodingGene

# Testing an object for class ProteinCodingGene with invalid EnsemblID
test_that("invalid EnsemblID of the gene for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG012345",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ),"The gene enemblID is not correct")
})


# Tests regarding the validity of object of class ProteinCodingGene

# Testing an object for class ProteinCodingGene with invalid HUGO gene symbol
test_that("invalid gene symbol for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "gene1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ),"Invalid gene symbol")
})

# Testing an object for class ProteinCodingGene with invalid full gene name
test_that("invalid full gene name for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210",
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "gene one",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ),"Invalid full gene name symbol")
})

# Testing an object for class ProteinCodingGene with invalid proteins
test_that("invalid protein objetc for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            id = "ENSP01234567891", 
            protein_sequence = AAString("123"), 
            description = "It has a specific function"),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The elements of proteins must be of class 'Protein'")

})

# Testing an object for class ProteinCodingGene with invalid EnsemblID of the protein
test_that("invalid EnsemblID of the protein for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            protein(
                id = "ENS123",
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The protein enemblID is not correct")

})

# Testing an object for class ProteinCodingGene with invalid aminoacid sequence
test_that("invalid protein sequence for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("123"), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "In the protein sequence there are one or more invalid aminoacids")

})

# Testing an object for class ProteinCodingGene with a replicate proteins
test_that("replicates in proteins for ProteinCodingGene object", {

    expect_error(proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."),
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "No replicates in proteins")

})

# Tests of all the method of the class ProteinCodingGene
# Before doing the test, an object of class ProteinCodingGene was defined
gene1 <- proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )

# Tests of all the getter methods

# Testing getGeneId method
test_that("method getGeneId executed correctly", {
    expect_equal(getGeneId(gene1), "ENSG01234567891")
})

# Testing getSymbol method
test_that("method getSymbol executed correctly", {
    expect_equal(getSymbol(gene1), "GENE1")
})

# Testing getFullName method
test_that("method getSymbol executed correctly", {
    expect_equal(getFullName(gene1), "GENE ONE")
})

# Testing getDescription method
test_that("method getDescription executed correctly", {
    expect_equal(getDescription(gene1), "This gene has a specific function")
})

# Testing getStructure method
test_that("method getStructure executed correctly", {
    expect_equal(getStructure(gene1), GRanges("chr1:1-1000"))
})

# Testing getProteins method
test_that("method getProduct executed correctly", {
    expect_equal(getProteins(gene1), list(
        protein(
            id = "ENSP01234567891", 
            protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
            description = "It has a specific function"),
        protein(
            id = "ENSP19876543210", 
            protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
            description = "It has a specific function")))
})

# Tests of all the setter methods

# Testing setSymbol method
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "GENE2"),
        proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE2",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with invalid new gene symbol
test_that("method setSymbol with invalid gene symbol", {
    expect_error(
        setSymbol(gene1, "gene1"), "Invalid gene symbol"
        )
})

# Testing setFullName method

test_that("method setFullName executed correctly", {
    expect_equal(
        setFullName(gene1, "GENE TWO"),
        proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE TWO",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setFullName method with invalid new full name gene
test_that("method setFullName with invalid new full name gene", {
    expect_error(
        setFullName(gene1, "gene TWO"), "Invalid full gene name"
        )
})

# Testing setDescription method
test_that("method setDescription executed correctly", {
    expect_equal(setDescription(gene1, "This gene has no specific function"),
    proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has no specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setDescription method with invalid description
test_that("method setDescription with invalid description", {
    expect_error(setDescription(gene1, 123), "The description is not correct"
    )
})

# Testing setStructure method
test_that("method setStructure executed correctly", {
    expect_equal(setStructure(gene1, GRanges("chr1:100-200")),
    proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP01234567891", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:100-200")
        )
    )
})

# Testing setStructure method with invalid new structure
test_that("method setStructure with invalid structure", {
    expect_error(setStructure(gene1, "chr1:100-200"), 
    "The new defined structure of the gene is not correct"
    )
})

# Testing setProtein method for object of class ProteinCodingGene
test_that("method setProtein executed correctly", {
    expect_equal(
        setProtein(
            gene1, 
            "ENSP01234567891", 
            protein(
                id = "ENSP00000000001", 
                protein_sequence = AAString("ARNDCQ"), 
                description = "It has no specific function")),
        proteinCodingGene(
        proteins = list(
            protein(
                id = "ENSP00000000001", 
                protein_sequence = AAString("ARNDCQ"), 
                description = "It has no specific function"),
            protein(
                id = "ENSP19876543210", 
                protein_sequence = AAString("ARNDCQEGHILKMFPSTWYVUOBJZX*-+."), 
                description = "It has a specific function")),
        geneID = "ENSG01234567891",
        gene_symbol = "GENE1",
        full_gene_name = "GENE ONE",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setProtein method with invalid new protein
test_that("method setProtein with invalid new protein", {
    expect_error(
        setProtein(
            gene1, 
            "ENSP01234567891", 
            list(
                id = "ENSP00000000001", 
                protein_sequence = AAString("ARNDCQ"), 
                description = "It has no specific function")), 
                "Invalid protein"
        )
})

# Tests for the function lengthProductPCG

# Testing the function lengthProductPCG, invalid input
test_that("invalid input for lengthProductPCG", {
    gene2 <- "This is not a gene"
    expect_error(lengthProductPCG(gene2),
    "Invalid function object. It should be of class 'ProteinCodingGene'")

})

# Testing the function lengthProductPCG, expected output
test_that("function lengthProductPCG executed correctly", {

    expected_output <- vapply(gene1@proteins, function(x) {length(x@protein_sequence)}, integer(1))
    names(expected_output) <- c("ENSP01234567891", "ENSP19876543210")
    expect_equal(lengthProductPCG(gene1), expected_output)

})