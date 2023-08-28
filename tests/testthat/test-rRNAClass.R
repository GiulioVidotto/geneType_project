# Testing if a rRNA object is define in a correct way
test_that("correct rRNA object", {

   gene1 <- rRNAgene(
    rRNA = rRNA(
        id = "ENST01234567891",
        rRNA_sequence = RNAString("AUUUGA"),
        ribosomal_subunit = "28S"),
    geneID = "ENSG01234567891",
    gene_symbol = "RNA28S",
    full_gene_name = "RNA, 28S ribosomal",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Testing if a rRNA object with invalid gene enemblID
test_that("invalid gene enemblID for rRNAgene object", {

    expect_error(
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG012345",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The gene enemblID is not correct")

})

# Testing if a rRNA object with invalid gene symbol
test_that("invalid gene symbol for rRNAgene object", {

    expect_error(
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA45S",
        full_gene_name = "RNA, 28S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid gene symbol")

})


# Testing if a rRNA object with invalid full name
test_that("invalid full gene name for rRNAgene object", {

    expect_error(
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 45S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid full gene name symbol")

})

# Testing an object with invalid enemblID rRNA
test_that("invalid rRNA enemblID for rRNAgene object", {

    expect_error(
        rRNAgene(
        rRNA = rRNA(
            id = "ENSP01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The rRNA enemblID is not correct")

})

# Testing an object with invalid aminoacid
test_that("invalid aminoacid for rRNAgene object", {

    expect_error(
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "2S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The ribosomal subunit is not correct")

})

# Define the object of class rRNA
gene1 <- rRNAgene(
    rRNA = rRNA(
        id = "ENST01234567891",
        rRNA_sequence = RNAString("AUUUGA"),
        ribosomal_subunit = "28S"),
    geneID = "ENSG01234567891",
    gene_symbol = "RNA28S",
    full_gene_name = "RNA, 28S ribosomal",
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
    expect_equal(getSymbol(gene1), "RNA28S")
})

# Testing getFullName method
test_that("method getSymbol executed correctly", {
    expect_equal(getFullName(gene1), "RNA, 28S ribosomal")
})

# Testing getDescription method
test_that("method getDescription executed correctly", {
    expect_equal(getDescription(gene1), "This gene has a specific function")
})

# Testing getStructure method
test_that("method getStructure executed correctly", {
    expect_equal(getStructure(gene1), GRanges("chr1:1-1000"))
})

# Testing getRRNA method
test_that("method getRRNA executed correctly", {
    expect_equal(getRRNA(gene1), rRNA(
        id = "ENST01234567891",
        rRNA_sequence = RNAString("AUUUGA"),
        ribosomal_subunit = "28S"))
})

# Tests of all the setter methods

# Testing setSymbol method
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "RNA28S"),
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with invalid new gene symbol
test_that("method setSymbol with invalid gene symbol", {
    expect_error(
        setSymbol(gene1, "RNA59S"), "Invalid gene symbol"
        )
})

# Testing setFullName method

test_that("method setFullName executed correctly", {
    expect_equal(
        setFullName(gene1, "RNA, 28S ribosomal"),
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setFullName method with invalid new full name gene
test_that("method setFullName with invalid new full name gene", {
    expect_error(
        setFullName(gene1, "RNA, 88S ribosomal"), "Invalid full gene name"
        )
})


# Testing setRibosomalSub method
test_that("method setRibosomalSub executed correctly", {

    expect_equal(
        setRibosomalSub(gene1, "45S"),
        rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "45S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA45S",
        full_gene_name = "RNA, 45S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setAminoacid method with invalid new aminoacid
test_that("method setAminoacid with invalid new aminoacid", {
    expect_error(
        setRibosomalSub(gene1, "M"), "Invalid Ribosomal Subunit"
        )
})

# Testing setDescription method
test_that("method setDescription executed correctly", {
    expect_equal(setDescription(gene1, "This gene has no specific function"),
    rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
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
    rRNAgene(
        rRNA = rRNA(
            id = "ENST01234567891",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "28S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA28S",
        full_gene_name = "RNA, 28S ribosomal",
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

# Testing setRRNA method for object of class ProteinCodingGene
test_that("method setRRNA executed correctly", {
    expect_equal(
        setRRNA(gene1,
            rRNA(id = "ENST00000000001",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "45S")),
        rRNAgene(
        rRNA = rRNA(
            id = "ENST00000000001",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "45S"),
        geneID = "ENSG01234567891",
        gene_symbol = "RNA45S",
        full_gene_name = "RNA, 45S ribosomal",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setRRNA method with invalid new rRNA
test_that("method setRRNA with invalid new rRNA", {
    expect_error(
        setRRNA(gene1, list(id = "ENST00000000001",
            rRNA_sequence = RNAString("AUUUGA"),
            ribosomal_subunit = "45S")),
            "Invalid rRNA"
        )
})

# Tests for the function lengthProductRRG

# Testing the function lengthProductRRG, invalid input
test_that("invalid input for lengthProductRRG", {
    gene2 <- "This is not a gene"
    expect_error(lengthProductRRG(gene2),
    "Invalid function object. It should be of class 'rRNAgene'")

})

# Testing the function lengthProductRRG, expected output
test_that("function lengthProductRRG executed correctly", {

    expected_output <- length(gene1@rRNA@rRNA_sequence)
    names(expected_output) <- c("ENST01234567891")
    expect_equal(lengthProductRRG(gene1), expected_output)

})