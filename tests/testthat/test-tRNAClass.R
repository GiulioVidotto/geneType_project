# Testing if a tRNA object is define in a correct way
test_that("correct tRNA object", {

    gene1 <- tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"),
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Testing if a tRNA object with invalid gene enemblID
test_that("invalid gene enemblID for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG0123456",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The gene enemblID is not correct")

})

# Testing if a tRNA object with invalid gene symbol
test_that("invalid gene symbol for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AUU1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid gene symbol")

})


# Testing if a tRNA object with invalid full name
test_that("invalid full gene name for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAUU 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid full gene name symbol")

})

# Testing an object with invalid enemblID tRNA
test_that("invalid tRNA enemblID for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST019876", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "The tRNA enemblID is not correct")

})

# Testing an object with invalid aminoacid
test_that("invalid aminoacid for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("AAA"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid aminoacid")

})

# Testing an object with invalid aminoacid
test_that("invalid aminoacid for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("2"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid aminoacid")

})

# Testing an object with invalid anticodon
test_that("invalid anticodon for tRNAgene object", {

    expect_error(
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
            aminoacid = AAString("A"),
            anticodon = RNAString("AGCCGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
    ), "Invalid anticodon")

})

# Define the object of class tRNA
gene1 <- tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321",
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."),
            aminoacid = AAString("A"),
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000"))


# Tests of all the getter methods

# Testing getGeneId method
test_that("method getGeneId executed correctly", {
    expect_equal(getGeneId(gene1), "ENSG01234567891")
})

# Testing getSymbol method
test_that("method getSymbol executed correctly", {
    expect_equal(getSymbol(gene1), "TRA-AGC1-1")
})

# Testing getFullName method
test_that("method getSymbol executed correctly", {
    expect_equal(getFullName(gene1), "tRNA-ALAAGC 1-1")
})

# Testing getDescription method
test_that("method getDescription executed correctly", {
    expect_equal(getDescription(gene1), "This gene has a specific function")
})

# Testing getStructure method
test_that("method getStructure executed correctly", {
    expect_equal(getStructure(gene1), GRanges("chr1:1-1000"))
})

# Testing getTRNA method
test_that("method getTRNA executed correctly", {
    expect_equal(getTRNA(gene1), tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")))
})

# Tests of all the setter methods

# Testing setSymbol method
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "TRM-AUU1-1"),
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("M"),
            anticodon = RNAString("AUU")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRM-AUU1-1",
        full_gene_name = "tRNA-METAUU 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with invalid new gene symbol
test_that("method setSymbol with invalid gene symbol", {
    expect_error(
        setSymbol(gene1, "TR3-BBB1-1"), "Invalid gene symbol"
        )
})

# Testing setFullName method

test_that("method setFullName executed correctly", {
    expect_equal(
        setFullName(gene1, "tRNA-METAUU 1-1"),
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("M"),
            anticodon = RNAString("AUU")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRM-AUU1-1",
        full_gene_name = "tRNA-METAUU 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setFullName method with invalid new full name gene
test_that("method setFullName with invalid new full name gene", {
    expect_error(
        setFullName(gene1, "tRNA-123AB 1-1"), "Invalid full gene name"
        )
})

# Testing setAminoacid method
test_that("method setAminoacid executed correctly", {

    expect_equal(
        setAminoacid(gene1, AAString("M")),
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("M"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRM-AGC1-1",
        full_gene_name = "tRNA-METAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setAminoacid method with invalid new aminoacid
test_that("method setAminoacid with invalid new aminoacid", {
    expect_error(
        setAminoacid(gene1, "M"), "Invalid aminoacid"
        )
})

# Testing setAnticodon method
test_that("method setAnticodon executed correctly", {

    expect_equal(
        setAnticodon(gene1, RNAString("AUU")),
        tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", tRNA_sequence = RNAString("AUUCGUMRWSYKVHDBN-+."), aminoacid = AAString("A"), anticodon = RNAString("AUU")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AUU1-1",
        full_gene_name = "tRNA-ALAAUU 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setAnticodon method with invalid new anticodon
test_that("method setAnticodon with invalid new anticodon", {
    expect_error(
        setAnticodon(gene1, "AUU"), "Invalid anticodon"
        )
})

# Testing setDescription method
test_that("method setDescription executed correctly", {
    expect_equal(setDescription(gene1, "This gene has no specific function"),
    tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has no specific function",
        structure = GRanges("chr1:1-1000"))
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
    tRNAgene(
        tRNA = tRNA(
            id = "ENST01987654321", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("A"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRA-AGC1-1",
        full_gene_name = "tRNA-ALAAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:100-200"))
    )
})

# Testing setStructure method with invalid new structure
test_that("method setStructure with invalid structure", {
    expect_error(setStructure(gene1, "chr1:100-200"), 
    "The new defined structure of the gene is not correct"
    )
})

# Testing setTRNA method for object of class ProteinCodingGene
test_that("method setTRNA executed correctly", {
    expect_equal(
        setTRNA(
            gene1, 
            tRNA(
                id = "ENST00000000001", 
                tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
                aminoacid = AAString("M"), 
                anticodon = RNAString("AGC"))),
        tRNAgene(
        tRNA = tRNA(
            id = "ENST00000000001", 
            tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
            aminoacid = AAString("M"), 
            anticodon = RNAString("AGC")),
        geneID = "ENSG01234567891",
        gene_symbol = "TRM-AGC1-1",
        full_gene_name = "tRNA-METAGC 1-1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000"))
    )
})

# Testing setProtein method with invalid new tRNA
test_that("method setProtein with invalid new tRNA", {
    expect_error(
        setTRNA(
            gene1,
            list(
                id = "ENST00000000001", 
                tRNA_sequence = RNAString("AGCCGUMRWSYKVHDBN-+."), 
                aminoacid = AAString("M"), 
                anticodon = RNAString("AGC"))), 
                "Invalid tRNA"
        )
})

# Tests for the function lengthProductTRG

# Testing the function lengthProductTRG, invalid input
test_that("invalid input for lengthProductTRG", {
    gene2 <- "This is not a gene"
    expect_error(lengthProductTRG(gene2), 
    "Invalid function object. It should be of class 'tRNAgene'")

})

# Testing the function lengthProductTRG, expected output
test_that("function lengthProductTRG executed correctly", {

    expected_output <- length(gene1@tRNA@tRNA_sequence)
    names(expected_output) <- c("ENST01987654321")
    expect_equal(lengthProductTRG(gene1), expected_output)

})