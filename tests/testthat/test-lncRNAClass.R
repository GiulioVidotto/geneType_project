# Testing if a LongNonCodingRNAgene object is define in a correct way
test_that("correct LongNonCodingRNAgene object", {

    gene1 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(
            id = "ENST01234567891", 
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210", 
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "LINC01018",
    full_gene_name = "long intergenic non-protein coding RNA 01018",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Testing if a LongNonCodingRNAgene object is define in a correct way
test_that("correct LongNonCodingRNAgene object", {

    gene1 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(
            id = "ENST01234567891",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210", 
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "ABCF1-DT",
    full_gene_name = "ABCF1 divergent transcript",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Testing if a LongNonCodingRNAgene object is define in a correct way
test_that("correct LongNonCodingRNAgene object", {

    gene1 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(
            id = "ENST01234567891", 
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210", 
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "AOAH-IT1",
    full_gene_name = "AOAH intronic transcript 1",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Testing if a LongNonCodingRNAgene object is define in a correct way
test_that("correct LongNonCodingRNAgene object", {

    gene1 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(
            id = "ENST01234567891",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210", 
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "C5-OT1",
    full_gene_name = "C5 3' UTR overlapping transcript 1",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Testing if a LongNonCodingRNAgene object is define in a correct way
test_that("correct LongNonCodingRNAgene object", {

    gene1 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(
            id = "ENST01234567891",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "FAS-AS1",
    full_gene_name = "FAS antisense RNA 1",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

    expect_no_error(gene1)

})

# Tests regarding the validity of object of class LongNonCodingRNAgene

# Testing an object for class LongNonCodingRNAgene with invalid EnsemblID of the gene 
test_that("invalid EnsemblID of the gene for LongNonCodingRNAgene object", {

    expect_error(
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG0123",
        gene_symbol = "LINC1018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ),"The gene enemblID is not correct")
    }
)

# Testing an object for class longNonCodingRNAGene with invalid HUGO gene symbol
test_that("invalid gene symbol for longNonCodingRNAGene object", {

    expect_error(
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "gene 1",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ),"Invalid gene symbol")
    }
)

# Testing an object for class longNonCodingRNAGene with invalid full gene name
test_that("invalid full gene name for longNonCodingRNAGene object", {

    expect_error(
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ),"The full gene name of LncRNAgenes that are intergenic with respect to protein coding genes should start with 'long intergenic non-protein coding RNA'")
    }
)

# Testing an object for class longNonCodingRNAGene with invalid proteins
test_that("invalid lncRNA objetc for longNonCodingRNAGene object", {

    expect_error(
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            id = "ENST01234567891",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ), "The elements of the LongNonCodingRNAs list must be of class 'LongNonCodingRNA'")
    }
)

# Testing an object for class longNonCodingRNAGene with invalid EnsemblID of the protein
test_that("invalid EnsemblID of the longNonCodingRNA for longNonCodingRNAGene object", {

    expect_error(
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENSP01234567891",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ), "The long non codingRNA enemblID is not correct")
    }
)

# Testing an object for class longNonCodingRNAGene with a replicate proteins
test_that("replicates in longNonCodingRNAs for longNonCodingRNAGene object", {

    expect_error(
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST01234567891",
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ), "No replicates in LongNonCodingRNAs")
    }
)

# Tests of all the method of the class longNonCodingRNAGene
# Before doing the test, an object of class longNonCodingRNAGene was defined
gene1 <- longNonCodingRNAGene(
    LongNonCodingRNAs = list(
        longNonCodingRNA(
            id = "ENST01234567891",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism")),
    geneID = "ENSG01234567891",
    gene_symbol = "LINC01018",
    full_gene_name = "long intergenic non-protein coding RNA 01018",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000"))

# Tests of all the getter methods

# Testing getGeneId method
test_that("method getGeneId executed correctly", {
    expect_equal(getGeneId(gene1), "ENSG01234567891")
})

# Testing getSymbol method
test_that("method getSymbol executed correctly", {
    expect_equal(getSymbol(gene1), "LINC01018")
})

# Testing getFullName method
test_that("method getSymbol executed correctly", {
    expect_equal(getFullName(gene1),
    "long intergenic non-protein coding RNA 01018")
})

# Testing getDescription method
test_that("method getDescription executed correctly", {
    expect_equal(getDescription(gene1), "This gene has a specific function")
})

# Testing getStructure method
test_that("method getStructure executed correctly", {
    expect_equal(getStructure(gene1), GRanges("chr1:1-1000"))
})

# Testing getLongNonCodingRNA method
test_that("method getLongNonCodingRNA executed correctly", {
    expect_equal(getLongNonCodingRNA(gene1), list(
        longNonCodingRNA(
            id = "ENST01234567891",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism"),
        longNonCodingRNA(
            id = "ENST19876543210",
            long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."),
            regulatory_mechanism = "Description of the specific regulatory mechanism")))
})

# Tests of all the setter methods

# Testing setSymbol method for long intergenic non-protein coding RNA
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "LINC10000"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC10000",
        full_gene_name = "long intergenic non-protein coding RNA 10000",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setSymbol method with for LncRNAs that are antisense to the genomic span of a protein coding gene
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "FAS-AS1"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "FAS-AS1",
        full_gene_name = "FAS antisense RNA 1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with for LncRNAs that are divergent
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "ABCF1-DT"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "ABCF1-DT",
        full_gene_name = "ABCF1 divergent transcript",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with for LncRNAs that are contained within an intron of a protein coding gene
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "AOAH-IT1"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "AOAH-IT1",
        full_gene_name = "AOAH intronic transcript 1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with for LncRNAs that overlap a protein coding gene
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "C5-OT1"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "C5-OT1",
        full_gene_name = "C5 3' UTR overlapping transcript 1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with invalid new gene symbol
test_that("method setSymbol with invalid gene symbol", {
    expect_error(
        setSymbol(gene1, "gene 2"), "Invalid gene symbol"
        )
})


# Testing setFullName method for long intergenic non-protein coding RNA
test_that("method setFullName executed correctly", {

    expect_equal(
        setFullName(gene1, "long intergenic non-protein coding RNA 10000"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC10000",
        full_gene_name = "long intergenic non-protein coding RNA 10000",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setFullName method with for LncRNAs that are antisense to the genomic span of a protein coding gene
test_that("method setFullName executed correctly", {

    expect_equal(
        setFullName(gene1, "FAS antisense RNA 1"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "FAS-AS1",
        full_gene_name = "FAS antisense RNA 1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setFullName method with for LncRNAs that are divergent
test_that("method setFullName executed correctly", {

    expect_equal(
        setFullName(gene1, "ABCF1 divergent transcript"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "ABCF1-DT",
        full_gene_name = "ABCF1 divergent transcript",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setFullName method with for LncRNAs that are contained within an intron of a protein coding gene
test_that("method setFullName executed correctly", {

    expect_equal(
        setFullName(gene1, "AOAH intronic transcript 1"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "AOAH-IT1",
        full_gene_name = "AOAH intronic transcript 1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setFullName method with for LncRNAs that overlap a protein coding gene
test_that("method setFullName executed correctly", {

    expect_equal(
        setFullName(gene1, "C5 3' UTR overlapping transcript 1"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "C5-OT1",
        full_gene_name = "C5 3' UTR overlapping transcript 1",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setFullName method with invalid full gene name
test_that("method setFullName with invalid full gene name", {
    expect_error(
        setFullName(gene1, "new full name"), "Invalid full gene name"
        )
})


# Testing setDescription method
test_that("method setDescription executed correctly", {
    expect_equal(
        setDescription(gene1, "This gene has no specific function"),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
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
    expect_equal(
        setStructure(gene1, GRanges("chr1:100-200")),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST01234567891", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:100-200")
        )
    )
})

# Testing setStructure method with invalid new structure
test_that("method setStructure with invalid structure", {
    expect_error(setStructure(gene1, "chr1:100-200"), "The new defined structure of the gene is not correct"
    )
})

# Testing setlongNonCodingRNA method
test_that("method setlongNonCodingRNA executed correctly", {
    expect_equal(
        setlongNonCodingRNA(
            gene1, 
            "ENST01234567891", 
            longNonCodingRNA(
                id = "ENST00000000001", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        longNonCodingRNAGene(
        LongNonCodingRNAs = list(
            longNonCodingRNA(
                id = "ENST00000000001", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism"),
            longNonCodingRNA(
                id = "ENST19876543210", 
                long_non_coding_RNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                regulatory_mechanism = "Description of the specific regulatory mechanism")),
        geneID = "ENSG01234567891",
        gene_symbol = "LINC01018",
        full_gene_name = "long intergenic non-protein coding RNA 01018",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setlongNonCodingRNA method with invalid new Long Non Coding RNA
test_that("method setlongNonCodingRNA with invalid new Long Non Coding RNA", {
    expect_error(
        setlongNonCodingRNA(
            gene1, 
            "ENST01234567891", 
            list(
                id = "ENSP00000000001", 
                microRNA_sequence = RNAString("ACG"), 
                silencing_mechanism = "Its silencing mechanism consists of")), 
                "Invalid longNonCodingRNA"
        )
})

# Tests for the function lengthProductLNCRG

# Testing the function lengthProductLNCRG, invalid input
test_that("invalid input for lengthProductLNCRG", {
    gene2 <- "This is not a gene"
    expect_error(lengthProductLNCRG(gene2), 
    "Invalid function object. It should be of class 'longNonCodingRNAGene'")

})

# Testing the function lengthProductLNCRG, expected output
test_that("function lengthProductLNCRG executed correctly", {

    expected_output <- vapply(gene1@LongNonCodingRNAs, function(x) {length(x@long_non_coding_RNA_sequence)}, integer(1))
    names(expected_output) <- c("ENST01234567891", "ENST19876543210")
    expect_equal(lengthProductLNCRG(gene1), expected_output)

})