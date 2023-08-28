# Testing if a microRNAGene gene object is define in a correct way
test_that("correct microRNAGene object", {

    gene1 <- microRNAGene(
    microRNAs = list(
        microRNA(
            id = "ENST01234567891", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has .."),
        microRNA(
            id = "ENST19876543210", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has ..")),
    geneID = "ENSG01234567891",
    gene_symbol = "MIR 12",
    full_gene_name = "microRNA 23",
    description = "This gene has a specific function",
    structure = GRanges("chr1:1-1000")
    )

    expect_no_error(gene1)

})

# Tests regarding the validity of object of class microRNAGene

# Testing an object for class microRNAGene with invalid EnsemblID of the gene 
test_that("invalid EnsemblID of the gene for microRNAGene object", {

    expect_error(
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ),"The gene enemblID is not correct")
    }
)

# Testing an object for class microRNAGene with invalid HUGO gene symbol
test_that("invalid gene symbol for microRNAGene object", {

    expect_error(
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "mir 12",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ),"Invalid gene symbol")
    }
)

# Testing an object for class microRNAGene with invalid full gene name
test_that("invalid gene symbol for microRNAGene object", {

    expect_error(
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microrna 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ),'Invalid full gene name symbol')
    }
)

# Testing an object for class microRNAGene with invalid proteins
test_that("invalid protein objetc for microRNAGene object", {

    expect_error(
        microRNAGene(
        microRNAs = list(
            id = "ENST01234567891", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has .."),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ), "The elements of the microRNAs list must be of class 'microRNA'")
    }
)

# Testing an object for class microRNAGene with invalid EnsemblID of the protein
test_that("invalid EnsemblID of the microRNA for microRNAGene object", {

    expect_error(
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENSP01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ), "The microRNA enemblID is not correct")
    }
)

# Testing an object for class microRNAGene with a replicate proteins
test_that("replicates in microRNAs for microRNAGene object", {

    expect_error(
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        ), "No replicates in microRNAs")
    }
)

# Tests of all the method of the class microRNAGene
# Before doing the test, an object of class microRNAGene was defined
gene1 <- microRNAGene(
    microRNAs = list(
        microRNA(
            id = "ENST01234567891", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has .."),
        microRNA(
            id = "ENST19876543210", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has ..")),
    geneID = "ENSG01234567891",
    gene_symbol = "MIR 12",
    full_gene_name = "microRNA 23",
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
    expect_equal(getSymbol(gene1), "MIR 12")
})

# Testing getFullName method
test_that("method getSymbol executed correctly", {
    expect_equal(getFullName(gene1), "microRNA 23")
})

# Testing getDescription method
test_that("method getDescription executed correctly", {
    expect_equal(getDescription(gene1), "This gene has a specific function")
})

# Testing getStructure method
test_that("method getStructure executed correctly", {
    expect_equal(getStructure(gene1), GRanges("chr1:1-1000"))
})

# Testing getMicroRNA method
test_that("method getMicroRNA executed correctly", {
    expect_equal(getMicroRNA(gene1), list(
        microRNA(
            id = "ENST01234567891", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has .."),
        microRNA(
            id = "ENST19876543210", 
            microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
            silencing_mechanism = "The microRNA has ..")))
})

# Tests of all the setter methods

# Testing setSymbol method
test_that("method setSymbol executed correctly", {

    expect_equal(
        setSymbol(gene1, "MIR 2"),
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210",
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 2",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})


# Testing setSymbol method with invalid new gene symbol
test_that("method setSymbol with invalid gene symbol", {
    expect_error(
        setSymbol(gene1, "mir 2"), "Invalid gene symbol"
        )
})

# Testing setFullName method

test_that("method setFullName executed correctly", {
    expect_equal(
        setFullName(gene1, "microRNA 2"),
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 2",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setFullName method with invalid new full name gene
test_that("method setFullName with invalid new full name gene", {
    expect_error(
        setFullName(gene1, "microrna 2"), "Invalid full gene name"
        )
})

# Testing setDescription method
test_that("method setDescription executed correctly", {
    expect_equal(
        setDescription(gene1, "This gene has no specific function"),
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
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
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST01234567891", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has .."),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
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

# Testing setMicroRNA method
test_that("method setMicroRNA executed correctly", {
    expect_equal(
        setMicroRNA(
            gene1, 
            "ENST01234567891", 
            microRNA(
                id = "ENST00000000001", 
                microRNA_sequence = RNAString("ACG"), 
                silencing_mechanism = "Its silencing mechanism consists of")),
        microRNAGene(
        microRNAs = list(
            microRNA(
                id = "ENST00000000001", 
                microRNA_sequence = RNAString("ACG"), 
                silencing_mechanism = "Its silencing mechanism consists of"),
            microRNA(
                id = "ENST19876543210", 
                microRNA_sequence = RNAString("ACGUMRWSYKVHDBN-+."), 
                silencing_mechanism = "The microRNA has ..")),
        geneID = "ENSG01234567891",
        gene_symbol = "MIR 12",
        full_gene_name = "microRNA 23",
        description = "This gene has a specific function",
        structure = GRanges("chr1:1-1000")
        )
    )
})

# Testing setMicroRNA method with invalid new microRNA
test_that("method setMicroRNA with invalid new microRNA", {
    expect_error(
        setMicroRNA(
            gene1, 
            "ENST01234567891", 
            list(
                id = "ENSP00000000001",
                microRNA_sequence = RNAString("ACG"), 
                silencing_mechanism = "Its silencing mechanism consists of")), 
                "Invalid microRNA"
        )
})

# Tests for the function lengthProductMRG

# Testing the function lengthProductMRG, invalid input
test_that("invalid input for lengthProductMRG", {
    gene2 <- "This is not a gene"
    expect_error(lengthProductMRG(gene2), 
    "Invalid function object. It should be of class 'MicroRNAGene'")

})

# Testing the function lengthProductMRG, expected output
test_that("function lengthProductMRG executed correctly", {

    expected_output <- vapply(gene1@microRNAs, function(x) {length(x@microRNA_sequence)}, integer(1))
    names(expected_output) <- c("ENST01234567891", "ENST19876543210")
    expect_equal(lengthProductMRG(gene1), expected_output)

})