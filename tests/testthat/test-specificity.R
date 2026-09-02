sam_line <- function(name, flag, chrom, pos, nm, md, length = 26L) {
  bases <- strrep("A", length)
  paste(name, flag, chrom, pos, 42L, paste0(length, "M"), "*", 0L, 0L, bases, strrep("I", length),
        "AS:i:0", paste0("NM:i:", nm), paste0("MD:Z:", md), sep = "\t")
}

test_that("MD tags give mismatch offsets along the reference", {
  expect_equal(md_mismatch_offsets("26"), integer())
  expect_equal(md_mismatch_offsets("0T25"), 1L)
  expect_equal(md_mismatch_offsets("10A5^AC6"), 11L)
  expect_equal(md_mismatch_offsets("24G0C"), c(25L, 26L))
})

test_that("SAM hits record orientation, Bismark chromosome names and 3' mismatches", {
  sam <- c("@HD\tVN:1.0",
           sam_line("pair1_left", 0L, "chr1_CT_converted", 101L, 0L, "26"),
           sam_line("pair1_left", 0L, "chr2_CT_converted", 101L, 1L, "24A1"),
           sam_line("pair1_right", 16L, "chr1_CT_converted", 400L, 1L, "10A15"),
           sam_line("pair1_right", 16L, "chr3_CT_converted", 400L, 1L, "1A24"),
           paste("pair1_right", 4L, "*", 0L, 0L, "*", "*", 0L, 0L, "A", "I", sep = "\t"))
  hits <- parse_bowtie2_sam(sam, "CT")
  expect_equal(nrow(hits), 4L)
  expect_equal(hits$chrom, c("chr1", "chr2", "chr1", "chr3"))
  expect_equal(hits$reverse, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(hits$end, c(126L, 126L, 425L, 425L))
  # Forward 3' end is the alignment end; reverse 3' end is the alignment start.
  expect_equal(hits$three_prime_mismatch, c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(nrow(parse_bowtie2_sam(character(), "GA")), 0L)
})

test_that("products need a forward hit followed by a nearby reverse hit", {
  hits <- data.frame(primer = c("L", "R", "L", "R", "R"), index = c("CT", "CT", "CT", "CT", "GA"),
                     chrom = c("c1", "c1", "c2", "c2", "c1"), start = c(101L, 400L, 100L, 9000L, 420L),
                     end = c(126L, 425L, 125L, 9025L, 445L), reverse = c(FALSE, TRUE, FALSE, TRUE, TRUE),
                     mismatches = 0L, three_prime_mismatch = FALSE)
  products <- predict_products(hits, max_product = 3000L)
  expect_equal(nrow(products), 1L)
  expect_equal(c(products$start, products$end, products$length), c(101L, 425L, 325L))
  expect_equal(c(products$forward_primer, products$reverse_primer), c("L", "R"))
  expect_equal(nrow(predict_products(hits, max_product = 300L)), 0L)
})

test_that("Bowtie2 finds the on-target product and a planted off-target copy", {
  skip_if_no_primer3()
  skip_if_no_bismark()
  design <- design_bsp_primers(demo_promoter("BSPPLUS"), "sense")
  prepare <- function(extra = NULL) {
    folder <- tempfile("bismark")
    dir.create(folder)
    genome <- Biostrings::readDNAStringSet(system.file("extdata", "bsp_demo.fa", package = "promoterPrimerMapper"))
    if (!is.null(extra)) genome <- c(genome, Biostrings::DNAStringSet(extra))
    Biostrings::writeXStringSet(genome, file.path(folder, "genome.fa"))
    status <- system2(Sys.which("bismark_genome_preparation"), c("--bowtie2", folder), stdout = FALSE, stderr = FALSE)
    expect_equal(status, 0L)
    folder
  }
  unique_genome <- prepare()
  checked <- check_bsp_specificity(design, unique_genome, top_n = nrow(design$pairs))
  expect_true(all(checked$pairs$on_target))
  expect_true(all(checked$pairs$specific))
  expect_true(all(checked$specificity$products$index == "CT"))

  set.seed(11)
  flank <- function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
  copy <- c(chrCopy = paste0(flank(500L), bsp_amplicons(design, 1L)[[1L]], flank(500L)))
  copied <- check_bsp_specificity(design, prepare(copy), top_n = 1L)
  expect_true(copied$pairs$on_target[1L])
  expect_equal(copied$pairs$n_offtarget_products[1L], 1L)
  expect_false(copied$pairs$specific[1L])
  # top_n counts per target: the best pair of every target is checked, nothing else.
  checked_rows <- which(!is.na(copied$pairs$specific))
  expect_equal(checked_rows, match(unique(design$pairs$target_id), design$pairs$target_id))
  expect_true(all(copied$pairs$specific[checked_rows[-1L]]))
  off_target <- copied$specificity$products[copied$specificity$products$chrom == "chrCopy", ]
  expect_equal(c(off_target$start, off_target$end), c(501L, 500L + design$pairs$product_size[1L]))
})

test_that("a missing Bismark index is reported", {
  expect_error(bismark_indexes(tempdir()), "bismark_genome_preparation")
})
