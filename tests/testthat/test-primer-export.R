export_design <- function() {
  pairs <- data.frame(rank = c(3L, 7L), target_id = c("CGI1", "CGI2"), strand = c("sense", "antisense"),
    assay_strand = c("OT", "OB"), scaffold = c('scaffold,"one"', "scaffold2"),
    left_sequence = c("ACGTT", "TTGCA"), right_sequence = c("TAGCA", "GGCAT"),
    left_length = 5L, right_length = 5L,
    left_genomic_start = c(101L, 901L), left_genomic_end = c(105L, 905L),
    right_genomic_start = c(196L, 806L), right_genomic_end = c(200L, 810L),
    left_tm = c(56.1, 57.2), right_tm = c(55.3, 56.4), left_gc = c(40, 60), right_gc = c(60, 80),
    left_self_any = c(3.1, 4.2), right_self_any = c(5.3, 6.4),
    left_self_end = c(1.1, 2.2), right_self_end = c(3.3, 4.4),
    template_left_start = c(1L, 101L), template_right_end = c(100L, 200L),
    product_size = 100L, product_tm = c(72.5, 74.5), cpg_in_amplicon = c(4L, 6L))
  structure(list(pairs = pairs, params = bsp_primer_params(), extracted = list(gene_name = "GENE", scaffold = "scaffold1",
    region_start = 1L, region_end = 1000L, strand = "+", tss = 500L),
    # Two CpGs overlap primer footprints. A CpG starting at the final product
    # base extends beyond it and must not be counted as a complete site.
    templates = list(sense = list(cpg_c = seq_len(201L) %in% c(2L, 20L, 30L, 40L, 50L, 98L, 100L)),
      antisense = list(cpg_c = seq_len(201L) %in% c(101L, 108L, 120L, 130L, 140L, 150L, 170L, 199L, 200L))),
    targets = data.frame(target_id = c("CGI1", "CGI2"), tss_start = c(-100L, 100L), tss_end = c(-1L, 200L)),
    gui_provenance = list(completed = "today\ntest run", extraction = list(annotation = list(path = "annotation\r\nsecond line.gtf"),
      fasta = list(path = "reference.fa"), index = list(path = "reference.fa.fai")))), class = "bsp_design")
}

test_that("primer exports have commented metadata then the exact CSV schema", {
  design <- export_design()
  path <- tempfile(fileext = ".csv")
  expect_identical(write_order_sheet(design, path), path)
  lines <- readLines(path)
  header <- which(!startsWith(lines, "#"))[1L]
  expect_gt(header, 3L)
  expect_true(all(startsWith(lines[seq_len(header - 1L)], "# ")))
  expect_match(lines[header], '^"Target","Primer","Sequence', fixed = FALSE)
  expect_false(any(startsWith(lines[header:length(lines)], "#")))
  expect_true("# second line.gtf" %in% lines)
  expect_true("# test run" %in% lines)
  table <- utils::read.csv(path, comment.char = "#", check.names = FALSE)
  expect_named(table, c("Target", "Primer", "Sequence (5'->3')", "Primer Length", "Scaffold", "Primer Start", "Primer Stop",
    "Tm", "GC%", "Self complementarity", "Self 3' complementarity", "Product length", "Product Tm", "CpGs in Product"))
  expect_equal(table$Target, c("CGI1", "CGI1", "CGI2", "CGI2"))
  expect_equal(table$Primer, c("GENE_CGI1_BSP3_F", "GENE_CGI1_BSP3_R", "GENE_CGI2_BSP7_F", "GENE_CGI2_BSP7_R"))
  expect_equal(table[["Sequence (5'->3')"]], c("ACGTT", "TAGCA", "TTGCA", "GGCAT"))
  expect_equal(table$Scaffold, rep(design$pairs$scaffold, each = 2L))
  expect_equal(table[["Primer Start"]], c(101L, 196L, 901L, 806L))
  expect_equal(table[["Primer Stop"]], c(105L, 200L, 905L, 810L))
  expect_equal(table[["Primer Stop"]] - table[["Primer Start"]] + 1L, table[["Primer Length"]])
  expect_equal(table$Tm, c(56.1, 55.3, 57.2, 56.4))
  expect_equal(table[["GC%"]], c(40, 60, 60, 80))
  expect_equal(table[["Self complementarity"]], c(3.1, 5.3, 4.2, 6.4))
  expect_equal(table[["Self 3' complementarity"]], c(1.1, 3.3, 2.2, 4.4))
  expect_equal(table[["Product length"]], rep(100L, 4L))
  expect_equal(table[["Product Tm"]], c(72.5, 72.5, 74.5, 74.5))
  expect_equal(table[["CpGs in Product"]], c(6L, 6L, 8L, 8L))
})

test_that("primer exports respect requested pair order and handle legacy product Tm", {
  design <- export_design()
  path <- tempfile(fileext = ".csv")
  write_order_sheet(design, path, ranks = c(7L, 3L))
  table <- utils::read.csv(path, comment.char = "#", check.names = FALSE)
  expect_equal(table$Primer, c("GENE_CGI2_BSP7_F", "GENE_CGI2_BSP7_R", "GENE_CGI1_BSP3_F", "GENE_CGI1_BSP3_R"))
  design$pairs$product_tm <- NULL
  write_order_sheet(design, path, ranks = 7L)
  table <- utils::read.csv(path, comment.char = "#", check.names = FALSE)
  expect_equal(nrow(table), 2L)
  expect_true(all(is.na(table[["Product Tm"]])))
  expect_error(write_order_sheet(design, path, ranks = 999L), "not in the design")
})
