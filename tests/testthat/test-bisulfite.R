test_that("bisulfite conversion can keep CpG cytosines", {
  expect_equal(bisulfite_convert("CCGACGTC"), "TTGATGTT")
  expect_equal(bisulfite_convert("CCGACGTC", cpg = "methylated"), "TCGACGTT")
  expect_equal(bisulfite_convert("cgcg", cpg = "methylated"), "CGCG")
})

test_that("TSS-relative positions skip zero", {
  expect_equal(tss_relative_position(1:5, 3L), c(-2L, -1L, 1L, 2L, 3L))
  expect_equal(sense_position_from_tss(c(-2L, -1L, 1L, 2L, 3L), 3L), 1:5)
  expect_error(sense_position_from_tss(0L, 3L), "0 is not a valid")
})

test_that("CpG sites carry sense, TSS-relative and genomic positions on both strands", {
  plus <- cpg_sites(promoter("ACGTTCGA", "+", tss = 102L))
  expect_equal(plus$position, c(2L, 6L))
  expect_equal(plus$tss_relative, c(-1L, 4L))
  expect_equal(plus$genomic_start, c(101L, 105L))
  minus <- cpg_sites(promoter("ACGTTCGA", "-", tss = 106L))
  expect_equal(minus$position, c(2L, 6L))
  expect_equal(minus$tss_relative, c(1L, 5L))
  expect_equal(minus$genomic_start, c(105L, 101L))
  expect_equal(minus$genomic_end, c(106L, 102L))
  expect_equal(cpg_summary(promoter("ACGTTCGA"))$cpg_count, 2L)
})

test_that("bisulfite templates flag converted and CpG cytosines per strand", {
  x <- promoter("TCAGCGTTCC", "+")
  sense <- bisulfite_template(x, "sense")
  expect_equal(sense$assay_strand, "OT")
  expect_equal(sense$converted, "TTAGTGTTTT")
  expect_equal(which(sense$cpg_c), 5L)
  # The trailing C may pair with a G outside the window, so it is masked and not counted.
  expect_equal(which(sense$converted_c), c(2L, 9L))
  expect_equal(which(sense$cpg_mask), c(5L, 10L))
  antisense <- bisulfite_template(x, "antisense")
  expect_equal(antisense$assay_strand, "OB")
  expect_equal(antisense$original, "GGAACGCTGA")
  expect_equal(which(antisense$cpg_c), 5L)
  expect_equal(which(antisense$converted_c), 7L)
  expect_equal(bisulfite_template(promoter("TCAGCGTTCC", "-"), "sense")$assay_strand, "OB")
})

test_that("CpG islands are found around the fixture TSSs", {
  for (gene in c("BSPPLUS", "BSPMINUS")) {
    islands <- cpg_island_scan(fixture_promoter(gene))
    expect_true(any(islands$tss_start < 0L & islands$tss_end > 0L))
    expect_true(all(islands$gc_percent >= 50 & islands$cpg_obs_exp >= 0.6))
  }
  expect_equal(nrow(cpg_island_scan(promoter(strrep("AT", 200)))), 0L)
})
