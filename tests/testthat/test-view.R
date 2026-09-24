test_that("brush ranges round and clamp to sequence positions", {
  expect_equal(brush_to_sense(10.4, 20.6, 100L), c(10L, 21L))
  expect_equal(brush_to_sense(20.6, 10.4, 100L), c(10L, 21L))
  expect_equal(brush_to_sense(-5, 250, 100L), c(1L, 100L))
})

test_that("clicks select the CpG island under the cursor", {
  islands <- data.frame(start = c(100L, 400L), end = c(300L, 700L))
  expect_equal(island_at(islands, 150, 0), 1L)
  expect_equal(island_at(islands, 700.4, -0.2), 2L)
  expect_true(is.na(island_at(islands, 350, 0)))
  expect_true(is.na(island_at(islands, 150, 0.9)))
  expect_true(is.na(island_at(islands[0, ], 150, 0)))
})

test_that("selected regions report TSS-relative, sequence and genomic positions", {
  for (gene in c("BSPPLUS", "BSPMINUS")) {
    x <- fixture_promoter(gene)
    tss <- tss_sense_position(x)
    region <- describe_region(x, tss - 100L, tss + 49L)
    expect_equal(region$tss_relative, "-100..+50")
    expect_equal(region$sequence, sprintf("%s..%s", tss - 100L, tss + 49L))
    genomic <- sense_to_genomic(x, tss - 100L, tss + 49L)
    expect_equal(region$genomic, sprintf("%s:%s-%s (%s)", x$scaffold, genomic$genomic_start, genomic$genomic_end, x$strand))
    expect_equal(region$length, 150L)
    cpg <- cpg_sites(x)
    expect_equal(region$cpg, sum(cpg$tss_relative >= -100L & cpg$tss_relative < 50L))
    expect_match(region$cpg_islands, "CGI")
    expect_equal(genomic$genomic_end - genomic$genomic_start + 1L, 150L)
    # The TSS itself maps to the annotated TSS coordinate.
    cursor <- describe_position(x, tss)
    expect_match(cursor, "^Cursor: \\+1 ")
    expect_true(grepl(sprintf("%s:%s,", x$scaffold, format(x$tss, big.mark = ",")), cursor, fixed = TRUE))
  }
  x <- fixture_promoter("BSPPLUS")
  expect_equal(describe_region(x, 1L, 50L)$cpg_islands, "none")
  expect_error(describe_region(x, 10L, 5000L), "beyond the extracted sequence")
  expect_equal(describe_position(x, -3), "")
})

test_that("island table and choices label every island", {
  x <- fixture_promoter("BSPPLUS")
  table <- island_table(x)
  expect_equal(table$island, c("CGI1", "CGI2"))
  expect_equal(table$design_window, format_tss_interval(island_targets(x)$tss_start, island_targets(x)$tss_end))
  choices <- island_choices(table)
  expect_equal(unname(choices), c("1", "2"))
  expect_match(names(choices)[2L], "^CGI2: -151\\.\\.\\+440 \\(591 bp")
  expect_length(island_choices(table[0, ]), 0L)
})

test_that("the CpG map draws every layer", {
  x <- fixture_promoter("BSPMINUS")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot_cpg_map(x))
  targets <- rbind(island_targets(x), non_island_targets(x), resolve_bsp_target(x))
  expect_null(plot_cpg_map(x, targets = targets, selection = c(800L, 1100L), island = 1L))
  # The x axis is in sequence positions, so brush coordinates map directly.
  expect_equal(graphics::par("usr")[1:2], c(1, nchar(x$sequence)) + c(-1, 1) * 0.04 * (nchar(x$sequence) - 1))
})
