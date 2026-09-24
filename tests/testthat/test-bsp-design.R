test_that("parameter defaults, input specs and validation agree", {
  defaults <- bsp_primer_params()
  expect_setequal(names(bsp_param_spec()), names(defaults))
  expect_identical(validate_bsp_params(list()), defaults)
  expect_equal(validate_bsp_params(list(tm_opt = 60))$tm_opt, 60)
  expect_error(validate_bsp_params(list(tm_min = 61, tm_opt = 60)), "min <= optimal <= max")
  expect_error(validate_bsp_params(list(primer_size_min = 22.5)), "whole number")
  expect_error(validate_bsp_params(list(nonsense = 1)), "Unknown BSP parameter")
  expect_error(validate_bsp_params(list(tm_opt = NA)), "must be a number")
  expect_error(validate_bsp_params(list(), sequence_length = 150L), "exceeds the extracted sequence")
})

test_that("targets resolve from TSS-relative bounds and CpG islands", {
  x <- fixture_promoter("BSPPLUS")
  default <- resolve_bsp_target(x)
  expect_equal(c(default$tss_start, default$tss_end), c(-50L, 50L))
  expect_equal(c(default$start, default$end), c(951L, 1050L))
  expect_equal(c(default$target_id, default$kind), c("R1", "range"))
  explicit <- resolve_bsp_target(x, c(100, -300))
  expect_equal(c(explicit$tss_start, explicit$tss_end), c(-300L, 100L))
  island <- resolve_bsp_target(x, "cpg_island")
  expect_equal(island$target_id, "CGI2")
  expect_lte(island$length, 265L)
  expect_true(island$tss_start < 0L && island$tss_end > 0L)
  expect_error(resolve_bsp_target(x, c(-2000, 0)), "0 is not a valid")
  expect_error(resolve_bsp_target(x, c(-2000, 10)), "outside the extracted sequence")
  expect_error(resolve_bsp_target(x, c(-500, 10)), "can be at most 406 bp")
})

test_that("each CpG island gets one window, centred on the TSS or the island", {
  x <- fixture_promoter("BSPPLUS", upstream = 1500L)
  islands <- cpg_island_scan(x)
  expect_equal(nrow(islands), 2L)
  targets <- island_targets(x)
  expect_equal(targets$target_id, c("CGI1", "CGI2"))
  expect_equal(targets$island, 1:2)
  expect_equal(targets$length, c(265L, 265L))
  # CGI1 lies upstream of the TSS: its window is centred on the island.
  expect_equal(targets$start[1L], (islands$start[1L] + islands$end[1L]) %/% 2L - 132L)
  # CGI2 contains the TSS (sense position 1501): its window is centred on it.
  expect_equal(c(targets$start[2L], targets$end[2L]), c(1501L - 132L, 1501L + 132L))
  expect_true(all(targets$start >= islands$start & targets$end <= islands$end))
  expect_equal(targets$n_cpg, vapply(seq_len(2L), function(i) sum(cpg_sites(x)$position >= targets$start[i] &
                                                                    cpg_sites(x)$position < targets$end[i]), integer(1)))
  # An island that fits inside the window is targeted whole.
  wide <- island_targets(x, list(product_size_min = 400L, product_size_max = 700L), islands = 1L)
  expect_equal(c(wide$start, wide$end), c(islands$start[1L], islands$end[1L]))
  expect_error(island_targets(x, islands = 3L), "CpG island 3 does not exist")
  minus <- island_targets(fixture_promoter("BSPMINUS", upstream = 1500L))
  expect_equal(c(minus$start, minus$end), c(1501L - 132L, 1501L + 132L))
})

test_that("non-island targets pick the densest CpG clusters outside islands", {
  x <- fixture_promoter("BSPPLUS", upstream = 1500L)
  islands <- cpg_island_scan(x)
  targets <- non_island_targets(x)
  expect_gt(nrow(targets), 0L)
  expect_lte(nrow(targets), bsp_primer_params()$fallback_targets)
  expect_true(all(targets$kind == "non_island"))
  expect_true(all(targets$n_cpg >= bsp_primer_params()$min_cpg_in_amplicon))
  expect_false(is.unsorted(rev(targets$n_cpg)))
  for (i in seq_len(nrow(targets))) {
    expect_false(any(targets$start[i] <= islands$end & targets$end[i] >= islands$start))
  }
  expect_equal(nrow(non_island_targets(x, list(fallback_targets = 1L))), 1L)
  # Equal clusters: the one nearer the TSS comes first.
  background <- strsplit(strrep("ATTAGATTCA", 100L), "")[[1L]]
  plant <- function(bases, at) {
    bases[at:(at + 7L)] <- strsplit("ACGACGCG", "")[[1L]]
    bases
  }
  bases <- plant(plant(background, 100L), 700L)
  synthetic <- promoter(paste(bases, collapse = ""), region_start = 1L, tss = 800L)
  expect_equal(nrow(cpg_island_scan(synthetic)), 0L)
  clusters <- non_island_targets(synthetic)
  expect_equal(clusters$target_id, c("NC1", "NC2"))
  expect_equal(clusters$start, c(701L, 101L))
  expect_equal(clusters$n_cpg, c(3L, 3L))
})

test_that("target tables are validated", {
  x <- fixture_promoter("BSPPLUS")
  params <- bsp_primer_params()
  table <- rbind(resolve_bsp_target(x), island_targets(x))
  expect_equal(bsp_targets(x, table, params)$target_id, c("R1", "CGI1", "CGI2"))
  expect_error(bsp_targets(x, rbind(table, table), params), "unique")
  bad <- table
  bad$end[1L] <- 99999L
  expect_error(bsp_targets(x, bad, params), "inside the extracted sequence")
  expect_error(bsp_targets(x, table[, c("target_id", "start")], params), "lacks column")
})

test_that("candidate primers obey the bisulfite rules", {
  x <- fixture_promoter("BSPMINUS")
  params <- bsp_primer_params()
  target <- resolve_bsp_target(x, NULL, params)
  for (strand in c("sense", "antisense")) {
    template <- bisulfite_template(x, strand)
    window <- sense_to_template(template, target$start, target$end)
    left <- bsp_candidates(template, window, params, "left")
    right <- bsp_candidates(template, window, params, "right")
    expect_gt(nrow(left), 0L)
    expect_gt(nrow(right), 0L)
    expect_true(all(left$end < window$start))
    expect_true(all(right$start > window$end))
    for (side in list(list(c = left, three = function(s, e) (e - 2L):e), list(c = right, three = function(s, e) s:(s + 2L)))) {
      candidates <- side$c
      for (i in seq_len(nrow(candidates))) {
        span <- candidates$start[i]:candidates$end[i]
        expect_false(any(template$cpg_mask[span]))
        expect_gte(sum(template$converted_c[span]), params$min_converted_c)
        expect_true(any(template$converted_c[side$three(candidates$start[i], candidates$end[i])]))
        expect_false(grepl("([ACG])\\1{5,}", substr(template$converted, candidates$start[i], candidates$end[i]), perl = TRUE))
      }
    }
    expect_equal(left$sequence, substring(template$converted, left$start, left$end))
    expect_equal(right$sequence, reverse_complement_text(substring(template$converted, right$start, right$end)))
  }
})

test_that("candidate filtering counts the footprints each bisulfite rule removes", {
  x <- fixture_promoter("BSPMINUS")
  params <- bsp_primer_params()
  target <- resolve_bsp_target(x, NULL, params)
  template <- bisulfite_template(x, "sense")
  window <- sense_to_template(template, target$start, target$end)
  left <- bsp_candidates(template, window, params, "left")
  rejected <- attr(left, "rejected")
  expect_named(rejected, unname(candidate_rules))
  first <- max(1L, window$end + params$primer_size_min - params$product_size_max + 1L)
  considered <- sum(vapply(params$primer_size_min:params$primer_size_max, function(size) max(0L, window$start - size - first + 1L), integer(1)))
  expect_equal(sum(rejected) + nrow(left), considered)
  expect_gt(rejected[["CpG under primer"]], 0L)
  # Rules that remove everything still return an empty table with the tallies.
  none <- bsp_candidates(template, window, utils::modifyList(params, list(min_converted_c = 30L)), "left")
  expect_equal(nrow(none), 0L)
  expect_equal(sum(attr(none, "rejected")), considered)
  expect_equal(format_rejections(c(a = 1, b = 5, c = 0, d = 3), top = 2L), "b 5, d 3")
  expect_equal(format_rejections(integer()), "")
})

test_that("pairs whose primer site repeats inside the amplicon are dropped", {
  # A tandem repeat: the right primer site (CCC...GGG) also occurs at the start of the amplicon.
  unit <- "GGTTTGGTTAATAGGTGTAGGGG"
  template <- list(converted_methylated = paste0("AAAA", unit, strrep("A", 200L), unit, "AAAA"))
  pairs <- data.frame(left_start = c(5L, 5L), left_end = c(26L, 18L), right_start = c(227L, 227L), right_end = c(249L, 249L))
  pairs$left_primer3_sequence <- substring(template$converted_methylated, pairs$left_start, pairs$left_end)
  expect_equal(count_matches(c("AAAA", "ATATA"), c("AA", "ATA")), c(3L, 2L))
  expect_equal(unique_primer_sites(template, pairs), c(FALSE, FALSE))
  template$converted_methylated <- paste0("AAAA", unit, strrep("A", 200L), "CCCCTACACCTATTAACCAAACC", "AAAA")
  expect_equal(unique_primer_sites(template, pairs), c(TRUE, TRUE))
})

test_that("hints name the main reason targets got no pair", {
  design <- list(
    targets = data.frame(target_id = c("CGI1", "CGI2"), n_pairs = c(0L, 5L)),
    diagnostics = data.frame(target_id = c("CGI1", "CGI2"), left_primer3 = c(0L, 10L), right_primer3 = c(3L, 10L), pairs_in_range = 0L),
    rejections = data.frame(target_id = c("CGI1", "CGI1", "CGI1", "CGI2"),
                            reason = c("CpG under primer", "GC out of range", "Tm too low", "too few converted C"),
                            n = c(60L, 30L, 10L, 500L)))
  hints <- design_hints(design)
  expect_length(hints, 2L)
  expect_match(hints[1L], "^60% of rejected candidate primers had a CpG under them\\. Allow 1-2 CpGs per primer")
  expect_match(hints[2L], "^40% of rejected candidate primers were too AT-rich")
  design$targets$n_pairs <- c(1L, 5L)
  expect_length(design_hints(design), 0L)
  design$targets$n_pairs <- c(0L, 0L)
  design$diagnostics$left_primer3 <- c(0L, 4L)
  expect_match(design_hints(design), "no pair fits the product size", all = FALSE)
})

test_that("allowed CpGs become degenerate bases in the 5' half", {
  x <- fixture_promoter("BSPPLUS")
  params <- validate_bsp_params(list(max_cpg_in_primer = 1L))
  target <- resolve_bsp_target(x, NULL, params)
  template <- bisulfite_template(x, "sense")
  left <- bsp_candidates(template, sense_to_template(template, target$start, target$end), params, "left")
  with_cpg <- left[left$cpg_in_primer == 1L, , drop = FALSE]
  expect_gt(nrow(with_cpg), 0L)
  expect_true(all(grepl("Y", with_cpg$sequence)))
  expect_true(all(regexpr("Y", with_cpg$sequence) <= with_cpg$length - with_cpg$length %/% 2L))
  expect_false(any(grepl("Y", with_cpg$primer3_sequence)))
})

test_that("diverse selection limits primer-site reuse", {
  params <- utils::modifyList(bsp_primer_params(), list(max_primer_reuse = 1L))
  rows <- diverse_rows(rep("sense", 4), c(100L, 102L, 100L, 200L), c(500L, 600L, 700L, 800L), c(1, 2, 3, 4), 10L, params)
  expect_equal(rows, c(1L, 4L))
})

test_that("designed pairs satisfy BSP invariants and round-trip through the FASTA", {
  skip_if_no_primer3()
  fasta <- fixture_reference("fa")
  for (gene in c("BSPPLUS", "BSPMINUS")) {
    x <- fixture_promoter(gene)
    design <- design_bsp_primers(x, "both")
    pairs <- design$pairs
    expect_s3_class(design, "bsp_design")
    expect_false(design$fallback)
    expect_gt(nrow(pairs), 0L)
    expect_equal(pairs$rank, seq_len(nrow(pairs)))
    expect_true(all(design$targets$kind == "island"))
    expect_setequal(unique(pairs$target_id), design$targets$target_id)
    expect_equal(design$targets$n_pairs, as.vector(table(factor(pairs$target_id, levels = design$targets$target_id))))
    expect_true(all(design$targets$n_pairs <= design$params$num_return))
    expect_setequal(unique(pairs$strand), c("sense", "antisense"))
    params <- design$params
    expect_true(all(pairs$product_size >= params$product_size_min & pairs$product_size <= params$product_size_max))
    expect_true(all(pairs$cpg_in_amplicon >= params$min_cpg_in_amplicon))
    expect_true(all(pairs$left_converted_c >= params$min_converted_c & pairs$right_converted_c >= params$min_converted_c))
    expect_true(all(pairs$left_cpg == 0L & pairs$right_cpg == 0L))
    expect_true(all(pairs$tm_diff <= params$max_tm_diff))
    for (i in seq_len(nrow(pairs))) {
      p <- pairs[i, ]
      bounds <- as.integer(strsplit(gsub("\\+", "", p$amplicon_tss), "..", fixed = TRUE)[[1L]])
      target <- design$targets[design$targets$target_id == p$target_id, ]
      expect_true(bounds[1L] <= target$tss_start && bounds[2L] >= target$tss_end)
      expect_equal(p$amplicon_end - p$amplicon_start + 1L, p$product_size)
      expect_equal(length(strsplit(p$cpg_tss_positions, ",")[[1L]]), p$cpg_in_amplicon)
      # Re-read each primer site from the genome, orient it to the assayed strand and convert it.
      oriented <- function(start, end) {
        genomic <- reference_sequence(fasta, p$scaffold, start, end)
        if (p$assay_strand == "OT") genomic else reverse_complement_text(genomic)
      }
      expect_equal(bisulfite_convert(oriented(p$left_genomic_start, p$left_genomic_end)), p$left_sequence)
      expect_equal(reverse_complement_text(bisulfite_convert(oriented(p$right_genomic_start, p$right_genomic_end))), p$right_sequence)
    }
    expect_equal(unique(pairs$assay_strand[pairs$strand == "sense"]), if (x$strand == "+") "OT" else "OB")
  }
})

test_that("an infeasible design reports diagnostics instead of failing", {
  skip_if_no_primer3()
  design <- design_bsp_primers(fixture_promoter("BSPPLUS"), "sense", params = list(min_converted_c = 15L))
  expect_null(design$pairs)
  expect_true(design$fallback)
  expect_true(all(design$targets$n_pairs == 0L))
  expect_setequal(design$diagnostics$target_id, design$targets$target_id)
  expect_true(all(design$diagnostics$left_candidates == 0L))
  expect_true(all(grepl("too few converted C", design$diagnostics$left_rejected)))
  expect_setequal(unique(design$rejections$target_id), design$targets$target_id)
  expect_output(print(design), "0 primer pairs")
  expect_output(print(design), "Hint: [0-9]+% of rejected candidate primers had too few non-CpG cytosines")
})

test_that("auto design falls back to non-island CpGs", {
  skip_if_no_primer3()
  fasta <- fixture_reference("fa")
  flank <- function(start, end) reference_sequence(fasta, "chrTest1", start, end)
  # No island at all.
  plain <- promoter(flank(1L, 1400L), region_start = 1L, tss = 700L)
  expect_equal(nrow(cpg_island_scan(plain)), 0L)
  design <- design_bsp_primers(plain, "both")
  expect_true(design$fallback)
  expect_true(all(design$targets$kind == "non_island"))
  expect_gt(nrow(design$pairs), 0L)
  expect_true(all(startsWith(design$pairs$target_id, "NC")))
  # An island too CpG-dense for any primer (a CpG every 4 bp) falls back to its flanks.
  dense <- promoter(paste0(flank(1L, 800L), strrep("ACGCTCGA", 125L), flank(801L, 1600L)), region_start = 1L, tss = 1300L)
  design <- design_bsp_primers(dense, "both")
  expect_equal(design$targets$kind[1L], "island")
  expect_equal(design$targets$n_pairs[1L], 0L)
  expect_true(design$fallback)
  expect_gt(nrow(design$pairs), 0L)
  expect_true(all(design$pairs$target_kind == "non_island"))
  expect_output(print(design), "non-island CpG targets were designed")
  # An explicit target never falls back.
  expect_false(design_bsp_primers(dense, "both", target = "cpg_island")$fallback)
})

test_that("design exports record parameters, oligos and amplicon references", {
  skip_if_no_primer3()
  design <- design_bsp_primers(fixture_promoter("BSPPLUS"), "sense")
  pairs_file <- tempfile(fileext = ".tsv")
  write_bsp_pairs(design, pairs_file)
  lines <- readLines(pairs_file)
  expect_match(lines[1L], "^# design: BSPPLUS chrTest1:")
  expect_equal(lines[2L], "# targets: CGI1=-510..-246; CGI2=-132..+133")
  expect_match(lines[3L], "^# params: product_size_min=200; ")
  expect_true(all(vapply(names(bsp_primer_params()), function(n) grepl(paste0(n, "="), lines[3L], fixed = TRUE), logical(1))))
  table <- utils::read.delim(pairs_file, comment.char = "#")
  expect_equal(nrow(table), nrow(design$pairs))

  order_file <- tempfile(fileext = ".tsv")
  write_order_sheet(design, order_file, ranks = 1L)
  sheet <- utils::read.delim(order_file)
  expect_equal(sheet$name, paste0("BSPPLUS_", design$pairs$target_id[1L], "_BSP1_", c("F", "R")))
  expect_equal(sheet$sequence, c(design$pairs$left_sequence[1L], design$pairs$right_sequence[1L]))

  amplicons <- bsp_amplicons(design, 1L)
  expect_length(amplicons, 3L)
  expect_true(all(startsWith(names(amplicons), paste0("BSPPLUS_", design$pairs$target_id[1L], "_BSP1|"))))
  expect_equal(unname(amplicons[2L]), bisulfite_convert(amplicons[[1L]]))
  expect_equal(unname(amplicons[3L]), bisulfite_convert(amplicons[[1L]], cpg = "methylated"))
  expect_equal(nchar(amplicons[[1L]]), design$pairs$product_size[1L])
  expect_true(startsWith(amplicons[[2L]], design$pairs$left_sequence[1L]))
  fasta_file <- tempfile(fileext = ".fa")
  write_amplicon_references(design, fasta_file, ranks = 1L)
  expect_equal(sum(startsWith(readLines(fasta_file), ">")), 3L)
  expect_error(write_order_sheet(design, order_file, ranks = 999L), "not in the design")
})
