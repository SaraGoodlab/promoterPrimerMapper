# Bisulfite-sequencing PCR (BSP) primer design.

#' Default BSP primer design parameters.
#'
#' Sanger/cloning-oriented defaults. Bisulfite-converted DNA is fragmented and
#' has a three-letter alphabet, so primers are longer and lower in GC than
#' ordinary PCR primers, and amplicons are kept short.
#'
#' * `product_size_min`, `product_size_max`: amplicon length (bp).
#' * `primer_size_min`, `primer_size_opt`, `primer_size_max`: primer length (nt).
#' * `tm_min`, `tm_opt`, `tm_max`: primer melting temperature (Primer3 nearest-neighbour, C).
#' * `max_tm_diff`: maximum Tm difference within a pair.
#' * `gc_min`, `gc_max`: primer GC percentage.
#' * `max_poly_x`: longest homopolymer allowed by Primer3 (T runs are common after conversion).
#' * `max_poly_non_t`: longest run of one non-T base in the converted template footprint.
#' * `min_converted_c`: minimum non-CpG cytosines under each primer; these
#'   discriminate converted from unconverted DNA.
#' * `three_prime_window`: a converted cytosine must lie within this many
#'   bases of each primer's 3' end (0 disables the rule).
#' * `max_cpg_in_primer`: CpG cytosines allowed under a primer; allowed CpGs
#'   must sit in the 5' half and are written as degenerate `Y`/`R` bases.
#' * `min_cpg_in_amplicon`: minimum CpGs between the primers.
#' * `candidates_per_side`: best-scoring primers per side carried into pairing.
#' * `num_return`: number of primer pairs reported per target.
#' * `fallback_targets`: how many non-island CpG targets to design when no
#'   CpG island yields a primer pair (see [non_island_targets()]).
#' * `max_primer_reuse`: how many reported pairs may share one primer site
#'   (primers whose 3' ends are within half the minimum primer length).
#' @return A named list.
#' @export
bsp_primer_params <- function() {
  list(
    product_size_min = 200L, product_size_max = 450L,
    primer_size_min = 22L, primer_size_opt = 26L, primer_size_max = 30L,
    tm_min = 52, tm_opt = 58, tm_max = 64, max_tm_diff = 3,
    gc_min = 15, gc_max = 55,
    max_poly_x = 8L, max_poly_non_t = 5L,
    min_converted_c = 4L, three_prime_window = 3L,
    max_cpg_in_primer = 0L, min_cpg_in_amplicon = 3L,
    candidates_per_side = 200L, num_return = 10L, max_primer_reuse = 2L,
    fallback_targets = 3L
  )
}

# Input specification for each parameter: GUI group, label, bounds and step.
# Shared by validate_bsp_params() and the Shiny inputs so both apply the same limits.
bsp_param_spec <- function() {
  spec <- function(group, label, min, max, step, integer = TRUE, help = NULL) {
    list(group = group, label = label, min = min, max = max, step = step, integer = integer, help = help)
  }
  list(
    product_size_min = spec("Amplicon", "Min product size (bp)", 50, 5000, 10),
    product_size_max = spec("Amplicon", "Max product size (bp)", 50, 5000, 10),
    min_cpg_in_amplicon = spec("Amplicon", "Min CpGs in amplicon", 0, 500, 1),
    primer_size_min = spec("Primer size & Tm", "Min primer length (nt)", 15, 40, 1),
    primer_size_opt = spec("Primer size & Tm", "Optimal primer length (nt)", 15, 40, 1),
    primer_size_max = spec("Primer size & Tm", "Max primer length (nt)", 15, 40, 1),
    tm_min = spec("Primer size & Tm", "Min Tm (\u00b0C)", 40, 80, 0.5, FALSE),
    tm_opt = spec("Primer size & Tm", "Optimal Tm (\u00b0C)", 40, 80, 0.5, FALSE),
    tm_max = spec("Primer size & Tm", "Max Tm (\u00b0C)", 40, 80, 0.5, FALSE),
    max_tm_diff = spec("Primer size & Tm", "Max Tm difference (\u00b0C)", 0, 20, 0.5, FALSE),
    gc_min = spec("Composition", "Min GC (%)", 0, 100, 1, FALSE),
    gc_max = spec("Composition", "Max GC (%)", 0, 100, 1, FALSE),
    max_poly_x = spec("Composition", "Max homopolymer (Primer3)", 1, 20, 1),
    max_poly_non_t = spec("Composition", "Max non-T run", 1, 20, 1, help = "Longest run of A, C or G in the converted template under a primer."),
    min_converted_c = spec("Bisulfite rules", "Min converted C per primer", 0, 20, 1, help = "Non-CpG cytosines under the primer; selects for fully converted DNA."),
    three_prime_window = spec("Bisulfite rules", "Converted C within 3' window (nt)", 0, 10, 1, help = "0 disables the 3'-end rule."),
    max_cpg_in_primer = spec("Bisulfite rules", "Max CpGs per primer", 0, 3, 1, help = "Allowed CpGs must be in the 5' half and become degenerate Y/R bases."),
    candidates_per_side = spec("Search", "Candidates per side", 10, 2000, 10, help = "Best Primer3-scored primers per side carried into pairing."),
    num_return = spec("Search", "Pairs reported per target", 1, 200, 1),
    max_primer_reuse = spec("Search", "Max pairs sharing a primer site", 1, 50, 1, help = "Primers whose 3' ends are within half the min primer length count as one site."),
    fallback_targets = spec("Search", "Non-island targets (fallback)", 1, 20, 1, help = "Densest non-island CpG regions designed when no CpG island yields a primer pair.")
  )
}

#' Validate and complete BSP design parameters.
#' @param params Named list; missing entries are taken from [bsp_primer_params()].
#' @param sequence_length Optional extracted sequence length to check against.
#' @return The completed parameter list with integer fields coerced.
#' @export
validate_bsp_params <- function(params = list(), sequence_length = NULL) {
  defaults <- bsp_primer_params()
  params <- as.list(params)
  unknown <- setdiff(names(params), names(defaults))
  if (length(unknown)) stop(sprintf("Unknown BSP parameter(s): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
  params <- utils::modifyList(defaults, params)[names(defaults)]
  spec <- bsp_param_spec()
  problems <- character()
  for (name in names(params)) {
    value <- params[[name]]
    s <- spec[[name]]
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
      problems <- c(problems, sprintf("%s must be a number.", s$label))
      next
    }
    if (s$integer && value != round(value)) problems <- c(problems, sprintf("%s must be a whole number.", s$label))
    if (value < s$min || value > s$max) problems <- c(problems, sprintf("%s must be between %s and %s.", s$label, s$min, s$max))
    params[[name]] <- if (s$integer) as.integer(round(value)) else as.numeric(value)
  }
  if (!length(problems)) {
    check <- function(ok, message) if (!ok) problems <<- c(problems, message)
    check(params$primer_size_min <= params$primer_size_opt && params$primer_size_opt <= params$primer_size_max,
          "Primer lengths must satisfy min <= optimal <= max.")
    check(params$tm_min <= params$tm_opt && params$tm_opt <= params$tm_max, "Tm values must satisfy min <= optimal <= max.")
    check(params$gc_min <= params$gc_max, "Min GC must not exceed max GC.")
    check(params$product_size_min <= params$product_size_max, "Min product size must not exceed max product size.")
    check(params$product_size_max >= 2L * params$primer_size_min, "Max product size must fit two primers.")
    check(params$three_prime_window <= params$primer_size_min, "The 3' window must not exceed the min primer length.")
    check(params$min_converted_c <= params$primer_size_max, "Min converted C cannot exceed the max primer length.")
    if (!is.null(sequence_length)) {
      check(params$product_size_min <= sequence_length,
            sprintf("Min product size (%s bp) exceeds the extracted sequence (%s bp).", params$product_size_min, sequence_length))
    }
  }
  if (length(problems)) stop(paste(problems, collapse = " "), call. = FALSE)
  params
}

# Targets ----------------------------------------------------------------------
#
# A design covers one or more targets: sense intervals that every amplicon of
# the target must contain. Targets are kept as a data frame with one row per
# target: `target_id`, `kind` (island, range, non_island), `island` (index into
# cpg_island_scan(), NA otherwise), sense `start`/`end`, TSS-relative
# `tss_start`/`tss_end`, `length` and `n_cpg` (CpGs fully inside).

# Widest target that still leaves room for both primers in a mid-sized product.
target_window_width <- function(params) {
  max(1L, (params$product_size_min + params$product_size_max) %/% 2L - 2L * params$primer_size_max)
}

make_targets <- function(extracted, start, end, kind, id, island = NA_integer_) {
  start <- as.integer(start)
  end <- as.integer(end)
  tss <- tss_sense_position(extracted)
  cpg <- cpg_sites(extracted)$position
  data.frame(
    target_id = as.character(id), kind = rep_len(kind, length(start)), island = rep_len(as.integer(island), length(start)),
    start = start, end = end,
    tss_start = tss_relative_position(start, tss), tss_end = tss_relative_position(end, tss),
    length = end - start + 1L,
    n_cpg = vapply(seq_along(start), function(i) sum(cpg >= start[i] & cpg + 1L <= end[i]), integer(1)),
    stringsAsFactors = FALSE
  )
}

# Index of the island containing the TSS, else the nearest one.
nearest_island <- function(islands, tss) {
  if (!nrow(islands)) return(NA_integer_)
  distance <- ifelse(tss < islands$start, islands$start - tss, ifelse(tss > islands$end, tss - islands$end, 0L))
  which.min(distance)
}

# One target window per island: the whole island when it fits, otherwise a
# window centred on the TSS (when the island contains it) or on the island centre.
island_window <- function(start, end, tss, width) {
  if (end - start + 1L <= width) return(c(start, end))
  centre <- if (tss >= start && tss <= end) tss else (start + end) %/% 2L
  first <- min(max(start, centre - width %/% 2L), end - width + 1L)
  c(first, first + width - 1L)
}

#' BSP targets for CpG islands.
#'
#' One target per island predicted by [cpg_island_scan()]. An island no longer
#' than the target width (mid product size minus two maximum-length primers,
#' 265 bp with the default parameters) is targeted whole; a longer island gets a
#' window of that width centred on the TSS when the island contains it, and on
#' the island centre otherwise.
#' @param extracted A `promoter_sequence` object.
#' @param params Parameter list from [bsp_primer_params()].
#' @param islands Island indexes (rows of [cpg_island_scan()]); all by default.
#' @return A targets data frame (`target_id` `CGI<k>`, `kind`, `island`, sense
#'   `start`/`end`, `tss_start`/`tss_end`, `length`, `n_cpg`); no rows when the
#'   sequence has no island.
#' @export
island_targets <- function(extracted, params = bsp_primer_params(), islands = NULL) {
  validate_extracted_sequence(extracted)
  params <- validate_bsp_params(params)
  found <- cpg_island_scan(extracted)
  if (is.null(islands)) islands <- seq_len(nrow(found))
  islands <- suppressWarnings(as.integer(islands))
  bad <- islands[is.na(islands) | islands < 1L | islands > nrow(found)]
  if (length(bad)) {
    stop(sprintf("CpG island %s does not exist; the sequence has %s predicted island%s.",
                 paste(bad, collapse = ", "), nrow(found), if (nrow(found) == 1L) "" else "s"), call. = FALSE)
  }
  if (!length(islands)) return(make_targets(extracted, integer(), integer(), "island", character()))
  tss <- tss_sense_position(extracted)
  width <- target_window_width(params)
  bounds <- vapply(islands, function(k) island_window(found$start[k], found$end[k], tss, width), integer(2))
  make_targets(extracted, bounds[1L, ], bounds[2L, ], "island", paste0("CGI", islands), islands)
}

#' BSP targets for CpG sites outside CpG islands.
#'
#' Used when no CpG island yields a primer pair. Windows of the target width
#' (see [island_targets()]) start at each CpG outside the predicted islands and
#' stop before the next island; each is trimmed to its first and last CpG and
#' kept when it holds at least `min_cpg_in_amplicon` CpGs (at least one). The
#' densest windows are chosen greedily, ties going to the window closest to the
#' TSS, skipping windows that overlap one already chosen, up to
#' `fallback_targets`.
#' @param extracted A `promoter_sequence` object.
#' @param params Parameter list from [bsp_primer_params()].
#' @return A targets data frame (`target_id` `NC<k>`, densest first).
#' @export
non_island_targets <- function(extracted, params = bsp_primer_params()) {
  metadata <- validate_extracted_sequence(extracted)
  params <- validate_bsp_params(params)
  n <- metadata$sequence_length
  islands <- cpg_island_scan(extracted)
  in_island <- logical(n)
  for (i in seq_len(nrow(islands))) in_island[islands$start[i]:islands$end[i]] <- TRUE
  # Last base of the island-free stretch containing each position.
  runs <- rle(in_island)
  run_end <- rep(cumsum(runs$lengths), runs$lengths)
  cpg <- cpg_sites(extracted)$position
  cpg <- cpg[!in_island[cpg] & !in_island[cpg + 1L]]
  empty <- make_targets(extracted, integer(), integer(), "non_island", character())
  if (!length(cpg)) return(empty)
  limit <- pmin(cpg + target_window_width(params) - 1L, run_end[cpg])
  count <- vapply(seq_along(cpg), function(i) sum(cpg >= cpg[i] & cpg + 1L <= limit[i]), integer(1))
  last <- vapply(seq_along(cpg), function(i) max(cpg[cpg >= cpg[i] & cpg + 1L <= limit[i]]) + 1L, integer(1))
  keep <- count >= max(1L, params$min_cpg_in_amplicon)
  if (!any(keep)) return(empty)
  start <- cpg[keep]
  end <- last[keep]
  count <- count[keep]
  tss <- tss_sense_position(extracted)
  ordered <- order(-count, abs((start + end) / 2 - tss), start)
  chosen <- integer()
  for (i in ordered) {
    if (any(start[i] <= end[chosen] & end[i] >= start[chosen])) next
    chosen <- c(chosen, i)
    if (length(chosen) >= params$fallback_targets) break
  }
  make_targets(extracted, start[chosen], end[chosen], "non_island", paste0("NC", seq_along(chosen)))
}

#' Resolve a single BSP target region.
#' @param extracted A `promoter_sequence` object.
#' @param target `NULL` (TSS -50..+50), a numeric pair of TSS-relative bounds
#'   such as `c(-300, 100)`, or `"cpg_island"`: the window of the CpG island
#'   containing or nearest the TSS (see [island_targets()]).
#' @param params Parameter list from [bsp_primer_params()].
#' @return A one-row targets data frame (see [island_targets()]); a range has
#'   `target_id` `R1`.
#' @export
resolve_bsp_target <- function(extracted, target = NULL, params = bsp_primer_params()) {
  metadata <- validate_extracted_sequence(extracted)
  params <- validate_bsp_params(params)
  tss <- tss_sense_position(extracted)
  if (is.null(target)) target <- c(-50L, 50L)
  if (is.character(target)) {
    if (!identical(target, "cpg_island")) stop("target must be NULL, two TSS-relative positions, or 'cpg_island'.", call. = FALSE)
    k <- nearest_island(cpg_island_scan(extracted), tss)
    if (is.na(k)) stop("No CpG island was predicted in the extracted sequence; give TSS-relative target positions instead.", call. = FALSE)
    return(island_targets(extracted, params, k))
  }
  if (length(target) != 2L || anyNA(target) || !is.numeric(target)) stop("target must be two TSS-relative positions.", call. = FALSE)
  target <- sort(as.integer(target))
  start <- sense_position_from_tss(target[1L], tss)
  end <- sense_position_from_tss(target[2L], tss)
  label <- format_tss_interval(target[1L], target[2L])
  if (start < 1L || end > metadata$sequence_length) {
    stop(sprintf("Target %s lies outside the extracted sequence; extract a larger window.", label), call. = FALSE)
  }
  widest <- params$product_size_max - 2L * params$primer_size_min
  if (end - start + 1L > widest) {
    stop(sprintf("Target %s is %s bp; with a max product of %s bp and primers of at least %s nt it can be at most %s bp. Narrow the range or raise the max product size.",
                 label, end - start + 1L, params$product_size_max, params$primer_size_min, widest), call. = FALSE)
  }
  make_targets(extracted, start, end, "range", "R1")
}

# Targets for design_bsp_primers(): "auto" (island targets; the caller adds the
# non-island fallback), "non_island", a targets data frame, or anything
# resolve_bsp_target() accepts.
bsp_targets <- function(extracted, target, params) {
  if (is.data.frame(target)) {
    required <- c("target_id", "kind", "start", "end")
    missing <- setdiff(required, names(target))
    if (length(missing)) stop(sprintf("Target table lacks column(s): %s.", paste(missing, collapse = ", ")), call. = FALSE)
    n <- validate_extracted_sequence(extracted)$sequence_length
    if (anyNA(target$start) || anyNA(target$end) || any(target$start < 1L | target$end > n | target$start > target$end)) {
      stop("Target intervals must lie inside the extracted sequence.", call. = FALSE)
    }
    if (anyDuplicated(target$target_id)) stop("Target IDs must be unique.", call. = FALSE)
    island <- if (is.null(target$island)) NA_integer_ else target$island
    return(make_targets(extracted, target$start, target$end, target$kind, target$target_id, island))
  }
  if (identical(target, "auto")) return(island_targets(extracted, params))
  if (identical(target, "non_island")) return(non_island_targets(extracted, params))
  resolve_bsp_target(extracted, target, params)
}

window_counts <- function(flags, start, end) {
  cs <- c(0L, cumsum(flags))
  cs[end + 1L] - cs[start]
}

degenerate_footprint <- function(template, start, end) {
  bases <- split_bases(template$converted_methylated)
  bases[template$cpg_c] <- "Y"
  substring(paste(bases, collapse = ""), start, end)
}

# Rejection reasons tallied while filtering candidates, in the order the
# bisulfite rules are applied (each footprint counts once, under the first rule
# it fails).
candidate_rules <- c(unknown = "N or ambiguous base", cpg = "CpG under primer", converted = "too few converted C",
                     three_prime = "no converted C at 3' end", run = "long non-T run")

# Enumerate primer footprints (template coordinates) that satisfy the
# bisulfite-specific rules. Left footprints end before the target, right
# footprints start after it, both within reach of the maximum product size.
# Returns a data frame (possibly empty) whose "rejected" attribute counts the
# footprints removed by each rule (names from candidate_rules).
bsp_candidates <- function(template, target, params, side = c("left", "right")) {
  side <- match.arg(side)
  n <- template$length
  bases <- split_bases(template$original)
  unknown <- !bases %in% c("A", "C", "G", "T")
  rejected <- stats::setNames(integer(length(candidate_rules)), unname(candidate_rules))
  out <- lapply(params$primer_size_min:params$primer_size_max, function(size) {
    if (side == "left") {
      first <- max(1L, target$end + params$primer_size_min - params$product_size_max + 1L)
      last <- target$start - size
    } else {
      first <- target$end + 1L
      last <- min(n, target$start - params$primer_size_min + params$product_size_max - 1L) - size + 1L
    }
    if (last < first) return(NULL)
    start <- first:last
    end <- start + size - 1L
    # Apply one rule at a time and record how many footprints it removes.
    apply_rule <- function(rule, pass) {
      rejected[[candidate_rules[[rule]]]] <<- rejected[[candidate_rules[[rule]]]] + sum(!pass)
      start <<- start[pass]
      end <<- end[pass]
    }
    apply_rule("unknown", window_counts(unknown, start, end) == 0L)
    if (params$max_cpg_in_primer == 0L) {
      apply_rule("cpg", window_counts(template$cpg_mask, start, end) == 0L)
    } else {
      half <- size %/% 2L
      three_prime_half <- if (side == "left") window_counts(template$cpg_mask, end - half + 1L, end) else window_counts(template$cpg_mask, start, start + half - 1L)
      edge <- window_counts(template$cpg_mask & !template$cpg_c, start, end)
      apply_rule("cpg", window_counts(template$cpg_c, start, end) <= params$max_cpg_in_primer & three_prime_half == 0L & edge == 0L)
    }
    apply_rule("converted", window_counts(template$converted_c, start, end) >= params$min_converted_c)
    if (params$three_prime_window > 0L) {
      w <- params$three_prime_window
      apply_rule("three_prime", if (side == "left") window_counts(template$converted_c, end - w + 1L, end) > 0L
                                else window_counts(template$converted_c, start, start + w - 1L) > 0L)
    }
    if (!length(start)) return(NULL)
    run_pattern <- sprintf("([ACG])\\1{%s,}", params$max_poly_non_t)
    apply_rule("run", !grepl(run_pattern, substring(template$converted_methylated, start, end), perl = TRUE))
    if (!length(start)) return(NULL)
    footprint <- substring(template$converted_methylated, start, end)
    degenerate <- degenerate_footprint(template, start, end)
    data.frame(
      start = start, end = end, length = end - start + 1L,
      primer3_sequence = if (side == "left") footprint else reverse_complement_text(footprint),
      sequence = if (side == "left") degenerate else reverse_complement_text(degenerate),
      converted_c = window_counts(template$converted_c, start, end),
      cpg_in_primer = window_counts(template$cpg_c, start, end),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    out <- data.frame(start = integer(), end = integer(), length = integer(), primer3_sequence = character(),
                      sequence = character(), converted_c = integer(), cpg_in_primer = integer(), stringsAsFactors = FALSE)
  }
  attr(out, "rejected") <- rejected
  out
}

# Score candidate oligos with Primer3 and keep the best per side. The
# "rejected" attribute adds Primer3's rejection reasons to those of the
# bisulfite rules.
score_candidates <- function(candidates, side, params, primer3) {
  rejected <- attr(candidates, "rejected")
  if (is.null(candidates) || !nrow(candidates)) return(structure(candidates[0, , drop = FALSE], rejected = rejected))
  scored <- primer3_check_oligos(candidates$primer3_sequence, side, params, primer3)
  rejected <- add_counts(rejected, attr(scored, "rejected"))
  if (is.null(scored) || !nrow(scored)) return(structure(candidates[0, , drop = FALSE], rejected = rejected))
  candidates <- cbind(candidates[scored$candidate, , drop = FALSE], scored[, setdiff(names(scored), "candidate"), drop = FALSE])
  candidates <- candidates[order(candidates$penalty, candidates$start), , drop = FALSE]
  structure(utils::head(candidates, params$candidates_per_side), rejected = rejected)
}

# Sum named count vectors.
add_counts <- function(a, b) {
  if (!length(b)) return(a)
  if (!length(a)) return(b)
  names <- union(names(a), names(b))
  stats::setNames(vapply(names, function(n) sum(a[n], b[n], na.rm = TRUE), numeric(1)), names)
}

# "CpG under primer 812, GC out of range 230" from rejection counts (largest first).
format_rejections <- function(counts, top = 3L) {
  counts <- counts[!is.na(counts) & counts > 0]
  if (!length(counts)) return("")
  counts <- utils::head(sort(counts, decreasing = TRUE), top)
  paste(sprintf("%s %s", names(counts), counts), collapse = ", ")
}

prefix_columns <- function(x, prefix) stats::setNames(x, paste0(prefix, names(x)))

# Combine left and right candidates into pairs that satisfy product size, Tm
# difference and amplicon CpG content.
pair_candidates <- function(template, left, right, params) {
  if (is.null(left) || is.null(right)) return(NULL)
  grid <- expand.grid(l = seq_len(nrow(left)), r = seq_len(nrow(right)))
  product <- right$end[grid$r] - left$start[grid$l] + 1L
  tm_diff <- abs(left$tm[grid$l] - right$tm[grid$r])
  ordered <- right$start[grid$r] > left$end[grid$l]
  inner_cpg <- ifelse(ordered, window_counts(template$cpg_c, left$end[grid$l] + 1L, pmax(right$start[grid$r] - 1L, left$end[grid$l])), 0L)
  keep <- ordered & product >= params$product_size_min & product <= params$product_size_max &
    tm_diff <= params$max_tm_diff & inner_cpg >= params$min_cpg_in_amplicon
  if (!any(keep)) return(NULL)
  grid <- grid[keep, , drop = FALSE]
  pairs <- cbind(
    prefix_columns(left[grid$l, , drop = FALSE], "left_"),
    prefix_columns(right[grid$r, , drop = FALSE], "right_"),
    data.frame(product_size = product[keep], tm_diff = round(tm_diff[keep], 2), cpg_in_amplicon = inner_cpg[keep])
  )
  rownames(pairs) <- NULL
  pairs$provisional <- pairs$left_penalty + pairs$right_penalty + pairs$tm_diff
  pairs[order(pairs$provisional), , drop = FALSE]
}

format_tss_interval <- function(a, b) sprintf("%s..%s", format_tss_relative(a), format_tss_relative(b))

# Occurrences of `pattern` in each `text`, counting overlapping matches.
count_matches <- function(text, pattern) {
  vapply(seq_along(text), function(i) {
    hits <- gregexpr(sprintf("(?=%s)", pattern[i]), text[i], perl = TRUE)[[1L]]
    sum(hits > 0L)
  }, integer(1))
}

# Pairs whose primer sites each occur exactly once in their converted amplicon.
# A repeated site (tandem repeats) would prime inside the amplicon, and makes
# primer3_core abort when it locates the primers in the template.
unique_primer_sites <- function(template, pairs) {
  amplicons <- substring(template$converted_methylated, pairs$left_start, pairs$right_end)
  right_sites <- substring(template$converted_methylated, pairs$right_start, pairs$right_end)
  count_matches(amplicons, pairs$left_primer3_sequence) == 1L & count_matches(amplicons, right_sites) == 1L
}

# Add sense, TSS-relative and genomic coordinates to template-coordinate pairs.
annotate_pairs <- function(extracted, template, pairs) {
  tss <- tss_sense_position(extracted)
  scaffold <- get_field(extracted, "scaffold")
  locate <- function(start, end) {
    sense <- template_to_sense(template, start, end)
    genomic <- sense_to_genomic(extracted, sense$start, sense$end)
    list(genomic = genomic,
         interval = sprintf("%s:%s-%s", scaffold, genomic$genomic_start, genomic$genomic_end),
         tss = format_tss_interval(tss_relative_position(sense$start, tss), tss_relative_position(sense$end, tss)))
  }
  left <- locate(pairs$left_start, pairs$left_end)
  right <- locate(pairs$right_start, pairs$right_end)
  amplicon <- locate(pairs$left_start, pairs$right_end)
  cpg_template <- which(template$cpg_c)
  cpg_positions <- vapply(seq_len(nrow(pairs)), function(i) {
    p <- cpg_template[cpg_template > pairs$left_end[i] & cpg_template < pairs$right_start[i]]
    sense_c <- if (template$strand == "sense") p else template$length - p
    paste(format_tss_relative(sort(tss_relative_position(sense_c, tss))), collapse = ",")
  }, character(1))
  data.frame(
    strand = template$strand, assay_strand = template$assay_strand,
    left_sequence = pairs$left_sequence, right_sequence = pairs$right_sequence,
    product_size = pairs$product_size, cpg_in_amplicon = pairs$cpg_in_amplicon,
    left_tm = round(pairs$left_tm, 1), right_tm = round(pairs$right_tm, 1), tm_diff = pairs$tm_diff,
    left_gc = round(pairs$left_gc_percent, 1), right_gc = round(pairs$right_gc_percent, 1),
    left_length = pairs$left_length, right_length = pairs$right_length,
    left_converted_c = pairs$left_converted_c, right_converted_c = pairs$right_converted_c,
    left_cpg = pairs$left_cpg_in_primer, right_cpg = pairs$right_cpg_in_primer,
    amplicon_interval = amplicon$interval, amplicon_tss = amplicon$tss,
    left_interval = left$interval, right_interval = right$interval,
    left_tss = left$tss, right_tss = right$tss,
    cpg_tss_positions = cpg_positions,
    score = round(pairs$pair_penalty + pairs$tm_diff, 3), pair_penalty = round(pairs$pair_penalty, 3),
    left_penalty = round(pairs$left_penalty, 3), right_penalty = round(pairs$right_penalty, 3),
    compl_any = pairs$compl_any, compl_end = pairs$compl_end,
    left_self_any = pairs$left_self_any, left_self_end = pairs$left_self_end, left_hairpin = pairs$left_hairpin,
    right_self_any = pairs$right_self_any, right_self_end = pairs$right_self_end, right_hairpin = pairs$right_hairpin,
    scaffold = scaffold,
    amplicon_start = amplicon$genomic$genomic_start, amplicon_end = amplicon$genomic$genomic_end,
    left_genomic_start = left$genomic$genomic_start, left_genomic_end = left$genomic$genomic_end,
    right_genomic_start = right$genomic$genomic_start, right_genomic_end = right$genomic$genomic_end,
    template_left_start = pairs$left_start, template_left_end = pairs$left_end,
    template_right_start = pairs$right_start, template_right_end = pairs$right_end,
    stringsAsFactors = FALSE
  )
}

# Design pairs for one target (a row of a targets table) on one template strand.
design_strand <- function(extracted, template, target, params, primer3) {
  window <- sense_to_template(template, target$start, target$end)
  counts <- c(left_candidates = 0L, right_candidates = 0L, left_primer3 = 0L, right_primer3 = 0L, pairs_in_range = 0L, pairs_checked = 0L)
  rejected <- c(left = "", right = "")
  reasons <- list(left = NULL, right = NULL)
  finish <- function(pairs = NULL) {
    rejections <- do.call(rbind, lapply(c("left", "right"), function(side) {
      n <- reasons[[side]]
      n <- n[!is.na(n) & n > 0]
      if (!length(n)) return(NULL)
      data.frame(target_id = target$target_id, strand = template$strand, side = side, reason = names(n), n = as.integer(n),
                 stringsAsFactors = FALSE)
    }))
    if (!is.null(pairs)) {
      pairs <- cbind(data.frame(target_id = target$target_id, target_kind = target$kind,
                                target_tss = format_tss_interval(target$tss_start, target$tss_end), stringsAsFactors = FALSE),
                     pairs)
    }
    list(pairs = pairs, diagnostics = data.frame(target_id = target$target_id, strand = template$strand,
                                                 assay_strand = template$assay_strand, t(counts),
                                                 left_rejected = rejected[["left"]], right_rejected = rejected[["right"]],
                                                 row.names = NULL, stringsAsFactors = FALSE),
         rejections = rejections)
  }
  left <- bsp_candidates(template, window, params, "left")
  right <- bsp_candidates(template, window, params, "right")
  counts[["left_candidates"]] <- nrow(left)
  counts[["right_candidates"]] <- nrow(right)
  # Score both sides even when one is empty, so the rejection reasons are complete.
  left <- score_candidates(left, "LEFT", params, primer3)
  right <- score_candidates(right, "RIGHT", params, primer3)
  reasons <- list(left = attr(left, "rejected"), right = attr(right, "rejected"))
  rejected <- c(left = format_rejections(reasons$left), right = format_rejections(reasons$right))
  counts[["left_primer3"]] <- nrow(left)
  counts[["right_primer3"]] <- nrow(right)
  if (!nrow(left) || !nrow(right)) return(finish())
  pairs <- pair_candidates(template, left, right, params)
  if (is.null(pairs)) return(finish())
  pairs <- pairs[unique_primer_sites(template, pairs), , drop = FALSE]
  if (!nrow(pairs)) return(finish())
  counts[["pairs_in_range"]] <- nrow(pairs)
  # Pre-select a diverse shortlist so the Primer3 pair check is not spent on shifted duplicates.
  shortlist <- diverse_rows(rep(template$strand, nrow(pairs)), pairs$left_end, pairs$right_start, pairs$provisional,
                            max(50L, 5L * params$num_return), params)
  pairs <- pairs[shortlist, , drop = FALSE]
  amplicons <- substring(template$converted_methylated, pairs$left_start, pairs$right_end)
  checked <- primer3_check_pairs(amplicons, pairs$left_primer3_sequence, pairs$right_primer3_sequence, params, primer3)
  if (is.null(checked) || !nrow(checked)) return(finish())
  counts[["pairs_checked"]] <- nrow(checked)
  pairs <- cbind(pairs[checked$candidate, , drop = FALSE], checked[, setdiff(names(checked), "candidate"), drop = FALSE])
  finish(annotate_pairs(extracted, template, pairs))
}

# Greedy diverse selection: walk pairs in score order and skip a pair when its
# left or right primer site is already used by `max_primer_reuse` chosen pairs.
# Primers whose 3' ends lie within half the minimum primer length of each other
# (same strand and side) count as one site, so shifted variants of one primer
# do not crowd out alternatives. Returns row indices in score order.
diverse_rows <- function(strand, left_three_prime, right_three_prime, score, n, params) {
  proximity <- max(1L, params$primer_size_min %/% 2L)
  ordered <- order(score)
  chosen <- integer()
  for (i in ordered) {
    same <- chosen[strand[chosen] == strand[i]]
    if (sum(abs(left_three_prime[same] - left_three_prime[i]) <= proximity) >= params$max_primer_reuse ||
        sum(abs(right_three_prime[same] - right_three_prime[i]) <= proximity) >= params$max_primer_reuse) next
    chosen <- c(chosen, i)
    if (length(chosen) >= n) break
  }
  chosen
}

# Up to num_return diverse pairs per target, targets in table order; ranks are
# unique across the whole design.
select_pairs <- function(pairs, params, target_ids = unique(pairs$target_id)) {
  if (is.null(pairs) || !nrow(pairs)) return(NULL)
  rows <- unlist(lapply(target_ids, function(id) {
    index <- which(pairs$target_id == id)
    index[diverse_rows(pairs$strand[index], pairs$template_left_end[index], pairs$template_right_start[index],
                       pairs$score[index], params$num_return, params)]
  }))
  pairs <- pairs[rows, , drop = FALSE]
  rownames(pairs) <- NULL
  cbind(rank = seq_len(nrow(pairs)), pairs)
}

# Design every target on every template; returns selected pairs and diagnostics.
design_targets <- function(extracted, templates, targets, params, primer3, progress = NULL) {
  steps <- nrow(targets) * length(templates)
  step <- 0L
  runs <- list()
  for (i in seq_len(nrow(targets))) {
    for (strand in names(templates)) {
      step <- step + 1L
      if (is.function(progress)) progress(step, steps, sprintf("%s, %s strand", targets$target_id[i], strand))
      runs[[length(runs) + 1L]] <- design_strand(extracted, templates[[strand]], targets[i, , drop = FALSE], params, primer3)
    }
  }
  pairs <- do.call(rbind, lapply(runs, `[[`, "pairs"))
  list(pairs = select_pairs(pairs, params, targets$target_id),
       diagnostics = do.call(rbind, lapply(runs, `[[`, "diagnostics")),
       rejections = do.call(rbind, lapply(runs, `[[`, "rejections")))
}

# Rejection reasons grouped by the sequence problem behind them, with advice.
rejection_groups <- list(
  cpg = list(reasons = "CpG under primer", label = "had a CpG under them",
             advice = "Allow 1-2 CpGs per primer ('Max CpGs per primer'; they become degenerate Y/R bases in the 5' half) and/or raise the max product size so primers can reach CpG-poor flanks."),
  at_rich = list(reasons = c("GC out of range", "Tm too low"), label = "were too AT-rich after conversion (GC or Tm too low)",
                 advice = "This is typical of C-rich, G-poor sequence, where nearly every base converts to T. Lower 'Min GC' and 'Min Tm', allow longer primers, or choose another region."),
  no_c = list(reasons = c("too few converted C", "no converted C at 3' end"), label = "had too few non-CpG cytosines to convert",
              advice = "This is typical of G-rich, C-poor sequence. Lower 'Min converted C per primer' or the 3' window, or choose another region."),
  runs = list(reasons = c("long non-T run", "homopolymer"), label = "had long runs of one base",
              advice = "Raise 'Max non-T run' or 'Max homopolymer', or choose another region.")
)

# Advice for targets without pairs, from the main groups of rejection reasons.
design_hints <- function(design) {
  failed <- design$targets$target_id[design$targets$n_pairs == 0L]
  if (!length(failed) || is.null(design$diagnostics)) return(character())
  hints <- character()
  rejections <- design$rejections
  if (!is.null(rejections)) rejections <- rejections[rejections$target_id %in% failed, , drop = FALSE]
  if (!is.null(rejections) && nrow(rejections)) {
    group <- vapply(rejections$reason, function(reason) {
      hit <- names(rejection_groups)[vapply(rejection_groups, function(g) reason %in% g$reasons, logical(1))]
      if (length(hit)) hit else "other"
    }, character(1), USE.NAMES = FALSE)
    share <- sort(tapply(rejections$n, group, sum) / sum(rejections$n), decreasing = TRUE)
    # The main problem, plus any other that removes at least a quarter of the sites.
    main <- names(share)[seq_along(share) == 1L | share >= 0.25]
    main <- setdiff(main, "other")
    hints <- vapply(main, function(g) sprintf("%s%% of rejected candidate primers %s. %s", round(100 * share[[g]]),
                                             rejection_groups[[g]]$label, rejection_groups[[g]]$advice), character(1), USE.NAMES = FALSE)
  }
  diagnostics <- design$diagnostics[design$diagnostics$target_id %in% failed, , drop = FALSE]
  if (any(diagnostics$left_primer3 > 0L & diagnostics$right_primer3 > 0L & diagnostics$pairs_in_range == 0L)) {
    hints <- c(hints, "Where primers exist on both sides, no pair fits the product size, Tm difference and minimum CpG count; widen the product size range or the max Tm difference.")
  }
  hints
}

#' Design bisulfite-sequencing PCR primers for an extracted promoter.
#'
#' Candidate primers are enumerated on the bisulfite-converted template of the
#' chosen strand(s) and filtered by bisulfite rules (no methylation-dependent
#' CpG cytosine under the primer, enough converted cytosines, a converted
#' cytosine near the 3' end, no long non-T runs). Primer3 then scores each
#' oligo (Tm, GC, hairpins, self-dimers); pairs are formed around each target
#' and re-checked by Primer3 for cross-dimers.
#'
#' By default (`target = "auto"`) every predicted CpG island is a target (see
#' [island_targets()]). When no island yields a primer pair, or the sequence
#' has no island, the densest non-island CpG regions are designed instead (see
#' [non_island_targets()]) and `fallback` is `TRUE`.
#' @param extracted A `promoter_sequence` object from [extract_tss_region()].
#' @param strand Template strand: `sense` (transcription orientation),
#'   `antisense`, or `both`.
#' @param target `"auto"`, `"non_island"`, a targets data frame from
#'   [island_targets()], [non_island_targets()] or [resolve_bsp_target()], or
#'   anything [resolve_bsp_target()] accepts (a TSS-relative pair such as
#'   `c(-300, 100)`, or `"cpg_island"`).
#' @param params Parameter list from [bsp_primer_params()]; `num_return` pairs
#'   are reported per target.
#' @param primer3 Path to `primer3_core`.
#' @param progress Optional function called as `progress(step, steps, detail)`
#'   before each target and strand is designed.
#' @return A `bsp_design` object: `pairs` (data frame ranked within targets,
#'   `NULL` when no pair passed), `targets` (targets table with `n_pairs`),
#'   `fallback`, `diagnostics` (candidate counts per target, strand and stage,
#'   with the main rejection reasons per side), `rejections` (counts of
#'   rejected primer sites per target, strand, side and reason), `templates`,
#'   `params`, `extracted`, `cpg` and `islands`.
#' @export
design_bsp_primers <- function(extracted, strand = c("sense", "antisense", "both"), target = "auto",
                               params = bsp_primer_params(), primer3 = primer3_path(), progress = NULL) {
  strand <- match.arg(strand)
  metadata <- validate_extracted_sequence(extracted)
  params <- validate_bsp_params(params, metadata$sequence_length)
  targets <- bsp_targets(extracted, target, params)
  strands <- if (strand == "both") c("sense", "antisense") else strand
  templates <- stats::setNames(lapply(strands, function(s) bisulfite_template(extracted, s)), strands)
  run <- design_targets(extracted, templates, targets, params, primer3, progress)
  fallback <- identical(target, "auto") && is.null(run$pairs)
  if (fallback) {
    extra <- non_island_targets(extracted, params)
    if (nrow(extra)) {
      second <- design_targets(extracted, templates, extra, params, primer3, progress)
      run <- list(pairs = second$pairs, diagnostics = rbind(run$diagnostics, second$diagnostics),
                  rejections = rbind(run$rejections, second$rejections))
      targets <- rbind(targets, extra)
    }
  }
  targets$n_pairs <- vapply(targets$target_id, function(id) sum(run$pairs$target_id == id), integer(1), USE.NAMES = FALSE)
  structure(list(
    pairs = run$pairs, targets = targets, fallback = fallback, diagnostics = run$diagnostics, rejections = run$rejections,
    templates = templates, params = params, extracted = extracted,
    cpg = cpg_sites(extracted), islands = cpg_island_scan(extracted)
  ), class = c("bsp_design", "list"))
}

# One line per target, e.g. "CGI2 (island, -132..+133): 10 pairs".
target_lines <- function(targets) {
  pairs <- if (is.null(targets$n_pairs)) "" else sprintf(": %s pair%s", targets$n_pairs, ifelse(targets$n_pairs == 1L, "", "s"))
  sprintf("%s (%s, %s)%s", targets$target_id, gsub("_", "-", targets$kind), format_tss_interval(targets$tss_start, targets$tss_end), pairs)
}

#' @export
print.bsp_design <- function(x, ...) {
  n <- if (is.null(x$pairs)) 0L else nrow(x$pairs)
  cat(sprintf("BSP design for %s: %s target%s, %s primer pair%s\n", design_label(x), nrow(x$targets),
              if (nrow(x$targets) == 1L) "" else "s", n, if (n == 1L) "" else "s"))
  if (nrow(x$targets)) cat(paste0("  ", target_lines(x$targets), "\n"), sep = "")
  if (isTRUE(x$fallback)) cat("No CpG island yielded a primer pair; non-island CpG targets were designed.\n")
  if (n) print(utils::head(x$pairs[, c("rank", "target_id", "strand", "left_sequence", "right_sequence", "product_size", "cpg_in_amplicon", "score")], 10L))
  else print(x$diagnostics)
  hints <- design_hints(x)
  for (hint in hints) cat(strwrap(hint, initial = "Hint: ", prefix = "      "), sep = "\n")
  invisible(x)
}

design_label <- function(design) {
  extracted <- design$extracted
  label <- get_field(extracted, "gene_name", NA_character_)
  if (is.na(label) || !nzchar(label)) label <- get_field(extracted, "gene_id", NA_character_)
  if (is.na(label) || !nzchar(label)) label <- "BSP"
  gsub("[^A-Za-z0-9._-]+", "_", label)
}

# Oligo/amplicon names such as GENE_CGI1_BSP3.
pair_names <- function(design, pairs) sprintf("%s_%s_BSP%s", design_label(design), pairs$target_id, pairs$rank)

design_pairs <- function(design, ranks = NULL) {
  if (!inherits(design, "bsp_design")) stop("design must be a bsp_design object from design_bsp_primers().", call. = FALSE)
  pairs <- design$pairs
  if (is.null(pairs) || !nrow(pairs)) stop("The design contains no primer pairs.", call. = FALSE)
  if (is.null(ranks)) return(pairs)
  missing <- setdiff(ranks, pairs$rank)
  if (length(missing)) stop(sprintf("Pair rank(s) not in the design: %s.", paste(missing, collapse = ", ")), call. = FALSE)
  pairs[match(ranks, pairs$rank), , drop = FALSE]
}

bsp_param_comments <- function(design) {
  extracted <- design$extracted
  targets <- design$targets
  c(sprintf("design: %s %s:%s-%s strand=%s TSS=%s", design_label(design), get_field(extracted, "scaffold"),
            get_field(extracted, "region_start"), get_field(extracted, "region_end"), get_field(extracted, "strand"),
            get_field(extracted, "tss")),
    paste0("targets: ", paste(sprintf("%s=%s", targets$target_id, format_tss_interval(targets$tss_start, targets$tss_end)), collapse = "; "),
           if (isTRUE(design$fallback)) " (no CpG island yielded a pair; non-island fallback)" else ""),
    paste0("params: ", paste(names(design$params), unlist(design$params), sep = "=", collapse = "; ")))
}

#' Write BSP primer pairs as a delimited table.
#'
#' The table is preceded by `#` comment lines recording the design region and
#' every parameter used, so the design can be reproduced.
#' @param design A `bsp_design` object.
#' @param file Output path.
#' @param sep Field separator.
#' @return The output path, invisibly.
#' @export
write_bsp_pairs <- function(design, file, sep = "\t") {
  write_primer_results(design_pairs(design), file, sep = sep, comments = bsp_param_comments(design))
}

#' Write an oligo order sheet.
#' @param design A `bsp_design` object.
#' @param file Output path.
#' @param ranks Pair ranks to include; all pairs by default.
#' @param sep Field separator.
#' @return The output path, invisibly.
#' @export
write_order_sheet <- function(design, file, ranks = NULL, sep = "\t") {
  pairs <- design_pairs(design, ranks)
  sheet <- data.frame(
    name = as.vector(rbind(paste0(pair_names(design, pairs), "_F"), paste0(pair_names(design, pairs), "_R"))),
    sequence = as.vector(rbind(pairs$left_sequence, pairs$right_sequence)),
    length = as.vector(rbind(pairs$left_length, pairs$right_length)),
    tm = as.vector(rbind(pairs$left_tm, pairs$right_tm)),
    stringsAsFactors = FALSE
  )
  write_primer_results(sheet, file, sep = sep)
}

#' Amplicon reference sequences for BSP read analysis.
#'
#' For each pair returns the amplicon from the left primer to the right primer
#' on the assayed strand as genomic (unconverted), fully converted
#' (unmethylated) and CpG-methylated converted sequence, the references
#' expected by tools such as BiQ Analyzer or QUMA.
#' @param design A `bsp_design` object.
#' @param ranks Pair ranks; all pairs by default.
#' @return A named character vector.
#' @export
bsp_amplicons <- function(design, ranks = NULL) {
  pairs <- design_pairs(design, ranks)
  names <- pair_names(design, pairs)
  out <- lapply(seq_len(nrow(pairs)), function(i) {
    p <- pairs[i, ]
    template <- design$templates[[p$strand]]
    header <- sprintf("%s|%s|assay=%s|strand=%s|TSS=%s", names[i], p$amplicon_interval, p$assay_strand, p$strand, p$amplicon_tss)
    span <- function(x) substr(x, p$template_left_start, p$template_right_end)
    stats::setNames(c(span(template$original), span(template$converted), span(template$converted_methylated)),
                    paste0(header, c("|genomic", "|converted_unmethylated", "|converted_methylated")))
  })
  unlist(out)
}

#' Write amplicon reference sequences as FASTA.
#' @inheritParams bsp_amplicons
#' @param file Output FASTA path.
#' @return The output path, invisibly.
#' @export
write_amplicon_references <- function(design, file, ranks = NULL) {
  amplicons <- bsp_amplicons(design, ranks)
  writeLines(as.vector(rbind(paste0(">", names(amplicons)), unname(amplicons))), file)
  invisible(file)
}
