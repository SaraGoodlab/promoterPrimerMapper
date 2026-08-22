# Visual summaries of an extracted promoter and a BSP design.

#' Plot CpG sites, CpG islands, targets and BSP amplicons around a TSS.
#'
#' The x axis is the position in the extracted (transcription-oriented)
#' sequence, labelled below with TSS-relative positions and above with genomic
#' coordinates. The plot is a single base-graphics panel whose x data
#' coordinates are sequence positions, so Shiny brush, click and hover events
#' map directly to positions.
#' @param extracted A `promoter_sequence` object.
#' @param design Optional `bsp_design`; its targets are shaded.
#' @param rank Optional pair rank whose amplicon and primers are drawn.
#' @param targets Optional targets table (see [island_targets()]) shaded instead
#'   of the design's targets.
#' @param selection Optional sense interval `c(start, end)` drawn as a dashed box.
#' @param island Optional island index to outline.
#' @return `NULL`, invisibly; draws with base graphics.
#' @export
plot_cpg_map <- function(extracted, design = NULL, rank = NULL, targets = NULL, selection = NULL, island = NULL) {
  metadata <- validate_extracted_sequence(extracted)
  n <- metadata$sequence_length
  tss <- tss_sense_position(extracted)
  cpg <- cpg_sites(extracted)
  islands <- cpg_island_scan(extracted)
  if (is.null(targets) && !is.null(design)) targets <- design$targets
  pair <- NULL
  if (!is.null(design) && !is.null(rank) && !is.null(design$pairs) && rank %in% design$pairs$rank) {
    pair <- design$pairs[design$pairs$rank == rank, , drop = FALSE]
  }
  old <- graphics::par(mar = c(4, 1, 3.6, 1))
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot(NA, xlim = c(1, n), ylim = c(-1.45, 1.5), axes = FALSE, xlab = "", ylab = "")
  # Bottom axis: TSS-relative positions (no position 0; the TSS is +1).
  relative <- pretty(tss_relative_position(c(1L, n), tss), n = 8L)
  relative[relative == 0] <- 1L
  at <- sense_position_from_tss(relative, tss)
  keep <- at >= 1L & at <= n
  graphics::axis(1, at = at[keep], labels = format_tss_relative(relative[keep]), cex.axis = 0.8)
  graphics::mtext("Position relative to TSS (bp; TSS = +1)", side = 1, line = 2.4, cex = 0.85)
  # Top axis: genomic coordinates, decreasing to the right on minus-strand genes.
  genomic <- pretty(c(metadata$start, metadata$end), n = 6L)
  genomic <- genomic[genomic >= metadata$start & genomic <= metadata$end]
  at <- if (metadata$strand == "+") genomic - metadata$start + 1L else metadata$end - genomic + 1L
  graphics::axis(3, at = at, labels = format(genomic, big.mark = ",", scientific = FALSE, trim = TRUE), cex.axis = 0.8)
  graphics::mtext(sprintf("%s genomic position (%s strand)", get_field(extracted, "scaffold"), metadata$strand), side = 3, line = 2.2, cex = 0.85)
  if (!is.null(targets) && nrow(targets)) {
    current <- if (is.null(pair)) rep(FALSE, nrow(targets)) else targets$target_id == pair$target_id
    graphics::rect(targets$start - 0.5, -0.9, targets$end + 0.5, 1.36, col = ifelse(current, "#d5d8dc", "#f2f2f2"), border = NA)
    # Island targets are named by their island label below the bar.
    named <- targets$kind != "island"
    if (any(named)) graphics::text(targets$start[named], 1.27, targets$target_id[named], pos = 4, offset = 0.2, cex = 0.7, col = "#566573")
  }
  if (nrow(islands)) {
    outlined <- seq_len(nrow(islands)) %in% island
    graphics::rect(islands$start - 0.5, -0.07, islands$end + 0.5, 0.07, col = "#27ae60",
                   border = ifelse(outlined, "#145a32", NA), lwd = ifelse(outlined, 2.5, 1))
    graphics::text((islands$start + islands$end) / 2, -0.17, paste0("CGI", seq_len(nrow(islands))), cex = 0.7, col = "#1e8449", font = 2)
  }
  graphics::segments(cpg$position, 0, cpg$position, 1, col = "#7f8c8d")
  graphics::points(cpg$position, rep(1, nrow(cpg)), pch = 21, bg = "#c0392b", col = "#c0392b", cex = 0.7)
  graphics::abline(h = 0, col = "#7f8c8d")
  graphics::abline(v = tss, lty = 2, col = "#2c3e50")
  graphics::text(tss, 1.45, "TSS", pos = 4, offset = 0.2, cex = 0.8)
  if (!is.null(selection) && length(selection) == 2L && !anyNA(selection)) {
    graphics::rect(min(selection) - 0.5, -0.9, max(selection) + 0.5, 1.36, border = "#2980b9", lty = 2, lwd = 1.5, col = NA)
  }
  if (!is.null(pair)) {
    template <- design$templates[[pair$strand]]
    amplicon <- template_to_sense(template, pair$template_left_start, pair$template_right_end)
    left <- template_to_sense(template, pair$template_left_start, pair$template_left_end)
    right <- template_to_sense(template, pair$template_right_start, pair$template_right_end)
    graphics::rect(amplicon$start - 0.5, -0.5, amplicon$end + 0.5, -0.34, col = "#eaf2fb", border = "#2980b9")
    graphics::rect(c(left$start, right$start) - 0.5, -0.68, c(left$end, right$end) + 0.5, -0.56, col = "#f1c40f", border = NA)
    graphics::text((amplicon$start + amplicon$end) / 2, -0.8,
                   sprintf("Pair %s (%s): %s bp, %s CpGs", pair$rank, pair$target_id, pair$product_size, pair$cpg_in_amplicon), cex = 0.8)
  }
  key <- c("CpG island" = "#27ae60", "Target" = "#f2f2f2", "Selected range" = NA, "Amplicon" = "#eaf2fb", "Primer" = "#f1c40f")
  shown <- c(nrow(islands) > 0L, !is.null(targets) && nrow(targets) > 0L, !is.null(selection), !is.null(pair), !is.null(pair))
  if (any(shown)) {
    graphics::legend(1, -0.97, legend = names(key)[shown], fill = key[shown], border = c("#7f8c8d", "#7f8c8d", "#2980b9", "#7f8c8d", "#7f8c8d")[shown],
                     bty = "n", cex = 0.7, horiz = TRUE)
  }
  invisible(NULL)
}

# Sense interval from the x range of a plot brush, rounded and clamped to the sequence.
brush_to_sense <- function(xmin, xmax, n) {
  bounds <- sort(round(c(xmin, xmax)))
  as.integer(pmin(pmax(bounds, 1), n))
}

# Index of the CpG island under a click on plot_cpg_map(): the click must fall
# on the island bar or its label (y within `band`). NA when there is none.
island_at <- function(islands, x, y = NULL, band = c(-0.3, 0.3)) {
  if (is.null(x) || !nrow(islands)) return(NA_integer_)
  if (!is.null(y) && (y < band[1L] || y > band[2L])) return(NA_integer_)
  hit <- which(islands$start - 0.5 <= x & islands$end + 0.5 >= x)
  if (length(hit)) hit[1L] else NA_integer_
}

# One-row description of a sense interval: TSS-relative, sequence and genomic
# positions, length, CpG count, GC content and overlapping CpG islands.
describe_region <- function(extracted, start, end) {
  metadata <- validate_extracted_sequence(extracted)
  start <- assert_scalar_integer(start, "Start", 1L)
  end <- assert_scalar_integer(end, "End", start)
  if (end > metadata$sequence_length) stop("The region extends beyond the extracted sequence.", call. = FALSE)
  tss <- tss_sense_position(extracted)
  genomic <- sense_to_genomic(extracted, start, end)
  bases <- split_bases(substr(get_field(extracted, "sequence"), start, end))
  cpg <- cpg_sites(extracted)$position
  islands <- cpg_island_scan(extracted)
  overlap <- which(islands$start <= end & islands$end >= start)
  data.frame(
    tss_relative = format_tss_interval(tss_relative_position(start, tss), tss_relative_position(end, tss)),
    sequence = sprintf("%s..%s", start, end),
    genomic = sprintf("%s:%s-%s (%s)", get_field(extracted, "scaffold"), genomic$genomic_start, genomic$genomic_end, metadata$strand),
    length = end - start + 1L,
    cpg = sum(cpg >= start & cpg + 1L <= end),
    gc_percent = round(100 * mean(bases %in% c("C", "G")), 1),
    cpg_islands = if (length(overlap)) paste0("CGI", overlap, collapse = ", ") else "none",
    stringsAsFactors = FALSE
  )
}

# Text for the position under the cursor.
describe_position <- function(extracted, x) {
  metadata <- validate_extracted_sequence(extracted)
  position <- as.integer(round(x))
  if (is.na(position) || position < 1L || position > metadata$sequence_length) return("")
  genomic <- sense_to_genomic(extracted, position, position)$genomic_start
  sprintf("Cursor: %s (sequence position %s, %s:%s, base %s)",
          format_tss_relative(tss_relative_position(position, tss_sense_position(extracted))), position,
          get_field(extracted, "scaffold"), format(genomic, big.mark = ",", scientific = FALSE),
          substr(get_field(extracted, "sequence"), position, position))
}

# Island table for display and the island select labels.
island_table <- function(extracted, params = bsp_primer_params()) {
  islands <- cpg_island_scan(extracted)
  if (!nrow(islands)) return(islands)
  genomic <- sense_to_genomic(extracted, islands$start, islands$end)
  windows <- island_targets(extracted, params)
  cpg <- cpg_sites(extracted)$position
  data.frame(
    island = paste0("CGI", seq_len(nrow(islands))),
    tss_relative = format_tss_interval(islands$tss_start, islands$tss_end),
    genomic = sprintf("%s:%s-%s", get_field(extracted, "scaffold"), genomic$genomic_start, genomic$genomic_end),
    length = islands$length,
    cpg = vapply(seq_len(nrow(islands)), function(i) sum(cpg >= islands$start[i] & cpg + 1L <= islands$end[i]), integer(1)),
    gc_percent = islands$gc_percent,
    cpg_obs_exp = islands$cpg_obs_exp,
    design_window = format_tss_interval(windows$tss_start, windows$tss_end),
    stringsAsFactors = FALSE
  )
}

island_choices <- function(table) {
  if (!nrow(table)) return(character(0))
  stats::setNames(as.character(seq_len(nrow(table))),
                  sprintf("%s: %s (%s bp, %s CpGs, GC %s%%)", table$island, table$tss_relative, table$length, table$cpg, table$gc_percent))
}

# HTML view of the converted template of one pair: primers highlighted, CpGs in
# bold, amplicon shaded, TSS outlined. Each line starts with the TSS-relative
# position of its first base.
sequence_view_html <- function(design, rank, line_width = 60L) {
  pair <- design_pairs(design, rank)
  template <- design$templates[[pair$strand]]
  bases <- split_bases(template$converted_methylated)
  n <- length(bases)
  classes <- rep("", n)
  add <- function(index, class) classes[index] <<- trimws(paste(classes[index], class))
  add(pair$template_left_start:pair$template_right_end, "amp")
  add(c(pair$template_left_start:pair$template_left_end, pair$template_right_start:pair$template_right_end), "primer")
  cpg <- which(template$cpg_c)
  add(c(cpg, cpg + 1L), "cpg")
  tss <- tss_sense_position(design$extracted)
  tss_template <- if (template$strand == "sense") tss else n - tss + 1L
  add(tss_template, "tss")
  starts <- seq(1L, n, by = line_width)
  sense_start <- if (template$strand == "sense") starts else n - starts + 1L
  labels <- format(format_tss_relative(tss_relative_position(sense_start, tss)), width = 6L, justify = "right")
  lines <- vapply(seq_along(starts), function(i) {
    index <- starts[i]:min(n, starts[i] + line_width - 1L)
    runs <- rle(classes[index])
    ends <- cumsum(runs$lengths)
    segments <- vapply(seq_along(runs$lengths), function(j) {
      text <- paste(bases[index[(ends[j] - runs$lengths[j] + 1L):ends[j]]], collapse = "")
      if (nzchar(runs$values[j])) sprintf("<span class=\"%s\">%s</span>", runs$values[j], text) else text
    }, character(1))
    paste0(labels[i], "  ", paste(segments, collapse = ""))
  }, character(1))
  legend <- paste0("<div class=\"bsp-legend small\"><span class=\"primer\">primer</span><span class=\"amp\">amplicon</span>",
                   "<span class=\"cpg\">CpG</span><span class=\"tss\">TSS</span> ",
                   sprintf("Converted %s template (%s), 5'&rarr;3'; CpG cytosines shown unconverted.", pair$strand, pair$assay_strand), "</div>")
  shiny::HTML(paste0(legend, "<pre class=\"bsp-seq\">", paste(lines, collapse = "\n"), "</pre>"))
}
