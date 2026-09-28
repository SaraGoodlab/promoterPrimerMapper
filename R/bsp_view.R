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
#' @param selection Optional sense interval `c(start, end)` drawn as a solid box.
#' @param island Optional island index to outline.
#' @param view Optional two sense positions to zoom to, with 50 bp of padding.
#' @return `NULL`, invisibly; draws with base graphics.
#' @export
plot_cpg_map <- function(extracted, design = NULL, rank = NULL, targets = NULL, selection = NULL, island = NULL, view = NULL) {
  metadata <- validate_extracted_sequence(extracted)
  n <- metadata$sequence_length
  bounds <- c(1L, n)
  if (!is.null(view)) {
    if (!is.numeric(view) || length(view) != 2L || any(!is.finite(view))) stop("view must contain two finite sense positions.", call. = FALSE)
    bounds <- c(max(1, min(view) - 50), min(n, max(view) + 50))
    if (bounds[1L] >= bounds[2L]) stop("view must overlap the extracted sequence.", call. = FALSE)
  }
  tss <- tss_sense_position(extracted)
  cpg <- cpg_sites(extracted)
  islands <- cpg_island_scan(extracted)
  if (is.null(targets) && !is.null(design)) targets <- design$targets
  pair <- NULL
  if (!is.null(design) && !is.null(rank) && !is.null(design$pairs) && rank %in% design$pairs$rank) {
    pair <- design$pairs[design$pairs$rank == rank, , drop = FALSE]
  }
  colors <- app_colors()
  chart_cex <- 11 / 12
  device_width <- grDevices::dev.size("in")[1L]
  target_fill <- grDevices::adjustcolor(colors[["primary"]], alpha.f = .09)
  cpg_color <- grDevices::adjustcolor(colors[["cpg"]], alpha.f = .70)
  old <- graphics::par(mar = c(4, 1, 3.6, 1), bg = colors[["surface"]], fg = colors[["ink"]],
                       col = colors[["ink"]], col.axis = colors[["muted"]], col.lab = colors[["muted"]],
                       ps = 12, cex = 1, family = "sans")
  on.exit(graphics::par(old), add = TRUE)
  # Narrow devices stack the legend into three rows; leave room above the axis.
  graphics::plot(NA, xlim = bounds, ylim = c(if (device_width < 6) -2.1 else -1.65, 1.55), axes = FALSE, xlab = "", ylab = "")
  # Bottom axis: TSS-relative positions (no position 0; the TSS is +1).
  relative <- pretty(tss_relative_position(bounds, tss), n = max(3L, min(8L, floor(device_width))))
  relative[relative == 0] <- 1L
  at <- sense_position_from_tss(relative, tss)
  keep <- at >= bounds[1L] & at <= bounds[2L]
  graphics::axis(1, at = at[keep], labels = format_tss_relative(relative[keep]), cex.axis = chart_cex, family = "mono", col = colors[["line"]])
  graphics::mtext("Position relative to TSS (bp; TSS = +1)", side = 1, line = 2.4, cex = 1)
  # Top axis: genomic coordinates, decreasing to the right on minus-strand genes.
  limits <- sense_to_genomic(extracted, as.integer(bounds[1L]), as.integer(bounds[2L]))
  genomic <- pretty(c(limits$genomic_start, limits$genomic_end), n = max(2L, min(6L, floor(device_width / 1.5))))
  genomic <- genomic[genomic >= metadata$start & genomic <= metadata$end]
  at <- if (metadata$strand == "+") genomic - metadata$start + 1L else metadata$end - genomic + 1L
  graphics::axis(3, at = at, labels = format(genomic, big.mark = ",", scientific = FALSE, trim = TRUE), cex.axis = chart_cex, family = "mono", col = colors[["line"]])
  graphics::mtext(sprintf("%s genomic position (%s strand)", get_field(extracted, "scaffold"), metadata$strand), side = 3, line = 2.2, cex = 1)
  if (!is.null(targets) && nrow(targets)) {
    current <- if (is.null(pair)) rep(FALSE, nrow(targets)) else targets$target_id == pair$target_id
    graphics::rect(targets$start - 0.5, -0.9, targets$end + 0.5, 1.36, col = target_fill,
                   border = colors[["primary"]], lty = 2, lwd = ifelse(current, 1.25, .8))
    # Island targets are named by their island label below the bar.
    named <- targets$kind != "island"
    if (any(named)) graphics::text(targets$start[named], 1.27, targets$target_id[named], pos = 4, offset = 0.2, cex = chart_cex, col = colors[["muted"]])
  }
  if (nrow(islands)) {
    outlined <- seq_len(nrow(islands)) %in% island
    if (any(outlined)) graphics::rect(islands$start[outlined] - .5, -.30, islands$end[outlined] + .5, .12,
                                     col = colors[["primary-soft"]], border = NA)
    graphics::rect(islands$start - 0.5, -0.07, islands$end + 0.5, 0.07, col = colors[["island"]],
                   border = ifelse(outlined, colors[["primary"]], NA), lwd = ifelse(outlined, 2.5, 1))
    graphics::text((islands$start + islands$end) / 2, -0.19, paste0("CGI", seq_len(nrow(islands))), cex = chart_cex,
                   col = ifelse(outlined, colors[["primary"]], colors[["ink"]]))
  }
  graphics::segments(cpg$position, 0, cpg$position, 1, col = cpg_color, lwd = .6)
  graphics::points(cpg$position, rep(1, nrow(cpg)), pch = 16, col = cpg_color, cex = .55)
  graphics::abline(h = 0, col = colors[["line"]])
  graphics::abline(v = tss, lty = 2, col = colors[["ink"]])
  if (tss >= bounds[1L] && tss <= bounds[2L]) graphics::text(tss, 1.45, "TSS", pos = 4, offset = 0.2, cex = chart_cex)
  if (!is.null(selection) && length(selection) == 2L && !anyNA(selection)) {
    graphics::rect(min(selection) - 0.5, -0.9, max(selection) + 0.5, 1.36, border = colors[["primary"]], lty = 1, lwd = 1.5, col = NA)
  }
  if (!is.null(pair)) {
    template <- design$templates[[pair$strand]]
    amplicon <- template_to_sense(template, pair$template_left_start, pair$template_right_end)
    left <- template_to_sense(template, pair$template_left_start, pair$template_left_end)
    right <- template_to_sense(template, pair$template_right_start, pair$template_right_end)
    graphics::rect(amplicon$start - 0.5, -0.5, amplicon$end + 0.5, -0.34, col = NA, border = colors[["primary"]], lwd = 1.5)
    graphics::rect(c(left$start, right$start) - 0.5, -0.68, c(left$end, right$end) + 0.5, -0.56, col = colors[["primer"]], border = NA)
    label <- sprintf("Pair %s (%s): %s bp, %s CpGs", pair$rank, pair$target_id, pair$product_size, pair$cpg_in_amplicon)
    half_width <- graphics::strwidth(label, cex = chart_cex) / 2
    label_x <- if (2 * half_width >= diff(bounds)) mean(bounds) else
      max(bounds[1L] + half_width, min((amplicon$start + amplicon$end) / 2, bounds[2L] - half_width))
    graphics::text(label_x, -0.8, label, cex = chart_cex)
  }
  key <- c("CpG site" = colors[["cpg"]], "CpG island" = colors[["island"]], "Target" = colors[["primary"]],
           "Selected range" = colors[["primary"]], "Amplicon" = colors[["primary"]], "Primer" = colors[["primer"]])
  shown <- c(nrow(cpg) > 0L, nrow(islands) > 0L, !is.null(targets) && nrow(targets) > 0L, !is.null(selection), !is.null(pair), !is.null(pair))
  if (any(shown)) {
    graphics::legend(bounds[1L], -1.02, legend = names(key)[shown], col = key[shown],
                     pch = c(16, 15, NA, NA, NA, 15)[shown], lty = c(NA, NA, 2, 1, 1, NA)[shown],
                     lwd = c(NA, NA, 1, 1.5, 1.5, NA)[shown],
                     bty = "n", cex = chart_cex, y.intersp = 1.2,
                     ncol = if (device_width < 6) 2L else if (device_width < 9) 3L else sum(shown))
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
sequence_view_html <- function(design, rank, line_width = 60L, full = TRUE, flank = 20L) {
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
  first <- if (full) 1L else max(1L, pair$template_left_start - flank)
  last <- if (full) n else min(n, pair$template_right_end + flank)
  starts <- seq(first, last, by = line_width)
  sense_start <- if (template$strand == "sense") starts else n - starts + 1L
  labels <- format(format_tss_relative(tss_relative_position(sense_start, tss)), width = 6L, justify = "right")
  lines <- vapply(seq_along(starts), function(i) {
    index <- starts[i]:min(last, starts[i] + line_width - 1L)
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
