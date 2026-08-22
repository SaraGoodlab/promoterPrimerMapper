# Bisulfite sequence model. Positions are 1-based and, unless stated otherwise,
# in the transcription-oriented ("sense") coordinates of a promoter_sequence.

split_bases <- function(sequence) strsplit(toupper(sequence), "", fixed = TRUE)[[1L]]

# Sense position of the TSS within an extracted sequence.
tss_sense_position <- function(extracted) {
  metadata <- validate_extracted_sequence(extracted)
  tss <- assert_scalar_integer(get_field(extracted, "tss"), "TSS", 1L)
  if (metadata$strand == "+") tss - metadata$start + 1L else metadata$end - tss + 1L
}

#' Convert between sequence positions and TSS-relative positions.
#'
#' TSS-relative positions follow the promoter convention: the TSS is `+1`, the
#' base immediately upstream is `-1`, and there is no position `0`.
#' @param position 1-based sense positions in the extracted sequence.
#' @param relative TSS-relative positions (non-zero).
#' @param tss_position Sense position of the TSS.
#' @return Integer vector.
#' @export
tss_relative_position <- function(position, tss_position) {
  offset <- as.integer(position) - as.integer(tss_position)
  ifelse(offset >= 0L, offset + 1L, offset)
}

#' @rdname tss_relative_position
#' @export
sense_position_from_tss <- function(relative, tss_position) {
  relative <- as.integer(relative)
  if (any(relative == 0L, na.rm = TRUE)) stop("TSS-relative positions start at +1 (TSS) and -1; 0 is not a valid position.", call. = FALSE)
  ifelse(relative > 0L, tss_position + relative - 1L, tss_position + relative)
}

format_tss_relative <- function(x) ifelse(x > 0L, paste0("+", x), as.character(x))

# Sense interval -> genomic interval; vectorised wrapper around map_relative_interval().
sense_to_genomic <- function(extracted, start, end) {
  metadata <- validate_extracted_sequence(extracted)
  out <- lapply(seq_along(start), function(i) map_relative_interval(start[i], end[i], metadata$start, metadata$end, metadata$strand, metadata$sequence_length))
  if (!length(out)) return(data.frame(genomic_start = integer(), genomic_end = integer()))
  do.call(rbind, out)
}

#' Locate CpG dinucleotides in an extracted sequence.
#' @param extracted A `promoter_sequence` object from [extract_tss_region()].
#' @return A data frame with one row per CpG: sense `position` of the C,
#'   `tss_relative` position of the C, and the genomic interval of the
#'   dinucleotide.
#' @export
cpg_sites <- function(extracted) {
  validate_extracted_sequence(extracted)
  bases <- split_bases(get_field(extracted, "sequence"))
  n <- length(bases)
  position <- if (n > 1L) which(bases[-n] == "C" & bases[-1L] == "G") else integer()
  genomic <- sense_to_genomic(extracted, position, position + 1L)
  data.frame(
    position = position,
    tss_relative = tss_relative_position(position, tss_sense_position(extracted)),
    scaffold = rep(get_field(extracted, "scaffold"), length(position)),
    genomic_start = genomic$genomic_start,
    genomic_end = genomic$genomic_end,
    stringsAsFactors = FALSE
  )
}

#' Summarise CpG content of an extracted sequence.
#' @param extracted A `promoter_sequence` object.
#' @return A one-row data frame with length, GC percentage, CpG count, and
#'   observed/expected CpG ratio.
#' @export
cpg_summary <- function(extracted) {
  validate_extracted_sequence(extracted)
  bases <- split_bases(get_field(extracted, "sequence"))
  n <- length(bases)
  c_count <- sum(bases == "C")
  g_count <- sum(bases == "G")
  cpg <- if (n > 1L) sum(bases[-n] == "C" & bases[-1L] == "G") else 0L
  data.frame(
    length = n,
    gc_percent = round(100 * (c_count + g_count) / n, 1),
    cpg_count = cpg,
    cpg_obs_exp = if (c_count * g_count > 0) round(cpg * n / (c_count * g_count), 2) else 0
  )
}

#' Predict CpG islands in an extracted sequence.
#'
#' Sliding-window scan with Gardiner-Garden & Frommer style criteria: every
#' window of `window` bases with GC fraction >= `min_gc` and observed/expected
#' CpG ratio >= `min_oe` is marked, overlapping windows are merged, and merged
#' segments shorter than `min_length` are dropped.
#' @param extracted A `promoter_sequence` object.
#' @param window Window width in bases.
#' @param min_gc Minimum GC fraction.
#' @param min_oe Minimum observed/expected CpG ratio.
#' @param min_length Minimum island length.
#' @return A data frame of islands in sense and TSS-relative coordinates with
#'   GC percentage and observed/expected ratio of each island.
#' @export
cpg_island_scan <- function(extracted, window = 100L, min_gc = 0.5, min_oe = 0.6, min_length = 200L) {
  validate_extracted_sequence(extracted)
  window <- assert_scalar_integer(window, "Window", 2L)
  min_length <- assert_scalar_integer(min_length, "Minimum island length", 1L)
  bases <- split_bases(get_field(extracted, "sequence"))
  n <- length(bases)
  empty <- data.frame(start = integer(), end = integer(), length = integer(), tss_start = integer(), tss_end = integer(),
                      gc_percent = numeric(), cpg_obs_exp = numeric())
  if (n < window) return(empty)
  window_sum <- function(x) {
    cs <- c(0L, cumsum(x))
    cs[(window + 1L):(n + 1L)] - cs[1L:(n - window + 1L)]
  }
  is_c <- bases == "C"
  is_g <- bases == "G"
  is_cpg <- c(is_c[-n] & is_g[-1L], FALSE)
  c_count <- window_sum(is_c)
  g_count <- window_sum(is_g)
  # Count a CpG in a window only when both of its bases fall inside it.
  cpg_count <- window_sum(is_cpg) - is_cpg[window:n]
  gc <- (c_count + g_count) / window
  oe <- ifelse(c_count * g_count > 0, cpg_count * window / (c_count * g_count), 0)
  passing <- which(gc >= min_gc & oe >= min_oe)
  if (!length(passing)) return(empty)
  covered <- logical(n)
  for (s in passing) covered[s:(s + window - 1L)] <- TRUE
  runs <- rle(covered)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  keep <- runs$values & runs$lengths >= min_length
  starts <- starts[keep]
  ends <- ends[keep]
  if (!length(starts)) return(empty)
  tss <- tss_sense_position(extracted)
  stats <- t(vapply(seq_along(starts), function(i) {
    b <- bases[starts[i]:ends[i]]
    m <- length(b)
    cc <- sum(b == "C")
    gg <- sum(b == "G")
    cg <- sum(b[-m] == "C" & b[-1L] == "G")
    c(100 * (cc + gg) / m, if (cc * gg > 0) cg * m / (cc * gg) else 0)
  }, numeric(2)))
  data.frame(
    start = starts, end = ends, length = ends - starts + 1L,
    tss_start = tss_relative_position(starts, tss), tss_end = tss_relative_position(ends, tss),
    gc_percent = round(stats[, 1L], 1), cpg_obs_exp = round(stats[, 2L], 2)
  )
}

#' Build a bisulfite PCR template for one strand of an extracted sequence.
#'
#' Bisulfite treatment makes the two genomic strands non-complementary, so a
#' BSP assay targets one converted strand. `sense` uses the extracted
#' (transcription-oriented) sequence as template, `antisense` its reverse
#' complement. The genomic strand of the template is reported as the Bismark
#' strand label: `OT` (original top, +) or `OB` (original bottom, -).
#' @param extracted A `promoter_sequence` object.
#' @param strand `sense` or `antisense`.
#' @return A `bisulfite_template` list holding the template-oriented original
#'   sequence, the unmethylated and CpG-methylated conversions, per-base
#'   `converted_c` (non-CpG cytosine, converted regardless of methylation),
#'   `cpg_c` (CpG cytosine) and `cpg_mask` (methylation-dependent bases: CpG
#'   cytosines plus a 3'-edge C whose CpG context lies outside the window),
#'   and the strand labels. Only cytosines of the template strand are
#'   converted, so the G of a CpG is not methylation-dependent.
#' @export
bisulfite_template <- function(extracted, strand = c("sense", "antisense")) {
  strand <- match.arg(strand)
  metadata <- validate_extracted_sequence(extracted)
  sense <- toupper(get_field(extracted, "sequence"))
  original <- if (strand == "sense") sense else reverse_complement_text(sense)
  bases <- split_bases(original)
  n <- length(bases)
  is_c <- bases == "C"
  cpg_c <- is_c & c(bases[-1L], "") == "G"
  # A C at the 3' edge may belong to a CpG whose G lies outside the window.
  edge_unknown <- seq_len(n) == n & is_c
  cpg_mask <- cpg_c | edge_unknown
  genomic_strand <- if (strand == "sense") metadata$strand else opposite_strand(metadata$strand)
  structure(list(
    strand = strand,
    genomic_strand = genomic_strand,
    assay_strand = assay_strand_label(genomic_strand),
    original = original,
    converted = bisulfite_convert(original),
    converted_methylated = bisulfite_convert(original, cpg = "methylated"),
    converted_c = is_c & !cpg_c & !edge_unknown,
    cpg_c = cpg_c,
    cpg_mask = cpg_mask,
    length = n
  ), class = c("bisulfite_template", "list"))
}

# Template interval -> sense interval.
template_to_sense <- function(template, start, end) {
  if (template$strand == "sense") return(list(start = start, end = end))
  list(start = template$length - end + 1L, end = template$length - start + 1L)
}

# Sense interval -> template interval (the mapping is its own inverse).
sense_to_template <- function(template, start, end) template_to_sense(template, start, end)
