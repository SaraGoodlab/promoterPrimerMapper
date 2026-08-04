#' Convert an external primer start to the package's internal 1-based start.
#' @param start Numeric external start.
#' @param convention Either `one-based` or `zero-based`.
#' @return A 1-based integer.
#' @export
normalize_primer_start <- function(start, convention = c("one-based", "zero-based")) {
  convention <- match.arg(convention)
  minimum <- if (convention == "zero-based") 0L else 1L
  start <- assert_scalar_integer(start, "Primer start", minimum)
  if (convention == "zero-based") start + 1L else start
}

#' Map a 1-based interval in transcription-oriented sequence space to reference coordinates.
#' @param relative_start 1-based sequence start.
#' @param relative_end 1-based sequence end.
#' @param region_start Genomic extraction interval start.
#' @param region_end Genomic extraction interval end.
#' @param strand Gene strand.
#' @param sequence_length Optional sequence length for boundary validation.
#' @return A two-column data frame with genomic start and end.
#' @export
map_relative_interval <- function(relative_start, relative_end, region_start, region_end, strand, sequence_length = NULL) {
  relative_start <- assert_scalar_integer(relative_start, "Relative start", 1L)
  relative_end <- assert_scalar_integer(relative_end, "Relative end", relative_start)
  region_start <- assert_scalar_integer(region_start, "Region start", 1L)
  region_end <- assert_scalar_integer(region_end, "Region end", region_start)
  strand <- validate_strand(strand)
  if (!is.null(sequence_length)) {
    sequence_length <- assert_scalar_integer(sequence_length, "Sequence length", 1L)
    if (relative_end > sequence_length) stop("Primer end is outside the extracted sequence.", call. = FALSE)
  }
  if (strand == "+") {
    genomic_start <- region_start + relative_start - 1L
    genomic_end <- region_start + relative_end - 1L
  } else {
    genomic_start <- region_end - relative_end + 1L
    genomic_end <- region_end - relative_start + 1L
  }
  if (genomic_start > genomic_end) stop("Mapped genomic interval is invalid.", call. = FALSE)
  data.frame(genomic_start = as.integer(genomic_start), genomic_end = as.integer(genomic_end), stringsAsFactors = FALSE)
}

primer_binding_strand <- function(gene_strand, primer_type) {
  gene_strand <- validate_strand(gene_strand)
  is_left <- grepl("^left", tolower(trimws(primer_type)))
  if (!is_left && !grepl("^right", tolower(trimws(primer_type)))) return(NA_character_)
  if (gene_strand == "+") if (is_left) "+" else "-" else if (is_left) "-" else "+"
}

#' Parse primer start and size values from pasted primer-tool output.
#' @param text Character scalar containing pasted output.
#' @return A data frame with `primer_type`, `start`, and `size`.
#' @export
parse_primer_output <- function(text) {
  if (length(text) != 1L || is.na(text) || !nzchar(trimws(text))) stop("Primer table is empty.", call. = FALSE)
  lines <- strsplit(text, "\\r?\\n")[[1L]]
  rows <- lapply(lines, function(line) {
    clean <- trimws(line)
    if (!nzchar(clean) || grepl("^(primer|product size|left primer[[:space:]]+start)", clean, ignore.case = TRUE)) return(NULL)
    label <- sub("^([A-Za-z][A-Za-z _-]*?)[[:space:]]+[0-9].*$", "\\1", clean)
    if (identical(label, clean)) return(NULL)
    numbers <- regmatches(clean, gregexpr("[0-9]+(?:\\.[0-9]+)?", clean, perl = TRUE))[[1L]]
    if (length(numbers) < 2L) return(NULL)
    type <- trimws(label)
    if (!grepl("^(left|right)", type, ignore.case = TRUE)) return(NULL)
    data.frame(primer_type = type, start = as.numeric(numbers[1L]), size = as.numeric(numbers[2L]), stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(result) || !nrow(result)) stop("No valid Left primer or Right primer rows were found.", call. = FALSE)
  result$start <- vapply(result$start, assert_scalar_integer, integer(1), name = "Primer start", minimum = 0L)
  result$size <- vapply(result$size, assert_scalar_integer, integer(1), name = "Primer size", minimum = 1L)
  result
}

#' Map primer rows to genomic coordinates without sequence alignment.
#' @param primers Data frame with `primer_type`, `start`, and `size`, or one such row as a list.
#' @param extracted_sequence A `promoter_sequence` object.
#' @param start_convention External convention: `one-based` or `zero-based`.
#' @return A data frame containing relative and genomic intervals.
#' @export
map_primer_coordinates <- function(primers, extracted_sequence, start_convention = c("one-based", "zero-based")) {
  start_convention <- match.arg(start_convention)
  metadata <- validate_extracted_sequence(extracted_sequence)
  if (is.list(primers) && !is.data.frame(primers) && !is.null(primers$start)) primers <- as.data.frame(primers, stringsAsFactors = FALSE)
  if (!is.data.frame(primers) || !all(c("start", "size") %in% names(primers))) stop("primers must contain start and size columns.", call. = FALSE)
  if (!"primer_type" %in% names(primers)) primers$primer_type <- paste0("Primer ", seq_len(nrow(primers)))
  if (!nrow(primers)) stop("No primer rows were supplied.", call. = FALSE)
  result <- lapply(seq_len(nrow(primers)), function(i) {
    relative_start <- normalize_primer_start(primers$start[i], start_convention)
    size <- assert_scalar_integer(primers$size[i], "Primer size", 1L)
    relative_end <- relative_start + size - 1L
    if (relative_end > metadata$sequence_length) stop(sprintf("Primer '%s' ends outside the extracted sequence.", primers$primer_type[i]), call. = FALSE)
    genomic <- map_relative_interval(relative_start, relative_end, metadata$start, metadata$end, metadata$strand, metadata$sequence_length)
    data.frame(
      primer_type = as.character(primers$primer_type[i]), relative_start = relative_start, relative_end = relative_end,
      size = size, scaffold = get_field(extracted_sequence, "scaffold"), strand = metadata$strand,
      binding_strand = primer_binding_strand(metadata$strand, primers$primer_type[i]),
      genomic_start = genomic$genomic_start, genomic_end = genomic$genomic_end,
      genomic_interval = sprintf("%s:%s-%s", get_field(extracted_sequence, "scaffold"), genomic$genomic_start, genomic$genomic_end),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

write_primer_results <- function(results, file, sep = "\t") {
  if (!is.data.frame(results)) stop("results must be a data frame.", call. = FALSE)
  utils::write.table(results, file, sep = sep, row.names = FALSE, quote = FALSE, na = "")
  invisible(file)
}
