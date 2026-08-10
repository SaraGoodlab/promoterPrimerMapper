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

is_right_primer <- function(primer_type) {
  grepl("^(right|reverse)", tolower(trimws(primer_type)))
}

# Interval covered by a primer given its reported start. Primer3 and MethPrimer
# report the 5' end of a right primer, which is its highest coordinate.
primer_interval <- function(start, size, is_right, right_start = c("five-prime", "leftmost")) {
  right_start <- match.arg(right_start)
  if (is_right && right_start == "five-prime") c(start - size + 1L, start) else c(start, start + size - 1L)
}

# Flip a 1-based interval between a sequence and its reverse complement.
flip_interval <- function(start, end, sequence_length) {
  c(sequence_length - end + 1L, sequence_length - start + 1L)
}

opposite_strand <- function(strand) if (validate_strand(strand) == "+") "-" else "+"

assay_strand_label <- function(genomic_strand) if (validate_strand(genomic_strand) == "+") "OT" else "OB"

primer_binding_strand <- function(template_strand, primer_type) {
  template_strand <- validate_strand(template_strand)
  is_left <- grepl("^(left|forward)", tolower(trimws(primer_type)))
  if (!is_left && !is_right_primer(primer_type)) return(NA_character_)
  if (is_left) template_strand else opposite_strand(template_strand)
}

#' Parse primer start and size values from pasted primer-tool output.
#'
#' Rows must start with a `Left`/`Right` (or `Forward`/`Reverse`) label followed
#' by Start and Size. A primer sequence column (the last token, when its length
#' equals Size) and a `Product size: N` line are captured when present so the
#' mapping can be verified.
#' @param text Character scalar containing pasted output.
#' @return A data frame with `primer_type`, `start`, `size`, and `sequence`; the
#'   reported product size, if any, is stored in `attr(, "product_size")`.
#' @export
parse_primer_output <- function(text) {
  if (length(text) != 1L || is.na(text) || !nzchar(trimws(text))) stop("Primer table is empty.", call. = FALSE)
  lines <- strsplit(text, "\\r?\\n")[[1L]]
  product_line <- grep("product size", lines, ignore.case = TRUE, value = TRUE)
  product_size <- if (length(product_line)) as.integer(sub("(?i)^.*product size[^0-9]*([0-9]+).*$", "\\1", product_line[1L], perl = TRUE)) else NULL
  rows <- lapply(lines, function(line) {
    clean <- trimws(line)
    if (!nzchar(clean) || grepl("^(primer|product size|left primer[[:space:]]+start)", clean, ignore.case = TRUE)) return(NULL)
    label <- sub("^([A-Za-z][A-Za-z _-]*?)[[:space:]]+[0-9].*$", "\\1", clean)
    if (identical(label, clean)) return(NULL)
    numbers <- regmatches(clean, gregexpr("[0-9]+(?:\\.[0-9]+)?", clean, perl = TRUE))[[1L]]
    if (length(numbers) < 2L) return(NULL)
    type <- trimws(label)
    if (!grepl("^(left|right|forward|reverse)", type, ignore.case = TRUE)) return(NULL)
    tokens <- strsplit(clean, "[[:space:]]+")[[1L]]
    last <- toupper(tokens[length(tokens)])
    sequence <- if (grepl("^[ACGTUNRYKMSWBDHV]+$", last) && nchar(last) == as.numeric(numbers[2L])) last else NA_character_
    data.frame(primer_type = type, start = as.numeric(numbers[1L]), size = as.numeric(numbers[2L]), sequence = sequence, stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(result) || !nrow(result)) stop("No valid Left primer or Right primer rows were found.", call. = FALSE)
  result$start <- vapply(result$start, assert_scalar_integer, integer(1), name = "Primer start", minimum = 0L)
  result$size <- vapply(result$size, assert_scalar_integer, integer(1), name = "Primer size", minimum = 1L)
  if (!is.null(product_size) && !is.na(product_size)) attr(result, "product_size") <- product_size
  result
}

# IUPAC-aware comparison of a primer with a target sequence of the same length.
iupac_match <- function(primer, target) {
  if (nchar(primer) != nchar(target)) return(FALSE)
  codes <- c(A = "A", C = "C", G = "G", T = "T", U = "T", R = "AG", Y = "CT", K = "GT", M = "AC",
             S = "CG", W = "AT", B = "CGT", D = "AGT", H = "ACT", V = "ACG", N = "ACGT")
  p <- strsplit(toupper(primer), "", fixed = TRUE)[[1L]]
  t <- strsplit(toupper(target), "", fixed = TRUE)[[1L]]
  allowed <- codes[p]
  if (anyNA(allowed)) return(FALSE)
  all(mapply(function(a, b) grepl(b, a, fixed = TRUE), allowed, t))
}

# Which representation of the template strand a pasted primer sequence matches.
classify_primer_sequence <- function(primer, base_sequence, start, end, is_right) {
  candidates <- c(
    "converted" = bisulfite_convert(base_sequence),
    "converted (methylated CpG)" = bisulfite_convert(base_sequence, cpg = "methylated"),
    "unconverted" = toupper(base_sequence)
  )
  for (label in names(candidates)) {
    footprint <- substr(candidates[[label]], start, end)
    expected <- if (is_right) reverse_complement_text(footprint) else footprint
    if (iupac_match(primer, expected)) return(label)
  }
  "mismatch"
}

#' Map primer rows to genomic coordinates without sequence alignment.
#' @param primers Data frame with `primer_type`, `start`, and `size` (and
#'   optionally `sequence`), or one such row as a list.
#' @param extracted_sequence A `promoter_sequence` object.
#' @param start_convention External convention: `one-based` or `zero-based`.
#' @param right_start How the tool reports right/reverse primer starts:
#'   `five-prime` (Primer3, MethPrimer: the primer's 5' end, i.e. its highest
#'   position) or `leftmost`.
#' @param sequence_variant The [sequence_variants()] entry the primers were
#'   designed on. Positions on the complementary variants are flipped back to
#'   transcription orientation before mapping.
#' @param product_size Product size reported by the primer tool; when supplied
#'   the mapped amplicon length must agree. Defaults to the value captured by
#'   [parse_primer_output()].
#' @return A data frame containing relative and genomic intervals. `relative_*`
#'   positions are in the coordinates of `sequence_variant`. `binding_strand` is
#'   the genomic strand whose 5'->3' sequence the primer reproduces (the primer
#'   anneals to the opposite strand); `assay_strand` is the original genomic
#'   strand used as template (`OT` = original top/+, `OB` = original bottom/-).
#'   `sequence_match` reports whether a pasted primer sequence matches the
#'   template at the mapped position.
#' @export
map_primer_coordinates <- function(primers, extracted_sequence, start_convention = c("one-based", "zero-based"),
                                   right_start = c("five-prime", "leftmost"),
                                   sequence_variant = c("requested", "complementary", "bisulfite_requested", "bisulfite_complementary"),
                                   product_size = attr(primers, "product_size")) {
  start_convention <- match.arg(start_convention)
  right_start <- match.arg(right_start)
  sequence_variant <- match.arg(sequence_variant)
  metadata <- validate_extracted_sequence(extracted_sequence)
  if (is.list(primers) && !is.data.frame(primers) && !is.null(primers$start)) primers <- as.data.frame(primers, stringsAsFactors = FALSE)
  if (!is.data.frame(primers) || !all(c("start", "size") %in% names(primers))) stop("primers must contain start and size columns.", call. = FALSE)
  if (!"primer_type" %in% names(primers)) primers$primer_type <- paste0("Primer ", seq_len(nrow(primers)))
  if (!nrow(primers)) stop("No primer rows were supplied.", call. = FALSE)
  complementary <- sequence_variant %in% c("complementary", "bisulfite_complementary")
  sense <- toupper(as.character(get_field(extracted_sequence, "sequence")))
  base_sequence <- if (complementary) reverse_complement_text(sense) else sense
  template_strand <- if (complementary) opposite_strand(metadata$strand) else metadata$strand
  scaffold <- get_field(extracted_sequence, "scaffold")
  result <- lapply(seq_len(nrow(primers)), function(i) {
    type <- as.character(primers$primer_type[i])
    size <- assert_scalar_integer(primers$size[i], "Primer size", 1L)
    interval <- primer_interval(normalize_primer_start(primers$start[i], start_convention), size, is_right_primer(type), right_start)
    if (interval[1L] < 1L || interval[2L] > metadata$sequence_length) stop(sprintf("Primer '%s' lies outside the extracted sequence.", type), call. = FALSE)
    sense_interval <- if (complementary) flip_interval(interval[1L], interval[2L], metadata$sequence_length) else interval
    genomic <- map_relative_interval(sense_interval[1L], sense_interval[2L], metadata$start, metadata$end, metadata$strand, metadata$sequence_length)
    pasted <- if ("sequence" %in% names(primers)) primers$sequence[i] else NA_character_
    data.frame(
      primer_type = type, relative_start = interval[1L], relative_end = interval[2L],
      size = size, scaffold = scaffold, strand = metadata$strand,
      binding_strand = primer_binding_strand(template_strand, type),
      assay_strand = assay_strand_label(template_strand),
      genomic_start = genomic$genomic_start, genomic_end = genomic$genomic_end,
      genomic_interval = sprintf("%s:%s-%s", scaffold, genomic$genomic_start, genomic$genomic_end),
      sequence_match = if (is.na(pasted) || !nzchar(pasted)) NA_character_ else classify_primer_sequence(pasted, base_sequence, interval[1L], interval[2L], is_right_primer(type)),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, result)
  if (!is.null(product_size) && !is.na(product_size)) {
    product_size <- assert_scalar_integer(product_size, "Product size", 1L)
    mapped <- max(result$relative_end) - min(result$relative_start) + 1L
    if (mapped != product_size) {
      stop(sprintf("Mapped amplicon is %s bp but the primer tool reported %s bp. Check the Start convention, right-primer convention, and sequence variant.",
                   mapped, product_size), call. = FALSE)
    }
  }
  result
}

write_primer_results <- function(results, file, sep = "\t", comments = NULL) {
  if (!is.data.frame(results)) stop("results must be a data frame.", call. = FALSE)
  if (length(comments)) writeLines(paste0("# ", comments), file)
  # write.table warns when appending column names after the comment lines.
  suppressWarnings(utils::write.table(results, file, sep = sep, row.names = FALSE, quote = FALSE, na = "", append = length(comments) > 0L))
  invisible(file)
}
