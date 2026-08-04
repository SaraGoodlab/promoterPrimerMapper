# Validation helpers

assert_scalar_integer <- function(x, name, minimum = NULL) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) || x != as.integer(x)) {
    stop(sprintf("%s must be one finite integer.", name), call. = FALSE)
  }
  value <- as.integer(x)
  if (!is.null(minimum) && value < minimum) {
    stop(sprintf("%s must be >= %s.", name, minimum), call. = FALSE)
  }
  value
}

get_field <- function(x, field, default = NULL) {
  if (is.data.frame(x) && field %in% names(x)) return(x[[field]][[1L]])
  if (is.list(x) && !is.null(x[[field]])) return(x[[field]][[1L]])
  if (!is.null(default)) return(default)
  NULL
}

validate_strand <- function(strand) {
  if (length(strand) != 1L || is.na(strand) || !strand %in% c("+", "-")) {
    stop("Strand must be '+' or '-'.", call. = FALSE)
  }
  strand
}

validate_extracted_sequence <- function(extracted_sequence) {
  required <- c("scaffold", "strand", "tss", "region_start", "region_end", "sequence")
  missing <- required[vapply(required, function(x) is.null(get_field(extracted_sequence, x)), logical(1))]
  if (length(missing)) stop(sprintf("Extracted sequence metadata is missing: %s.", paste(missing, collapse = ", ")), call. = FALSE)
  validate_strand(as.character(get_field(extracted_sequence, "strand")))
  assert_scalar_integer(get_field(extracted_sequence, "tss"), "TSS", 1L)
  start <- assert_scalar_integer(get_field(extracted_sequence, "region_start"), "Region start", 1L)
  end <- assert_scalar_integer(get_field(extracted_sequence, "region_end"), "Region end", start)
  sequence <- as.character(get_field(extracted_sequence, "sequence"))
  if (!length(sequence) || nchar(sequence) == 0L) stop("Extracted sequence is empty.", call. = FALSE)
  list(start = start, end = end, strand = as.character(get_field(extracted_sequence, "strand")), sequence_length = nchar(sequence))
}
