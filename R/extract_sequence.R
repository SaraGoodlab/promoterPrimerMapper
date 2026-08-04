reverse_complement_text <- function(sequence) {
  if (requireNamespace("Biostrings", quietly = TRUE)) return(as.character(Biostrings::reverseComplement(Biostrings::DNAString(sequence))))
  chars <- strsplit(toupper(sequence), "", fixed = TRUE)[[1L]]
  paste(rev(strsplit(chartr("ACGTN", "TGCAN", paste(chars, collapse = "")), "", fixed = TRUE)[[1L]]), collapse = "")
}

reference_sequence <- function(reference, scaffold, region_start, region_end) {
  if (inherits(reference, "FaFile")) {
    if (!Rsamtools::isOpen(reference)) base::open(reference)
    index <- Rsamtools::scanFaIndex(reference)
    available <- as.character(GenomicRanges::seqnames(index))
    if (!scaffold %in% available) stop(sprintf("FASTA scaffold '%s' does not match the GTF scaffold names.", scaffold), call. = FALSE)
    scaffold_length <- max(GenomicRanges::end(index)[available == scaffold])
    if (region_start < 1L || region_end > scaffold_length) stop(sprintf("Extraction interval %s:%s-%s exceeds scaffold boundaries (1-%s).", scaffold, region_start, region_end, scaffold_length), call. = FALSE)
    gr <- GenomicRanges::GRanges(scaffold, IRanges::IRanges(region_start, region_end), strand = "*")
    return(as.character(Rsamtools::scanFa(reference, gr)[[1L]]))
  }
  if (is.character(reference) && length(reference) == 1L && file.exists(reference)) {
    index_path <- paste0(reference, ".fai")
    if (!file.exists(index_path)) stop(sprintf("FASTA index is missing: %s", index_path), call. = FALSE)
    fa <- Rsamtools::FaFile(reference)
    on.exit(if (Rsamtools::isOpen(fa)) base::close(fa), add = TRUE)
    return(reference_sequence(fa, scaffold, region_start, region_end))
  }
  if (inherits(reference, "DNAStringSet") || (is.character(reference) && !is.null(names(reference)))) {
    if (!scaffold %in% names(reference)) stop(sprintf("FASTA scaffold '%s' does not match the GTF scaffold names.", scaffold), call. = FALSE)
    sequence <- as.character(reference[[scaffold]])
    if (region_start < 1L || region_end > nchar(sequence)) stop(sprintf("Extraction interval %s:%s-%s exceeds scaffold boundaries (1-%s).", scaffold, region_start, region_end, nchar(sequence)), call. = FALSE)
    return(substr(sequence, region_start, region_end))
  }
  stop("reference must be an indexed FASTA path, FaFile, named character sequence, or DNAStringSet.", call. = FALSE)
}

ensure_fasta_index <- function(fasta_path, index_path) {
  if (length(fasta_path) != 1L || !nzchar(fasta_path) || !file.exists(fasta_path)) {
    stop("Reference FASTA path does not exist.", call. = FALSE)
  }
  adjacent_index <- paste0(fasta_path, ".fai")
  supplied_index <- trimws(index_path)
  if (!nzchar(supplied_index) || !file.exists(supplied_index)) {
    stop("FASTA index path is required and does not exist.", call. = FALSE)
  }
  if (normalizePath(supplied_index, mustWork = TRUE) == normalizePath(adjacent_index, mustWork = FALSE)) return(fasta_path)
  Rsamtools::FaFile(fasta_path, index = supplied_index)
}

#' Extract a TSS-centered sequence in transcriptional orientation.
#' @param tss_record One row from [collapse_tss()] or a named list.
#' @param reference Indexed FASTA path, `Rsamtools::FaFile`, named sequence vector, or DNAStringSet.
#' @param upstream Number of bases before the TSS in transcriptional orientation.
#' @param downstream Number of bases after the TSS in transcriptional orientation.
#' @return A `promoter_sequence` object containing sequence and complete mapping metadata.
#' @export
extract_tss_region <- function(tss_record, reference, upstream, downstream) {
  upstream <- assert_scalar_integer(upstream, "Upstream", 0L)
  downstream <- assert_scalar_integer(downstream, "Downstream", 0L)
  scaffold <- as.character(get_field(tss_record, "scaffold"))
  strand <- validate_strand(as.character(get_field(tss_record, "strand")))
  tss <- assert_scalar_integer(get_field(tss_record, "tss"), "TSS", 1L)
  if (strand == "+") {
    region_start <- tss - upstream
    region_end <- tss + downstream
  } else {
    region_start <- tss - downstream
    region_end <- tss + upstream
  }
  if (region_start < 1L) stop("Extraction interval starts before the beginning of the reference scaffold.", call. = FALSE)
  genomic <- reference_sequence(reference, scaffold, region_start, region_end)
  sequence <- if (strand == "-") reverse_complement_text(genomic) else genomic
  metadata <- list(
    gene_id = as.character(get_field(tss_record, "gene_id", NA_character_)),
    gene_name = as.character(get_field(tss_record, "gene_name", NA_character_)),
    tss = tss,
    transcript_ids = as.character(get_field(tss_record, "transcript_ids", NA_character_)),
    transcript_id_list = get_field(tss_record, "transcript_id_list", list(NA_character_)),
    scaffold = scaffold,
    strand = strand,
    region_start = region_start,
    region_end = region_end,
    upstream = upstream,
    downstream = downstream,
    sequence_length = nchar(sequence),
    sequence = sequence
  )
  if (metadata$sequence_length != upstream + downstream + 1L) stop("Extracted sequence length does not match the requested interval.", call. = FALSE)
  class(metadata) <- c("promoter_sequence", "list")
  metadata
}

#' Build a FASTA header for an extracted sequence.
#' @param extracted_sequence A `promoter_sequence` object.
#' @return A FASTA header without the leading `>`.
#' @export
fasta_header <- function(extracted_sequence) {
  validate_extracted_sequence(extracted_sequence)
  gene <- get_field(extracted_sequence, "gene_id", NA_character_)
  if (is.na(gene) || !nzchar(gene)) gene <- get_field(extracted_sequence, "gene_name", "unknown_gene")
  sprintf("%s|%s:%s-%s|strand=%s|TSS=%s|upstream=%s|downstream=%s", gene,
          get_field(extracted_sequence, "scaffold"), get_field(extracted_sequence, "region_start"),
          get_field(extracted_sequence, "region_end"), get_field(extracted_sequence, "strand"),
          get_field(extracted_sequence, "tss"), get_field(extracted_sequence, "upstream", NA_integer_),
          get_field(extracted_sequence, "downstream", NA_integer_))
}

#' Convert a sequence to its bisulfite-treated representation.
#'
#' This models conversion of unmethylated cytosines as `C -> T`. Methylation
#' status is not available from the reference FASTA, so all cytosines are
#' converted.
#' @param sequence DNA sequence.
#' @return Uppercase bisulfite-treated sequence.
#' @export
bisulfite_convert <- function(sequence) {
  if (length(sequence) != 1L || is.na(sequence) || !nzchar(sequence)) {
    stop("sequence must be one non-empty DNA sequence.", call. = FALSE)
  }
  chartr("C", "T", toupper(as.character(sequence)))
}

#' Generate the four sequence FASTA variants.
#'
#' The complementary sequence is the reverse complement, reported 5' to 3'.
#' Bisulfite conversion is applied independently to each strand, so the four
#' outputs are requested, complementary, bisulfite-requested, and
#' bisulfite-complementary.
#' @param extracted_sequence A `promoter_sequence` object.
#' @return A named character vector containing four sequences.
#' @export
sequence_variants <- function(extracted_sequence) {
  validate_extracted_sequence(extracted_sequence)
  requested <- get_field(extracted_sequence, "sequence")
  complementary <- reverse_complement_text(requested)
  c(
    requested = requested,
    complementary = complementary,
    bisulfite_requested = bisulfite_convert(requested),
    bisulfite_complementary = bisulfite_convert(complementary)
  )
}

#' Write one sequence variant as a FASTA file.
#' @param extracted_sequence A `promoter_sequence` object.
#' @param variant One of `requested`, `complementary`, `bisulfite_requested`, or `bisulfite_complementary`.
#' @param file Output FASTA path.
#' @return The output path, invisibly.
#' @export
write_sequence_fasta <- function(extracted_sequence, variant = "requested", file) {
  variants <- sequence_variants(extracted_sequence)
  if (length(variant) != 1L || !variant %in% names(variants)) {
    stop(sprintf("variant must be one of: %s.", paste(names(variants), collapse = ", ")), call. = FALSE)
  }
  header <- paste0(fasta_header(extracted_sequence), "|variant=", variant)
  writeLines(c(paste0(">", header), unname(variants[[variant]])), file)
  invisible(file)
}

write_promoter_fasta <- function(extracted_sequence, file) {
  write_sequence_fasta(extracted_sequence, "requested", file)
}
