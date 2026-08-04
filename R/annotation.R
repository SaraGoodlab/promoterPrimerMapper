annotation_columns <- function(annotation) {
  if (!inherits(annotation, "GRanges")) stop("annotation must be a GenomicRanges::GRanges object.", call. = FALSE)
  data.frame(
    feature = if ("type" %in% names(GenomicRanges::mcols(annotation))) as.character(GenomicRanges::mcols(annotation)$type) else rep(NA_character_, length(annotation)),
    gene_id = if ("gene_id" %in% names(GenomicRanges::mcols(annotation))) as.character(GenomicRanges::mcols(annotation)$gene_id) else rep(NA_character_, length(annotation)),
    gene_name = if ("gene_name" %in% names(GenomicRanges::mcols(annotation))) as.character(GenomicRanges::mcols(annotation)$gene_name) else rep(NA_character_, length(annotation)),
    transcript_id = if ("transcript_id" %in% names(GenomicRanges::mcols(annotation))) as.character(GenomicRanges::mcols(annotation)$transcript_id) else rep(NA_character_, length(annotation)),
    scaffold = as.character(GenomicRanges::seqnames(annotation)),
    transcript_start = GenomicRanges::start(annotation),
    transcript_end = GenomicRanges::end(annotation),
    strand = as.character(GenomicRanges::strand(annotation)),
    stringsAsFactors = FALSE
  )
}

#' Import a GTF annotation.
#'
#' @param path Path to a GTF/GFF file. Compressed input supported by rtracklayer.
#' @return A GRanges object imported by rtracklayer.
#' @export
load_annotation <- function(path) {
  if (length(path) != 1L || !nzchar(path) || !file.exists(path)) {
    stop("GTF annotation file does not exist.", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) stop("Package 'rtracklayer' is required.", call. = FALSE)
  annotation <- rtracklayer::import(path)
  if (!inherits(annotation, "GRanges")) stop("The imported annotation is not a GRanges object.", call. = FALSE)
  annotation
}

#' Find all annotation features belonging to one gene.
#'
#' @param annotation GRanges returned by [load_annotation()].
#' @param query Exact gene ID or gene name.
#' @param by One of `auto`, `gene_id`, or `gene_name`.
#' @return A GRanges subset containing all features for the matching gene.
#' @export
find_gene <- function(annotation, query, by = c("auto", "gene_id", "gene_name")) {
  by <- match.arg(by)
  fields <- annotation_columns(annotation)
  if (length(query) != 1L || is.na(query) || !nzchar(query)) stop("Gene query must be a non-empty string.", call. = FALSE)
  if (by == "auto") {
    id_hits <- !is.na(fields$gene_id) & fields$gene_id == query
    name_hits <- !is.na(fields$gene_name) & fields$gene_name == query
    hits <- id_hits | name_hits
  } else {
    hits <- !is.na(fields[[by]]) & fields[[by]] == query
  }
  if (!any(hits)) stop(sprintf("Gene '%s' was not found.", query), call. = FALSE)
  candidate_ids <- fields$gene_id[hits & !is.na(fields$gene_id)]
  matched_ids <- unique(candidate_ids[nzchar(candidate_ids)])
  if (length(matched_ids) > 1L) stop(sprintf("Gene query '%s' matches multiple gene IDs: %s.", query, paste(matched_ids, collapse = ", ")), call. = FALSE)
  annotation[hits]
}

#' Calculate transcript-boundary TSS records for a gene.
#' @param annotation GRanges annotation.
#' @param gene Optional GRanges subset or gene query string.
#' @param gene_id Optional exact gene ID.
#' @param gene_name Optional exact gene name.
#' @return A data frame with one row per transcript.
#' @export
get_transcript_tss <- function(annotation, gene = NULL, gene_id = NULL, gene_name = NULL) {
  if (!is.null(gene_id) && !is.null(gene_name)) stop("Supply gene, gene_id, or gene_name, not multiple selectors.", call. = FALSE)
  selected <- if (inherits(gene, "GRanges")) gene else if (!is.null(gene)) find_gene(annotation, gene) else if (!is.null(gene_id)) find_gene(annotation, gene_id, "gene_id") else if (!is.null(gene_name)) find_gene(annotation, gene_name, "gene_name") else annotation
  fields <- annotation_columns(selected)
  is_transcript <- !is.na(fields$feature) & tolower(fields$feature) == "transcript"
  if (!any(is_transcript)) stop("No transcript features were found for the selected gene.", call. = FALSE)
  fields <- fields[is_transcript, , drop = FALSE]
  if (any(!fields$strand %in% c("+", "-"))) stop("Selected transcripts contain an invalid strand.", call. = FALSE)
  if (any(is.na(fields$transcript_id) | !nzchar(fields$transcript_id))) stop("One or more selected transcripts have a missing transcript_id.", call. = FALSE)
  scaffolds <- unique(fields$scaffold)
  if (length(scaffolds) > 1L) stop(sprintf("Selected gene has transcripts distributed across multiple scaffolds: %s.", paste(scaffolds, collapse = ", ")), call. = FALSE)
  fields$tss <- ifelse(fields$strand == "+", fields$transcript_start, fields$transcript_end)
  fields
}

#' Collapse transcripts sharing one transcription start site.
#' @param transcripts Data frame returned by [get_transcript_tss()].
#' @return One row per distinct gene/scaffold/strand/TSS.
#' @export
collapse_tss <- function(transcripts) {
  if (!is.data.frame(transcripts) || !all(c("gene_id", "gene_name", "transcript_id", "scaffold", "strand", "tss") %in% names(transcripts))) {
    stop("transcripts must be the data frame returned by get_transcript_tss().", call. = FALSE)
  }
  if (!nrow(transcripts)) return(data.frame())
  if (anyNA(transcripts$tss) || any(!transcripts$strand %in% c("+", "-"))) stop("Cannot collapse transcripts with missing TSS or invalid strand.", call. = FALSE)
  key <- interaction(transcripts$gene_id, transcripts$scaffold, transcripts$strand, transcripts$tss, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(transcripts)), key)
  result <- lapply(groups, function(ix) {
    row <- transcripts[ix[1L], c("gene_id", "gene_name", "scaffold", "strand", "tss"), drop = FALSE]
    ids <- unique(as.character(transcripts$transcript_id[ix]))
    row$transcript_ids <- paste(ids, collapse = ", ")
    row$transcript_id_list <- I(list(ids))
    row$n_transcripts <- length(ids)
    row
  })
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  out[order(out$scaffold, out$tss), , drop = FALSE]
}
