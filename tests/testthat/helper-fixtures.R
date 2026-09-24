# Shared fixtures and skips for the test suite.

promoter <- function(sequence, strand = "+", region_start = 100L, tss = NULL) {
  region_end <- region_start + nchar(sequence) - 1L
  if (is.null(tss)) tss <- if (strand == "+") region_start else region_end
  structure(list(gene_id = "g", scaffold = "s1", strand = strand, tss = tss, region_start = region_start,
                 region_end = region_end, sequence = sequence), class = c("promoter_sequence", "list"))
}

fixture_reference <- function(file) test_path("fixtures", paste0("bsp_reference.", file))

fixture_promoter <- function(gene, upstream = 1000L, downstream = 500L) {
  annotation <- load_annotation(fixture_reference("gtf"))
  record <- collapse_tss(get_transcript_tss(annotation, gene = gene))[1L, , drop = FALSE]
  extract_tss_region(record, fixture_reference("fa"), upstream, downstream)
}

tool_available <- function(locate) !inherits(try(locate(), silent = TRUE), "try-error")

skip_if_no_primer3 <- function() skip_if_not(tool_available(primer3_path), "primer3_core is not available")

skip_if_no_bismark <- function() {
  skip_if_not(tool_available(bowtie2_path), "bowtie2 is not available")
  skip_if_not(nzchar(Sys.which("bismark_genome_preparation")), "bismark_genome_preparation is not available")
}
