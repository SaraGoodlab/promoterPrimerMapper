make_annotation <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("scaffold1", 5),
    ranges = IRanges::IRanges(c(40, 50, 50, 55, 70), c(90, 70, 75, 80, 80)),
    strand = c("+", "+", "+", "+", "-"),
    type = c("gene", "transcript", "transcript", "transcript", "transcript"),
    gene_id = c("geneA", "geneA", "geneA", "geneA", "geneB"),
    gene_name = c("Gene A", "Gene A", "Gene A", "Gene A", "Gene B"),
    transcript_id = c(NA, "txA1", "txA2", "txA3", "txB1")
  )
  gr
}

test_that("TSS is inferred from transcript boundaries", {
  records <- get_transcript_tss(make_annotation(), gene_id = "geneA")
  expect_equal(records$tss, c(50L, 50L, 55L))
  minus <- get_transcript_tss(make_annotation(), gene_id = "geneB")
  expect_equal(minus$tss, 80L)
})

test_that("transcripts with identical TSS collapse and distinct TSSs remain", {
  collapsed <- collapse_tss(get_transcript_tss(make_annotation(), gene_id = "geneA"))
  expect_equal(nrow(collapsed), 2L)
  same <- collapsed[collapsed$tss == 50L, , drop = FALSE]
  expect_equal(same$n_transcripts, 2L)
  expect_equal(same$transcript_ids, "txA1, txA2")
  expect_equal(collapsed$tss, c(50L, 55L))
})

test_that("plus and minus extraction are transcription-oriented", {
  reference <- c(scaffold1 = paste(rep("ACGT", 25), collapse = ""))
  plus_record <- data.frame(gene_id = "geneA", gene_name = "Gene A", scaffold = "scaffold1", strand = "+", tss = 50L, transcript_ids = "txA1")
  minus_record <- data.frame(gene_id = "geneB", gene_name = "Gene B", scaffold = "scaffold1", strand = "-", tss = 80L, transcript_ids = "txB1")
  plus <- extract_tss_region(plus_record, reference, 2L, 3L)
  minus <- extract_tss_region(minus_record, reference, 2L, 3L)
  expect_equal(c(plus$region_start, plus$region_end, plus$sequence_length), c(48L, 53L, 6L))
  expect_equal(plus$sequence, substr(reference[[1]], 48, 53))
  expect_equal(c(minus$region_start, minus$region_end), c(77L, 82L))
  expect_equal(minus$sequence, as.character(Biostrings::reverseComplement(Biostrings::DNAString(substr(reference[[1]], 77, 82)))))
})

test_that("four FASTA sequence variants are generated", {
  sequence <- structure(list(scaffold = "scaffold1", strand = "+", tss = 50L, region_start = 48L,
                             region_end = 53L, upstream = 2L, downstream = 3L,
                             sequence = "ACCGTG", sequence_length = 6L), class = c("promoter_sequence", "list"))
  variants <- sequence_variants(sequence)
  expect_equal(names(variants), c("requested", "complementary", "bisulfite_requested", "bisulfite_complementary"))
  expect_equal(unname(variants), c("ACCGTG", "CACGGT", "ATTGTG", "TATGGT"))
  expect_equal(bisulfite_convert("CcN"), "TTN")
  output <- tempfile(fileext = ".fa")
  write_sequence_fasta(sequence, "bisulfite_complementary", output)
  expect_equal(readLines(output), c(">unknown_gene|scaffold1:48-53|strand=+|TSS=50|upstream=2|downstream=3|variant=bisulfite_complementary", "TATGGT"))
})

test_that("relative intervals map on both strands and respect boundaries", {
  expect_equal(map_relative_interval(1, 1, 100, 110, "+")$genomic_start, 100L)
  expect_equal(map_relative_interval(11, 11, 100, 110, "+")$genomic_end, 110L)
  expect_equal(map_relative_interval(1, 1, 100, 110, "-")$genomic_start, 110L)
  expect_equal(map_relative_interval(11, 11, 100, 110, "-")$genomic_end, 100L)
  expect_equal(map_relative_interval(3, 5, 100, 110, "-"), data.frame(genomic_start = 106L, genomic_end = 108L))
  expect_error(map_relative_interval(1, 12, 100, 110, "+", 11L), "outside")
})

test_that("primer starts normalize explicitly and primer intervals map without alignment", {
  expect_equal(normalize_primer_start(1, "one-based"), 1L)
  expect_equal(normalize_primer_start(0, "zero-based"), 1L)
  sequence <- structure(list(scaffold = "scaffold1", strand = "+", tss = 50L, region_start = 48L,
                             region_end = 53L, sequence = "ACGTAC", sequence_length = 6L), class = c("promoter_sequence", "list"))
  primers <- data.frame(primer_type = c("Left primer", "Right primer"), start = c(1L, 3L), size = c(2L, 3L))
  mapped <- map_primer_coordinates(primers, sequence)
  expect_equal(mapped$relative_end, c(2L, 5L))
  expect_equal(mapped$genomic_interval, c("scaffold1:48-49", "scaffold1:50-52"))
  zero <- map_primer_coordinates(data.frame(primer_type = "Left primer", start = 0L, size = 1L), sequence, "zero-based")
  expect_equal(zero$relative_start, 1L)
})

test_that("pasted primer output ignores optional columns", {
  pasted <- "Left primer  1074  22  59.36  54.55  4  TTATTT\nRight primer 1299  28  54.70  50.00  10 AAAA\nProduct size: 226"
  parsed <- parse_primer_output(pasted)
  expect_equal(parsed$primer_type, c("Left primer", "Right primer"))
  expect_equal(parsed$start, c(1074L, 1299L))
  expect_equal(parsed$size, c(22L, 28L))
})

test_that("synthetic GTF and FASTA are retained as an integration fixture", {
  gtf_path <- system.file("extdata", "synthetic.gtf", package = "promoterPrimerMapper")
  fasta_path <- system.file("extdata", "synthetic.fa", package = "promoterPrimerMapper")
  expect_true(file.exists(gtf_path))
  expect_true(file.exists(fasta_path))
  annotation <- load_annotation(gtf_path)
  expect_equal(nrow(collapse_tss(get_transcript_tss(annotation, gene_name = "Gene A"))), 2L)
  local_fasta <- tempfile(fileext = ".fa")
  file.copy(fasta_path, local_fasta)
  Rsamtools::indexFa(local_fasta)
  record <- collapse_tss(get_transcript_tss(annotation, gene_id = "geneB"))[1, , drop = FALSE]
  extracted <- extract_tss_region(record, local_fasta, 2L, 2L)
  expect_equal(extracted$sequence_length, 5L)
})
