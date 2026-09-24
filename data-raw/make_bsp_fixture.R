# Regenerates tests/testthat/fixtures/bsp_reference.{fa,fa.fai,gtf}: two synthetic promoters
# with a CpG island around the TSS and CpG-depleted flanks, one per strand.
# Run from the package root: Rscript data-raw/make_bsp_fixture.R
set.seed(42)

simulate <- function(n, gc, cpg_keep) {
  bases <- sample(c("A", "C", "G", "T"), n, replace = TRUE, prob = c(1 - gc, gc, gc, 1 - gc) / 2)
  # Deplete CpGs as in vertebrate genomes (methyl-C deamination to TpG).
  cpg <- which(bases[-n] == "C" & bases[-1L] == "G")
  lost <- cpg[stats::runif(length(cpg)) > cpg_keep]
  bases[lost] <- "T"
  paste(bases, collapse = "")
}

promoter <- function(length, island_start, island_end) {
  paste0(simulate(island_start - 1L, 0.40, 0.2),
         simulate(island_end - island_start + 1L, 0.64, 0.8),
         simulate(length - island_end, 0.40, 0.2))
}

scaffolds <- c(chrTest1 = promoter(6000L, 2500L, 3400L), chrTest2 = promoter(5000L, 2100L, 2900L))
fasta <- "tests/testthat/fixtures/bsp_reference.fa"
lines <- unlist(lapply(names(scaffolds), function(name) {
  s <- scaffolds[[name]]
  c(paste0(">", name), substring(s, seq(1L, nchar(s), 60L), pmin(seq(60L, nchar(s) + 59L, 60L), nchar(s))))
}))
writeLines(lines, fasta)
Rsamtools::indexFa(fasta)

gtf <- c(
  'chrTest1\tsynthetic\tgene\t3001\t5800\t.\t+\t.\tgene_id "BSPPLUS"; gene_name "BSPPLUS";',
  'chrTest1\tsynthetic\ttranscript\t3001\t5800\t.\t+\t.\tgene_id "BSPPLUS"; gene_name "BSPPLUS"; transcript_id "BSPPLUS-201";',
  'chrTest1\tsynthetic\ttranscript\t3001\t4900\t.\t+\t.\tgene_id "BSPPLUS"; gene_name "BSPPLUS"; transcript_id "BSPPLUS-202";',
  'chrTest2\tsynthetic\tgene\t300\t2500\t.\t-\t.\tgene_id "BSPMINUS"; gene_name "BSPMINUS";',
  'chrTest2\tsynthetic\ttranscript\t300\t2500\t.\t-\t.\tgene_id "BSPMINUS"; gene_name "BSPMINUS"; transcript_id "BSPMINUS-201";'
)
writeLines(gtf, "tests/testthat/fixtures/bsp_reference.gtf")
