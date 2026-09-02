# In-silico bisulfite PCR: align primers to a Bismark-prepared converted genome
# with Bowtie2 and predict every product the pair could amplify.
#
# Bismark's CT index holds the C->T converted top strand; its GA index holds the
# G->A converted top strand, i.e. the converted bottom strand in top-strand
# orientation. On either index a product needs a primer aligned forward
# followed, within the product limit, by a primer aligned in reverse.

#' Locate the bowtie2 executable.
#'
#' Uses `getOption("promoterPrimerMapper.bowtie2")` when set, otherwise searches
#' the PATH.
#' @return Path to `bowtie2`.
#' @export
bowtie2_path <- function() {
  find_executable("promoterPrimerMapper.bowtie2", "bowtie2",
                  "Install Bowtie2, e.g. 'sudo apt install bowtie2' or 'conda install -c bioconda bowtie2'.")
}

bismark_indexes <- function(bismark_genome) {
  if (length(bismark_genome) != 1L || !nzchar(bismark_genome) || !dir.exists(bismark_genome)) {
    stop("Bismark genome folder does not exist.", call. = FALSE)
  }
  base <- file.path(bismark_genome, "Bisulfite_Genome")
  indexes <- c(CT = file.path(base, "CT_conversion", "BS_CT"), GA = file.path(base, "GA_conversion", "BS_GA"))
  present <- file.exists(paste0(indexes, ".1.bt2")) | file.exists(paste0(indexes, ".1.bt2l"))
  if (!all(present)) {
    stop(sprintf("No Bowtie2 bisulfite indexes under %s. Run 'bismark_genome_preparation --bowtie2 %s' first.", base, bismark_genome), call. = FALSE)
  }
  indexes
}

# 1-based reference offsets of mismatches encoded in a SAM MD:Z tag.
md_mismatch_offsets <- function(md) {
  tokens <- regmatches(md, gregexpr("[0-9]+|\\^[A-Z]+|[A-Z]", md))[[1L]]
  position <- 0L
  offsets <- integer()
  for (token in tokens) {
    if (grepl("^[0-9]+$", token)) {
      position <- position + as.integer(token)
    } else if (startsWith(token, "^")) {
      position <- position + nchar(token) - 1L
    } else {
      position <- position + 1L
      offsets <- c(offsets, position)
    }
  }
  offsets
}

#' Parse Bowtie2 SAM output for primer alignments.
#' @param lines SAM lines.
#' @param index Label of the index the lines came from (`CT` or `GA`).
#' @param three_prime Number of 3'-terminal primer bases in which a mismatch
#'   is taken to prevent extension.
#' @return A data frame of aligned hits with chromosome (Bismark conversion
#'   suffix removed), 1-based interval, orientation, mismatch count, and
#'   whether a mismatch falls in the 3'-terminal bases.
#' @export
parse_bowtie2_sam <- function(lines, index, three_prime = 3L) {
  lines <- lines[nzchar(lines) & !startsWith(lines, "@")]
  empty <- data.frame(primer = character(), index = character(), chrom = character(), start = integer(), end = integer(),
                      reverse = logical(), mismatches = integer(), three_prime_mismatch = logical(), stringsAsFactors = FALSE)
  if (!length(lines)) return(empty)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  flag <- as.integer(vapply(fields, `[`, character(1), 2L))
  mapped <- bitwAnd(flag, 4L) == 0L
  if (!any(mapped)) return(empty)
  fields <- fields[mapped]
  flag <- flag[mapped]
  tag <- function(f, name) {
    hit <- grep(paste0("^", name, ":"), f[-(1:11)], value = TRUE)
    if (length(hit)) sub("^[A-Z0-9]+:[A-Za-z]:", "", hit[1L]) else NA_character_
  }
  length <- nchar(vapply(fields, `[`, character(1), 10L))
  start <- as.integer(vapply(fields, `[`, character(1), 4L))
  reverse <- bitwAnd(flag, 16L) > 0L
  three_prime_mismatch <- vapply(seq_along(fields), function(i) {
    offsets <- md_mismatch_offsets(tag(fields[[i]], "MD"))
    # The primer's 3' end is the alignment end when forward, the start when reverse.
    if (reverse[i]) any(offsets <= three_prime) else any(offsets > length[i] - three_prime)
  }, logical(1))
  data.frame(
    primer = vapply(fields, `[`, character(1), 1L),
    index = index,
    chrom = sub("_(CT|GA)_converted$", "", vapply(fields, `[`, character(1), 3L)),
    start = start, end = start + length - 1L, reverse = reverse,
    mismatches = as.integer(vapply(fields, tag, character(1), name = "NM")),
    three_prime_mismatch = three_prime_mismatch,
    stringsAsFactors = FALSE
  )
}

#' Predict PCR products from primer hits.
#' @param hits Hits from [parse_bowtie2_sam()], already filtered to primers
#'   able to prime.
#' @param max_product Maximum product length.
#' @return A data frame of products: index, chromosome, interval, length and
#'   the forward/reverse primers that generate it.
#' @export
predict_products <- function(hits, max_product = 3000L) {
  forward <- hits[!hits$reverse, , drop = FALSE]
  reverse <- hits[hits$reverse, , drop = FALSE]
  empty <- data.frame(index = character(), chrom = character(), start = integer(), end = integer(), length = integer(),
                      forward_primer = character(), reverse_primer = character(), stringsAsFactors = FALSE)
  if (!nrow(forward) || !nrow(reverse)) return(empty)
  joined <- merge(forward, reverse, by = c("index", "chrom"), suffixes = c(".f", ".r"))
  joined <- joined[joined$end.r > joined$end.f & joined$start.r >= joined$start.f &
                     joined$end.r - joined$start.f + 1L <= max_product, , drop = FALSE]
  if (!nrow(joined)) return(empty)
  data.frame(index = joined$index, chrom = joined$chrom, start = joined$start.f, end = joined$end.r,
             length = joined$end.r - joined$start.f + 1L, forward_primer = joined$primer.f, reverse_primer = joined$primer.r,
             stringsAsFactors = FALSE)
}

# Primer sequences for alignment: degenerate CpG bases take the unmethylated base.
alignable_sequence <- function(sequence) chartr("YR", "TA", toupper(sequence))

#' Check BSP primer pairs for off-target products in a bisulfite-converted genome.
#'
#' Primers are aligned with Bowtie2 to both converted strands of a genome
#' prepared with `bismark_genome_preparation --bowtie2`. A hit can prime when it
#' has at most `max_mismatch` mismatches and none in the 3'-terminal three
#' bases. Every forward/reverse combination of a pair's primers (including a
#' primer with itself) within `max_product` on the same converted strand is a
#' predicted product.
#' @param design A `bsp_design` object from [design_bsp_primers()].
#' @param bismark_genome Folder containing `Bisulfite_Genome/`.
#' @param top_n Number of best-ranked pairs to check per target.
#' @param max_mismatch Maximum mismatches for a priming hit.
#' @param max_product Longest product considered.
#' @param max_hits Maximum alignments reported per primer and index.
#' @param threads Bowtie2 threads.
#' @param bowtie2 Path to `bowtie2`.
#' @return The design with specificity columns added to `pairs`
#'   (`n_hits_left`, `n_hits_right`, `n_products`, `on_target`,
#'   `n_offtarget_products`, `specific`; `NA` for unchecked pairs) and a
#'   `specificity` element holding all hits and products.
#' @export
check_bsp_specificity <- function(design, bismark_genome, top_n = 10L, max_mismatch = 2L, max_product = 3000L,
                                  max_hits = 100L, threads = 1L, bowtie2 = bowtie2_path()) {
  pairs <- design_pairs(design)
  top_n <- assert_scalar_integer(top_n, "top_n", 1L)
  max_mismatch <- assert_scalar_integer(max_mismatch, "max_mismatch", 0L)
  max_product <- assert_scalar_integer(max_product, "max_product", 1L)
  indexes <- bismark_indexes(bismark_genome)
  group <- if (is.null(pairs$target_id)) rep("", nrow(pairs)) else pairs$target_id
  checked <- pairs[unlist(lapply(split(seq_len(nrow(pairs)), factor(group, levels = unique(group))), utils::head, top_n)), , drop = FALSE]
  fasta <- tempfile(fileext = ".fa")
  on.exit(unlink(fasta), add = TRUE)
  names <- c(sprintf("pair%s_left", checked$rank), sprintf("pair%s_right", checked$rank))
  writeLines(as.vector(rbind(paste0(">", names), alignable_sequence(c(checked$left_sequence, checked$right_sequence)))), fasta)
  args <- c("-f", "-U", fasta, "--end-to-end", "-k", max_hits, "-L", 10L, "-N", 1L, "-i", "S,1,0.5",
            "--mp", "6,6", "--rdg", "100,100", "--rfg", "100,100", "--score-min", sprintf("C,%s,0", -6L * max_mismatch),
            "--no-unal", "--no-hd", "-p", threads)
  hits <- do.call(rbind, lapply(names(indexes), function(index) {
    sam <- suppressWarnings(system2(bowtie2, c("-x", indexes[[index]], args), stdout = TRUE, stderr = FALSE))
    status <- attr(sam, "status")
    if (!is.null(status) && status != 0L) stop(sprintf("bowtie2 exited with status %s on the %s index.", status, index), call. = FALSE)
    parse_bowtie2_sam(sam, index)
  }))
  priming <- hits[hits$mismatches <= max_mismatch & !hits$three_prime_mismatch, , drop = FALSE]
  summaries <- lapply(seq_len(nrow(checked)), function(i) {
    p <- checked[i, ]
    primers <- sprintf("pair%s_%s", p$rank, c("left", "right"))
    products <- predict_products(priming[priming$primer %in% primers, , drop = FALSE], max_product)
    products$rank <- rep(p$rank, nrow(products))
    on_target <- products$chrom == p$scaffold & products$start == p$amplicon_start & products$end == p$amplicon_end
    list(summary = data.frame(rank = p$rank,
                              n_hits_left = sum(priming$primer == primers[1L]), n_hits_right = sum(priming$primer == primers[2L]),
                              n_products = nrow(products), on_target = any(on_target), n_offtarget_products = sum(!on_target)),
         products = products)
  })
  summary <- do.call(rbind, lapply(summaries, `[[`, "summary"))
  summary$specific <- summary$on_target & summary$n_offtarget_products == 0L
  columns <- setdiff(names(summary), "rank")
  previous <- design$pairs[setdiff(names(design$pairs), columns)]
  design$pairs <- cbind(previous, summary[match(previous$rank, summary$rank), columns])
  rownames(design$pairs) <- NULL
  design$specificity <- list(hits = hits, products = do.call(rbind, lapply(summaries, `[[`, "products")),
                             settings = list(max_mismatch = max_mismatch, max_product = max_product, top_n = top_n))
  design
}
