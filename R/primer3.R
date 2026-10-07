# Thin interface to the primer3_core executable (Boulder-IO over stdin/stdout).

find_executable <- function(option, program, hint) {
  configured <- getOption(option, "")
  if (length(configured) == 1L && nzchar(configured)) {
    if (!file.exists(configured)) stop(sprintf("%s was not found at '%s' (option %s).", program, configured, option), call. = FALSE)
    return(normalizePath(configured))
  }
  found <- Sys.which(program)
  if (!nzchar(found)) stop(sprintf("%s was not found on the PATH. %s Or set options(%s = \"/path/to/%s\").", program, hint, option, program), call. = FALSE)
  unname(found)
}

#' Locate the primer3_core executable.
#'
#' Uses `getOption("promoterPrimerMapper.primer3")` when set, otherwise searches
#' the PATH.
#' @return Path to `primer3_core`.
#' @export
primer3_path <- function() {
  find_executable("promoterPrimerMapper.primer3", "primer3_core",
                  "Install Primer3 (>= 2.5.0), e.g. 'sudo apt install primer3', 'conda install -c bioconda primer3' or 'brew install primer3'.")
}

# Format one Boulder-IO record from a named list of tags.
boulder_record <- function(tags) {
  values <- vapply(tags, function(x) paste(x, collapse = " "), character(1))
  c(paste0(names(tags), "=", values), "=")
}

#' Run Primer3 on one or more Boulder-IO records.
#' @param records A named list of tags (one record) or a list of such lists.
#' @param path Path to `primer3_core`.
#' @return A list of parsed output records (named character vectors).
#' @export
run_primer3 <- function(records, path = primer3_path()) {
  if (!is.null(names(records))) records <- list(records)
  input <- tempfile(fileext = ".boulder")
  on.exit(unlink(input), add = TRUE)
  writeLines(unlist(lapply(records, boulder_record), use.names = FALSE), input)
  output <- suppressWarnings(system2(path, stdin = input, stdout = TRUE, stderr = FALSE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop(sprintf("primer3_core exited with status %s.", status), call. = FALSE)
  parse_boulder(output)
}

# run_primer3() for scoring batches. Some inputs make primer3_core abort on an
# internal assertion instead of reporting an error, which would lose the whole
# batch; on a crash the records are re-run one by one and those that crash are
# returned as PRIMER_ERROR records so callers can drop them.
run_primer3_tolerant <- function(records, path = primer3_path()) {
  output <- tryCatch(run_primer3(records, path), error = function(e) e)
  if (!inherits(output, "error")) return(output)
  if (length(records) < 2L || !grepl("exited with status", conditionMessage(output))) stop(output)
  lapply(records, function(record) {
    single <- tryCatch(run_primer3(list(record), path)[[1L]], error = function(e) NULL)
    if (is.null(single)) c(SEQUENCE_ID = as.character(record$SEQUENCE_ID), PRIMER_ERROR = "primer3_core crashed on this record") else single
  })
}

# TRUE for records that primer3_core crashed on (see run_primer3_tolerant()).
primer3_crashed <- function(record) identical(record_value(record, "PRIMER_ERROR"), "primer3_core crashed on this record")

#' Parse Boulder-IO text into records.
#' @param lines Character vector of Boulder-IO lines.
#' @return A list of named character vectors, one per record.
#' @export
parse_boulder <- function(lines) {
  ends <- which(lines == "=")
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  lapply(seq_along(ends), function(i) {
    body <- if (ends[i] > starts[i]) lines[starts[i]:(ends[i] - 1L)] else character()
    split <- regexpr("=", body, fixed = TRUE)
    stats::setNames(substring(body, split + 1L), substring(body, 1L, split - 1L))
  })
}

record_value <- function(record, key, default = NA_character_) {
  if (key %in% names(record)) record[[key]] else default
}

record_number <- function(record, key) suppressWarnings(as.numeric(record_value(record, key)))

# One oligo ("LEFT" or "RIGHT") from a Primer3 output record. Positions follow
# PRIMER_FIRST_BASE_INDEX=1; right-primer positions are converted from the
# reported 5' end to a leftmost start with primer_interval().
primer3_oligo <- function(record, side, index = 0L) {
  key <- sprintf("PRIMER_%s_%s", side, index)
  position <- as.integer(strsplit(record_value(record, key, "NA,NA"), ",", fixed = TRUE)[[1L]])
  interval <- primer_interval(position[1L], position[2L], is_right = side == "RIGHT")
  data.frame(
    start = interval[1L], end = interval[2L], length = position[2L],
    sequence = record_value(record, paste0(key, "_SEQUENCE")),
    penalty = record_number(record, paste0(key, "_PENALTY")),
    tm = record_number(record, paste0(key, "_TM")),
    gc_percent = record_number(record, paste0(key, "_GC_PERCENT")),
    self_any = record_number(record, paste0(key, "_SELF_ANY_TH")),
    self_end = record_number(record, paste0(key, "_SELF_END_TH")),
    hairpin = record_number(record, paste0(key, "_HAIRPIN_TH")),
    end_stability = record_number(record, paste0(key, "_END_STABILITY")),
    stringsAsFactors = FALSE
  )
}

primer3_error <- function(record) {
  error <- record_value(record, "PRIMER_ERROR")
  if (!is.na(error)) stop(sprintf("Primer3 error: %s", error), call. = FALSE)
}

#' Translate BSP design parameters into Primer3 global tags.
#' @param params Parameter list from [bsp_primer_params()].
#' @return A named list of Primer3 tags.
#' @export
primer3_global_tags <- function(params = bsp_primer_params()) {
  params <- validate_bsp_params(params)
  tags <- list(
    PRIMER_FIRST_BASE_INDEX = 1L,
    PRIMER_PICK_INTERNAL_OLIGO = 0L,
    PRIMER_EXPLAIN_FLAG = 1L,
    PRIMER_MIN_SIZE = params$primer_size_min,
    PRIMER_OPT_SIZE = params$primer_size_opt,
    PRIMER_MAX_SIZE = params$primer_size_max,
    PRIMER_MIN_TM = params$tm_min,
    PRIMER_OPT_TM = params$tm_opt,
    PRIMER_MAX_TM = params$tm_max,
    PRIMER_PAIR_MAX_DIFF_TM = params$max_tm_diff,
    PRIMER_MIN_GC = params$gc_min,
    PRIMER_MAX_GC = params$gc_max,
    PRIMER_MAX_POLY_X = params$max_poly_x,
    PRIMER_MAX_NS_ACCEPTED = 0L,
    PRIMER_PRODUCT_SIZE_RANGE = sprintf("%s-%s", params$product_size_min, params$product_size_max)
  )
  config <- getOption("promoterPrimerMapper.primer3_config", "")
  if (nzchar(config)) tags$PRIMER_THERMODYNAMIC_PARAMETERS_PATH <- config
  tags
}

# Evaluate candidate oligos with PRIMER_TASK=check_primers. Primers that fail
# Primer3's constraints (Tm, GC, poly-X, self-complementarity, hairpin) are
# dropped; the rest gain Primer3's penalty and thermodynamic values.
primer3_check_oligos <- function(sequences, side = c("LEFT", "RIGHT"), params, path = primer3_path()) {
  side <- match.arg(side)
  if (!length(sequences)) return(NULL)
  globals <- primer3_global_tags(params)
  records <- lapply(seq_along(sequences), function(i) {
    oligo <- if (side == "LEFT") list(SEQUENCE_PRIMER = sequences[i]) else list(SEQUENCE_PRIMER_REVCOMP = sequences[i])
    c(list(SEQUENCE_ID = i, PRIMER_TASK = "check_primers",
           PRIMER_PICK_LEFT_PRIMER = as.integer(side == "LEFT"), PRIMER_PICK_RIGHT_PRIMER = as.integer(side == "RIGHT")),
      oligo, globals)
  })
  output <- run_primer3_tolerant(records, path)
  reasons <- character()
  rows <- lapply(output, function(record) {
    if (primer3_crashed(record)) {
      reasons <<- c(reasons, "Primer3 crash")
      return(NULL)
    }
    primer3_error(record)
    if (record_number(record, sprintf("PRIMER_%s_NUM_RETURNED", side)) < 1) {
      reasons <<- c(reasons, primer3_reject_reason(record_value(record, sprintf("PRIMER_%s_EXPLAIN", side))))
      return(NULL)
    }
    oligo <- primer3_oligo(record, side)
    oligo$candidate <- as.integer(record_value(record, "SEQUENCE_ID"))
    oligo[, c("candidate", "penalty", "tm", "gc_percent", "self_any", "self_end", "hairpin", "end_stability")]
  })
  out <- do.call(rbind, rows)
  counts <- table(reasons)
  structure(if (is.null(out)) data.frame() else out, rejected = stats::setNames(as.vector(counts), names(counts)))
}

# Readable reason from a Primer3 explain string such as
# "considered 1, GC content failed 1, ok 0" (Primer3 reports the first failed check).
primer3_reject_reason <- function(explain) {
  if (is.na(explain) || !nzchar(explain)) return("rejected by Primer3")
  parts <- trimws(strsplit(explain, ",", fixed = TRUE)[[1L]])
  parts <- parts[!grepl("^(considered|ok) [0-9]+$", parts)]
  if (!length(parts)) return("rejected by Primer3")
  reason <- sub(" [0-9]+$", "", parts[1L])
  labels <- c("GC content failed" = "GC out of range", "low tm" = "Tm too low", "high tm" = "Tm too high",
              "long poly-x seq" = "homopolymer", "high hairpin stability" = "hairpin", "high any compl" = "self-dimer",
              "high end compl" = "3' self-dimer", "high template mispriming score" = "template mispriming",
              "too many Ns" = "N in primer", "GC clamp failed" = "GC clamp")
  if (reason %in% names(labels)) unname(labels[[reason]]) else reason
}

# Check primer pairs (cross-dimers, Tm difference, product) against a template;
# `template` is one sequence or one per pair.
primer3_check_pairs <- function(template, left, right, params, path = primer3_path()) {
  globals <- primer3_global_tags(params)
  template <- rep_len(template, length(left))
  records <- lapply(seq_along(left), function(i) {
    c(list(SEQUENCE_ID = i, SEQUENCE_TEMPLATE = template[i], SEQUENCE_PRIMER = left[i], SEQUENCE_PRIMER_REVCOMP = right[i],
           PRIMER_TASK = "check_primers", PRIMER_PICK_LEFT_PRIMER = 1L, PRIMER_PICK_RIGHT_PRIMER = 1L,
           # Primer3 emits product Tm only with a non-default product Tm bound.
           # This nonbinding bound requests the value without filtering pairs.
           PRIMER_PRODUCT_MIN_TM = -999999),
      globals)
  })
  output <- run_primer3_tolerant(records, path)
  rows <- lapply(output, function(record) {
    if (primer3_crashed(record)) return(NULL)
    primer3_error(record)
    if (record_number(record, "PRIMER_PAIR_NUM_RETURNED") < 1) return(NULL)
    data.frame(
      candidate = as.integer(record_value(record, "SEQUENCE_ID")),
      pair_penalty = record_number(record, "PRIMER_PAIR_0_PENALTY"),
      compl_any = record_number(record, "PRIMER_PAIR_0_COMPL_ANY_TH"),
      compl_end = record_number(record, "PRIMER_PAIR_0_COMPL_END_TH"),
      product_tm = record_number(record, "PRIMER_PAIR_0_PRODUCT_TM")
    )
  })
  do.call(rbind, rows)
}
