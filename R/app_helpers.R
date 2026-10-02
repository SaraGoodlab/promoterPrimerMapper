# Small presentation and provenance helpers shared by the Shiny views.
input_or <- function(value, default) if (is.null(value)) default else value

reference_identity <- function(path) {
  path <- trimws(input_or(path, ""))
  if (!nzchar(path)) return(list(path = "", size = NA_real_, modified = NA_real_))
  info <- file.info(path)
  list(path = normalizePath(path, mustWork = FALSE), size = unname(info$size), modified = as.numeric(info$mtime))
}

# A cheap, deterministic fingerprint for polling; never include the current time.
# Index parts and access permissions matter even when the enclosing folder is unchanged.
watched_file_state <- function(input) {
  executable <- function(path, option, program) {
    custom <- optional_path(path)
    if (!is.null(custom)) return(custom)
    configured <- getOption(option, "")
    if (is.character(configured) && length(configured) == 1L && !is.na(configured) && nzchar(configured)) return(configured)
    unname(Sys.which(program))
  }
  paths <- list(annotation = input$gtf_path, fasta = input$fasta_path, index = input$fasta_index_path,
    genome = input$bismark_genome,
    primer3 = executable(input$primer3_path, "promoterPrimerMapper.primer3", "primer3_core"),
    bowtie2 = executable(input$bowtie2_path, "promoterPrimerMapper.bowtie2", "bowtie2"))
  genome <- optional_path(input$bismark_genome)
  if (!is.null(genome)) {
    prefixes <- file.path(genome, "Bisulfite_Genome", c("CT_conversion/BS_CT", "GA_conversion/BS_GA"))
    indexes <- as.vector(outer(prefixes, c(".1", ".2", ".3", ".4", ".rev.1", ".rev.2"), paste0))
    paths$bowtie2_indexes <- c(paste0(indexes, ".bt2"), paste0(indexes, ".bt2l"))
  }
  fingerprint <- function(path) {
    identity <- reference_identity(path)
    path <- identity$path
    identity$readable <- nzchar(path) && file.access(path, 4L) == 0L
    identity$executable <- nzchar(path) && file.access(path, 1L) == 0L
    identity$directory <- nzchar(path) && dir.exists(path)
    identity
  }
  lapply(paths, function(path) lapply(input_or(path, ""), fingerprint))
}

extraction_inputs <- function(input, records = NULL) {
  choice <- suppressWarnings(as.integer(input$tss_choice))
  record <- if (length(choice) == 1L && !is.na(choice) && !is.null(records) && choice >= 1L && choice <= nrow(records))
    records[choice, c("gene_id", "scaffold", "strand", "tss"), drop = FALSE] else NULL
  if (!is.null(record)) rownames(record) <- NULL
  list(annotation = reference_identity(input$gtf_path), fasta = reference_identity(input$fasta_path),
       index = reference_identity(input$fasta_index_path), gene = trimws(input_or(input$gene_query, "")),
       tss = record, upstream = input$upstream, downstream = input$downstream)
}

design_inputs <- function(input, mode, params) {
  target <- switch(mode, island = as.character(input$bsp_island),
                   range = sort(c(input$range_from, input$range_to)), "auto")
  list(mode = mode, target = target, strand = input_or(input$bsp_strand, "both"), params = params,
       primer3 = trimws(input_or(input$primer3_path, "")))
}

has_primer_pairs <- function(design) {
  inherits(design, "bsp_design") && !is.null(design$pairs) && nrow(design$pairs) > 0L
}

specificity_labels <- function(pairs) {
  labels <- rep("Not checked", nrow(pairs))
  if (is.null(pairs$specific)) return(labels)
  checked <- !is.na(pairs$specific)
  labels[checked] <- "On-target product not found"
  if (!is.null(pairs$n_offtarget_products)) labels[checked & !is.na(pairs$n_offtarget_products) & pairs$n_offtarget_products > 0] <- "Off-target products found"
  labels[checked & pairs$specific %in% TRUE] <- "Specific within checked criteria"
  labels
}

empty_state <- function(title, text = NULL) {
  shiny::div(class = "ppm-empty", shiny::strong(title), if (!is.null(text)) shiny::p(text, class = "mb-0"))
}

metric_list <- function(values) {
  shiny::tags$dl(class = "ppm-metrics", lapply(names(values), function(name)
    shiny::div(shiny::tags$dt(name), shiny::tags$dd(as.character(values[[name]])))))
}

file_problem <- function(path, label) {
  path <- trimws(input_or(path, ""))
  if (!nzchar(path)) return(paste("Enter", label))
  if (!file.exists(path) || dir.exists(path)) return(paste(label, "file is missing"))
  if (file.access(path, 4L) != 0L) return(paste(label, "file is not readable"))
  NULL
}

tool_problem <- function(path, locate, label) {
  tryCatch({
    path <- optional_path(path)
    found <- if (is.null(path)) locate() else path
    if (!file.exists(found) || dir.exists(found) || file.access(found, 1L) != 0L) return(paste(label, "is not executable"))
    NULL
  }, error = function(e) paste(label, "is unavailable. Set its path or install it."))
}

filter_pairs <- function(pairs, target = "all", strand = "all", sort = "rank") {
  if (is.null(pairs)) return(NULL)
  if (target != "all") pairs <- pairs[pairs$target_id == target, , drop = FALSE]
  if (strand != "all") pairs <- pairs[pairs$strand == strand, , drop = FALSE]
  if (!sort %in% c("rank", "score", "product_size", "cpg_in_amplicon")) sort <- "rank"
  value <- pairs[[sort]]
  if (sort == "cpg_in_amplicon") value <- -value
  pairs[order(value, pairs$rank), , drop = FALSE]
}

pair_table_ui <- function(pairs, rank = NULL) {
  if (is.null(pairs) || !nrow(pairs)) return(empty_state("No matching pairs", "Change the target or strand filter."))
  headers <- c("Pair", "Target", "Template", "Product (bp)", "CpGs", "Tm F/R (\u00b0C)", "Score", "Specificity")
  labels <- specificity_labels(pairs)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    p <- pairs[i, ]
    selected <- identical(as.integer(p$rank), as.integer(rank))
    shiny::tags$tr(class = if (selected) "ppm-selected" else NULL,
      shiny::tags$td(shiny::tags$button(type = "button", class = "btn btn-outline-primary btn-sm",
        `data-pair-rank` = p$rank, `aria-pressed` = if (selected) "true" else "false", paste("Pair", p$rank))),
      shiny::tags$td(p$target_id), shiny::tags$td(p$strand),
      shiny::tags$td(class = "ppm-number", p$product_size), shiny::tags$td(class = "ppm-number", p$cpg_in_amplicon),
      shiny::tags$td(class = "ppm-number", sprintf("%.1f / %.1f", p$left_tm, p$right_tm)),
      shiny::tags$td(class = "ppm-number", sprintf("%.2f", p$score)), shiny::tags$td(labels[i]))
  })
  shiny::tagList(shiny::p(class = "ppm-help", sprintf("%s matching pairs. Choose a pair to update the map, details and exports.", nrow(pairs))),
    if (length(rank) == 1L && !rank %in% pairs$rank) shiny::p(class = "ppm-help", "The selected pair is outside the current filters."),
    shiny::div(class = "ppm-table-scroll ppm-pair-scroll", role = "region", `aria-label` = "Primer pair comparison", tabindex = "0",
      shiny::tags$table(class = "table table-sm ppm-pair-table",
        shiny::tags$caption(class = "visually-hidden", "Primer pairs ranked by design score"),
        shiny::tags$thead(shiny::tags$tr(lapply(headers, function(h) shiny::tags$th(scope = "col", h)))),
        shiny::tags$tbody(rows))))
}

pair_detail_ui <- function(pair) {
  oligo <- function(prefix, label) {
    value <- function(field) pair[[paste0(prefix, "_", field)]]
    shiny::div(class = "ppm-oligo", shiny::h5(label, " \u00b7 5\u2032\u21923\u2032"), shiny::tags$code(value("sequence")),
      shiny::div(class = "ppm-actions", shiny::tags$button(type = "button", class = "btn btn-outline-secondary btn-sm",
        `data-copy-oligo` = value("sequence"), paste("Copy", tolower(label))),
        shiny::span(role = "status", `aria-live` = "polite", class = "ppm-help")),
      metric_list(list("Length (nt)" = value("length"), "Tm (\u00b0C)" = sprintf("%.1f", value("tm")),
        "GC (%)" = sprintf("%.1f", value("gc")), "Converted C" = value("converted_c"), "Primer CpGs" = value("cpg"))),
      shiny::p(class = "ppm-help", paste(value("interval"), "| TSS-relative", value("tss"))),
      shiny::tags$details(class = "ppm-details", shiny::tags$summary("Thermodynamics"),
        metric_list(list("Self any (\u00b0C)" = value("self_any"), "Self end (\u00b0C)" = value("self_end"), "Hairpin (\u00b0C)" = value("hairpin")))))
  }
  shiny::tagList(shiny::h5(sprintf("Pair %s \u00b7 %s \u00b7 %s (%s)", pair$rank, pair$target_id, pair$strand, pair$assay_strand)),
    metric_list(list("Product (bp)" = pair$product_size, "Amplicon CpGs" = pair$cpg_in_amplicon, "Score" = sprintf("%.2f", pair$score))),
    shiny::p(class = "ppm-help", paste("Specificity:", specificity_labels(pair))),
    oligo("left", "Forward"), oligo("right", "Reverse"))
}
