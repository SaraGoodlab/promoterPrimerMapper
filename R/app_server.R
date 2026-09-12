# Read every BSP parameter input back into a parameter list. Inputs that have
# not been rendered yet fall back to the bsp_primer_params() defaults.
collect_bsp_params <- function(input) {
  defaults <- bsp_primer_params()
  stats::setNames(lapply(names(defaults), function(name) {
    value <- input[[bsp_param_input_id(name)]]
    if (is.null(value)) defaults[[name]] else value
  }), names(defaults))
}

optional_path <- function(value) {
  value <- trimws(if (is.null(value)) "" else value)
  if (nzchar(value)) value else NULL
}

# Sense interval of a TSS-relative From/To range; errors with a readable message.
range_interval <- function(extracted, from, to) {
  if (!is.numeric(from) || !is.numeric(to) || length(from) != 1L || length(to) != 1L || is.na(from) || is.na(to)) {
    stop("Enter both ends of the range.", call. = FALSE)
  }
  tss <- tss_sense_position(extracted)
  bounds <- sort(sense_position_from_tss(round(c(from, to)), tss))
  n <- validate_extracted_sequence(extracted)$sequence_length
  if (bounds[1L] < 1L || bounds[2L] > n) {
    stop(sprintf("The range must lie within the extracted sequence (%s).",
                 format_tss_interval(tss_relative_position(1L, tss), tss_relative_position(n, tss))), call. = FALSE)
  }
  as.integer(bounds)
}

# Targets for a GUI target mode: all islands, one selected island, or the selected range.
mode_targets <- function(extracted, mode, params, island = NULL, from = NULL, to = NULL) {
  switch(mode,
    auto = island_targets(extracted, params),
    island = {
      k <- suppressWarnings(as.integer(island))
      if (length(k) != 1L || is.na(k)) stop("No CpG island is selected (the sequence may have none); choose another target.", call. = FALSE)
      island_targets(extracted, params, k)
    },
    range = {
      range_interval(extracted, from, to)
      resolve_bsp_target(extracted, c(from, to), params)
    },
    stop(sprintf("Unknown target mode '%s'.", mode), call. = FALSE))
}

# Status text for a finished design.
design_message <- function(design) {
  targets <- design$targets
  n <- if (is.null(design$pairs)) 0L else nrow(design$pairs)
  if (!nrow(targets)) {
    return("The sequence has no CpG island and no CpG cluster outside islands to target; select a range on the map instead.")
  }
  advice <- if (length(design_hints(design))) "See the hints under Primer pairs for what blocked the design." else
    "See the candidate counts under Primer pairs; relax the parameters, choose another target or widen the extraction window."
  if (!n) {
    no_fallback <- if (isTRUE(design$fallback) && !any(targets$kind == "non_island"))
      "No CpG cluster outside the islands qualified for the non-island fallback." else NULL
    return(paste(c(sprintf("No primer pair satisfied all constraints for %s.", paste(targets$target_id, collapse = ", ")), no_fallback, advice), collapse = " "))
  }
  found <- targets[targets$n_pairs > 0L, , drop = FALSE]
  text <- sprintf("Designed %s BSP primer pair%s for %s.", n, if (n == 1L) "" else "s", paste(target_lines(found), collapse = "; "))
  if (isTRUE(design$fallback)) text <- paste("No CpG island yielded a primer pair, so non-island CpG regions were targeted.", text)
  empty <- targets$target_id[targets$n_pairs == 0L & !(isTRUE(design$fallback) & targets$kind == "island")]
  if (length(empty)) text <- paste(text, sprintf("No pair for %s.", paste(empty, collapse = ", ")), advice)
  text
}

pair_choices <- function(pairs) {
  if (is.null(pairs) || !nrow(pairs)) return(character(0))
  stats::setNames(pairs$rank, sprintf("Pair %s \u2013 %s (%s bp, %s CpGs, %s strand)", pairs$rank, pairs$target_id,
                                      pairs$product_size, pairs$cpg_in_amplicon, pairs$strand))
}

app_server <- function(input, output, session) {
  annotation_cache <- shiny::reactiveVal(NULL)
  tss_records <- shiny::reactiveVal(NULL)
  selected_tss <- shiny::reactiveVal(NULL)
  extracted <- shiny::reactiveVal(NULL)
  design <- shiny::reactiveVal(NULL)
  message_text <- shiny::reactiveVal(NULL)

  output$status <- shiny::renderUI({
    value <- message_text()
    if (is.null(value)) return(NULL)
    shiny::div(class = "alert alert-info", value)
  })

  # Re-import the GTF only when its path or modification time changes.
  cached_annotation <- function(path) {
    key <- paste(normalizePath(path, mustWork = FALSE), file.mtime(path))
    cache <- annotation_cache()
    if (!is.null(cache) && identical(cache$key, key)) return(cache$value)
    value <- load_annotation(path)
    annotation_cache(list(key = key, value = value))
    value
  }

  shiny::observeEvent(input$find_gene, {
    shiny::req(input$gtf_path, input$gene_query)
    tryCatch({
      shiny::withProgress(message = "Finding gene and TSSs", value = 0, {
        shiny::incProgress(0.35, detail = "Loading GTF annotation")
        gtf <- cached_annotation(trimws(input$gtf_path))
        shiny::incProgress(0.35, detail = "Calculating transcript TSSs")
        records <- collapse_tss(get_transcript_tss(gtf, gene = input$gene_query))
        shiny::incProgress(0.2, detail = "Preparing TSS choices")
        tss_records(records)
        choices <- stats::setNames(seq_len(nrow(records)), sprintf("%s:%s (%s, %s)", records$scaffold, records$tss, records$strand, records$transcript_ids))
        shiny::updateSelectInput(session, "tss_choice", choices = choices, selected = if (nrow(records) == 1L) "1" else character(0))
        shiny::incProgress(0.1)
        message_text(sprintf("Found %s distinct TSS%s.", nrow(records), if (nrow(records) == 1L) "" else "s"))
      })
    }, error = function(e) message_text(conditionMessage(e)))
  })

  output$tss_table <- shiny::renderTable({
    records <- tss_records()
    if (is.null(records)) return(NULL)
    records[, c("tss", "scaffold", "strand", "n_transcripts", "transcript_ids"), drop = FALSE]
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  fasta_input <- shiny::reactive({
    shiny::req(input$fasta_path, input$fasta_index_path)
    ensure_fasta_index(trimws(input$fasta_path), trimws(input$fasta_index_path))
  })

  bsp_params <- shiny::reactive(collect_bsp_params(input))
  # Display helpers keep working with defaults while a parameter input is invalid.
  display_params <- shiny::reactive(tryCatch(validate_bsp_params(bsp_params()), error = function(e) bsp_primer_params()))
  target_mode <- shiny::reactive(if (is.null(input$bsp_target_mode)) "auto" else input$bsp_target_mode)

  # Design the current extraction; returns the status text.
  run_design <- function(value, mode) {
    tryCatch({
      params <- validate_bsp_params(bsp_params(), value$sequence_length)
      target <- if (identical(mode, "auto")) "auto" else mode_targets(value, mode, params, input$bsp_island, input$range_from, input$range_to)
      primer3 <- optional_path(input$primer3_path)
      strand <- if (is.null(input$bsp_strand)) "both" else input$bsp_strand
      result <- shiny::withProgress(message = "Designing BSP primers", value = 0, {
        design_bsp_primers(value, strand, target, params, primer3 = if (is.null(primer3)) primer3_path() else primer3,
                           progress = function(step, steps, detail) shiny::setProgress((step - 1) / steps, detail = detail))
      })
      design(result)
      shiny::updateSelectInput(session, "bsp_pair", choices = pair_choices(result$pairs))
      design_message(result)
    }, error = function(e) paste("Primer design failed:", conditionMessage(e)))
  }

  shiny::observeEvent(input$extract, {
    shiny::req(tss_records(), input$tss_choice, input$upstream, input$downstream)
    value <- tryCatch({
      shiny::withProgress(message = "Extracting reference sequence", value = 0, {
        shiny::incProgress(0.2, detail = "Opening FASTA and index")
        record <- tss_records()[as.integer(input$tss_choice), , drop = FALSE]
        reference <- fasta_input()
        shiny::incProgress(0.3, detail = "Reading genomic interval")
        value <- extract_tss_region(record, reference, input$upstream, input$downstream)
        selected_tss(record)
        value
      })
    }, error = function(e) {
      message_text(conditionMessage(e))
      NULL
    })
    if (is.null(value)) return()
    extracted(value)
    design(NULL)
    shiny::updateSelectInput(session, "bsp_pair", choices = character(0))
    islands <- island_table(value, display_params())
    nearest <- nearest_island(cpg_island_scan(value), tss_sense_position(value))
    shiny::updateSelectInput(session, "bsp_island", choices = island_choices(islands),
                             selected = if (is.na(nearest)) character(0) else as.character(nearest))
    shiny::updateRadioButtons(session, "bsp_target_mode", selected = "auto")
    extracted_text <- sprintf("Extracted %s bp from %s:%s-%s.", value$sequence_length, value$scaffold, value$region_start, value$region_end)
    message_text(paste(extracted_text, run_design(value, "auto")))
  })

  # CpG map and selection ---------------------------------------------------------

  # Selected range as a sense interval, or a message explaining why it is invalid.
  selection <- shiny::reactive({
    value <- extracted()
    if (is.null(value)) return(NULL)
    tryCatch(range_interval(value, input$range_from, input$range_to), error = function(e) conditionMessage(e))
  })

  planned_targets <- shiny::reactive({
    value <- extracted()
    shiny::req(value)
    tryCatch(mode_targets(value, target_mode(), display_params(), input$bsp_island, input$range_from, input$range_to),
             error = function(e) conditionMessage(e))
  })

  shiny::observeEvent(input$cpg_brush, {
    value <- extracted()
    brush <- input$cpg_brush
    if (is.null(value) || is.null(brush$xmin) || is.null(brush$xmax)) return()
    relative <- tss_relative_position(brush_to_sense(brush$xmin, brush$xmax, value$sequence_length), tss_sense_position(value))
    shiny::updateNumericInput(session, "range_from", value = relative[1L])
    shiny::updateNumericInput(session, "range_to", value = relative[2L])
  })

  shiny::observeEvent(input$cpg_click, {
    value <- extracted()
    click <- input$cpg_click
    if (is.null(value) || is.null(click$x)) return()
    k <- island_at(cpg_island_scan(value), click$x, click$y)
    if (is.na(k)) return()
    shiny::updateSelectInput(session, "bsp_island", selected = as.character(k))
    shiny::updateRadioButtons(session, "bsp_target_mode", selected = "island")
    message_text(sprintf("Selected CpG island CGI%s as the target; press 'Design BSP primers' to design it.", k))
  })

  shiny::observeEvent(input$use_selection, {
    interval <- selection()
    if (!is.numeric(interval)) {
      message_text(if (is.character(interval)) interval else "Extract a sequence and select a range first.")
      return()
    }
    shiny::updateRadioButtons(session, "bsp_target_mode", selected = "range")
    tss <- tss_sense_position(extracted())
    message_text(sprintf("Target set to the selected range %s; press 'Design BSP primers' to design it.",
                         format_tss_interval(tss_relative_position(interval[1L], tss), tss_relative_position(interval[2L], tss))))
  })

  output$region_summary <- shiny::renderUI({
    value <- extracted()
    if (is.null(value)) return(shiny::helpText("Find a gene and extract its promoter to see the CpG map."))
    summary <- cpg_summary(value)
    n_islands <- nrow(cpg_island_scan(value))
    shiny::p(class = "small mb-1", shiny::strong(design_label(list(extracted = value))),
             sprintf(" %s:%s-%s (%s strand), TSS %s:%s \u2022 %s bp, GC %s%%, %s CpGs, CpG o/e %s, %s CpG island%s",
                     value$scaffold, value$region_start, value$region_end, value$strand, value$scaffold, value$tss,
                     summary$length, summary$gc_percent, summary$cpg_count, summary$cpg_obs_exp, n_islands, if (n_islands == 1L) "" else "s"))
  })

  output$cpg_map <- shiny::renderPlot({
    value <- extracted()
    shiny::req(value)
    current <- design()
    rank <- suppressWarnings(as.integer(input$bsp_pair))
    targets <- NULL
    if (is.null(current)) {
      planned <- planned_targets()
      if (is.data.frame(planned)) targets <- planned
    }
    interval <- selection()
    island <- if (identical(target_mode(), "island")) suppressWarnings(as.integer(input$bsp_island)) else NULL
    plot_cpg_map(value, current, if (length(rank) == 1L && !is.na(rank)) rank else NULL,
                 targets = targets, selection = if (is.numeric(interval)) interval else NULL, island = island)
  })

  output$cursor_position <- shiny::renderText({
    value <- extracted()
    hover <- input$cpg_hover
    if (is.null(value)) return("")
    if (is.null(hover$x)) return("Hover over the map to read TSS-relative and genomic positions.")
    describe_position(value, hover$x)
  })

  output$selection_table <- shiny::renderTable({
    interval <- selection()
    if (is.null(interval)) return(NULL)
    if (is.character(interval)) return(data.frame(note = interval))
    region <- describe_region(extracted(), interval[1L], interval[2L])
    region$gc_percent <- sprintf("%.1f", region$gc_percent)
    stats::setNames(region, c("TSS-relative", "Sequence", "Genomic", "Length (bp)", "CpGs", "GC %", "CpG islands"))
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$cpg_islands <- shiny::renderTable({
    shiny::req(extracted())
    islands <- island_table(extracted(), display_params())
    if (!nrow(islands)) {
      return(data.frame(note = "No CpG island predicted (100 bp windows, GC >= 50%, CpG o/e >= 0.6, length >= 200 bp); non-island CpG regions will be targeted."))
    }
    islands$gc_percent <- sprintf("%.1f", islands$gc_percent)
    islands$cpg_obs_exp <- sprintf("%.2f", islands$cpg_obs_exp)
    stats::setNames(islands, c("Island", "TSS-relative", "Genomic", "Length (bp)", "CpGs", "GC %", "CpG o/e", "Design window"))
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$target_preview <- shiny::renderTable({
    shiny::req(extracted())
    planned <- planned_targets()
    if (is.character(planned)) return(data.frame(note = planned))
    if (!nrow(planned)) {
      return(data.frame(note = sprintf("No CpG island detected; up to %s non-island CpG regions will be designed.", display_params()$fallback_targets)))
    }
    data.frame(target = planned$target_id, `TSS-relative` = format_tss_interval(planned$tss_start, planned$tss_end),
               `length (bp)` = planned$length, CpGs = planned$n_cpg, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE, spacing = "s", caption = "Targets to design", caption.placement = "top")

  output$sequence_header <- shiny::renderText({
    value <- extracted()
    if (is.null(value)) return(NULL)
    paste0(">", fasta_header(value))
  })
  output$sequence_text <- shiny::renderText({ if (!is.null(extracted())) extracted()$sequence })

  output$download_fasta <- shiny::downloadHandler(
    filename = function() "requested_promoter.fa",
    content = function(file) write_sequence_fasta(extracted(), "requested", file)
  )
  output$download_complementary_fasta <- shiny::downloadHandler(
    filename = function() "complementary_promoter.fa",
    content = function(file) write_sequence_fasta(extracted(), "complementary", file)
  )
  output$download_bisulfite_requested_fasta <- shiny::downloadHandler(
    filename = function() "bisulfite_requested_promoter.fa",
    content = function(file) write_sequence_fasta(extracted(), "bisulfite_requested", file)
  )
  output$download_bisulfite_complementary_fasta <- shiny::downloadHandler(
    filename = function() "bisulfite_complementary_promoter.fa",
    content = function(file) write_sequence_fasta(extracted(), "bisulfite_complementary", file)
  )

  # BSP design ---------------------------------------------------------------

  shiny::observeEvent(input$bsp_reset, {
    defaults <- bsp_primer_params()
    for (name in names(defaults)) shiny::updateNumericInput(session, bsp_param_input_id(name), value = defaults[[name]])
  })

  shiny::observeEvent(input$design_primers, {
    if (is.null(extracted())) {
      message_text("Extract a sequence before designing primers.")
      return()
    }
    message_text(run_design(extracted(), target_mode()))
  })

  selected_rank <- shiny::reactive({
    value <- design()
    rank <- suppressWarnings(as.integer(input$bsp_pair))
    shiny::req(value, value$pairs, length(rank) == 1L, !is.na(rank), rank %in% value$pairs$rank)
    rank
  })

  output$bsp_targets <- shiny::renderTable({
    value <- design()
    if (is.null(value) || !nrow(value$targets)) return(NULL)
    t <- value$targets
    data.frame(target = t$target_id, kind = gsub("_", "-", t$kind), `TSS-relative` = format_tss_interval(t$tss_start, t$tss_end),
               `length (bp)` = t$length, CpGs = t$n_cpg, pairs = t$n_pairs, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE, spacing = "s", caption = "Designed targets", caption.placement = "top")

  output$bsp_table <- shiny::renderTable({
    value <- design()
    if (is.null(value) || is.null(value$pairs)) return(NULL)
    columns <- c("rank", "target_id", "strand", "assay_strand", "left_sequence", "right_sequence", "product_size", "cpg_in_amplicon",
                 "left_tm", "right_tm", "amplicon_tss", "score", "specific", "n_offtarget_products")
    value$pairs[, intersect(columns, names(value$pairs)), drop = FALSE]
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$bsp_hints <- shiny::renderUI({
    value <- design()
    if (is.null(value)) return(NULL)
    hints <- design_hints(value)
    if (!length(hints)) return(NULL)
    shiny::div(class = "alert alert-warning small", shiny::strong("Why targets got no pair"),
               shiny::tags$ul(class = "mb-0", lapply(hints, shiny::tags$li)))
  })

  output$bsp_diagnostics <- shiny::renderTable({
    value <- design()
    if (is.null(value) || is.null(value$diagnostics)) return(NULL)
    empty <- value$targets$target_id[value$targets$n_pairs == 0L]
    rows <- value$diagnostics[value$diagnostics$target_id %in% empty, , drop = FALSE]
    if (!nrow(rows)) return(NULL)
    stats::setNames(rows, c("Target", "Strand", "Assay", "Left: bisulfite rules", "Right: bisulfite rules", "Left: Primer3", "Right: Primer3",
                            "Pairs in range", "Pairs checked", "Left: main rejections", "Right: main rejections"))
  }, striped = TRUE, bordered = TRUE, spacing = "s", caption = "Candidate primers remaining after each stage, for targets without pairs", caption.placement = "top")

  output$bsp_pair_detail <- shiny::renderTable({
    pair <- design_pairs(design(), selected_rank())
    side <- function(prefix) c(pair[[paste0(prefix, "_sequence")]], pair[[paste0(prefix, "_length")]], pair[[paste0(prefix, "_tm")]],
                               pair[[paste0(prefix, "_gc")]], pair[[paste0(prefix, "_converted_c")]], pair[[paste0(prefix, "_cpg")]],
                               pair[[paste0(prefix, "_interval")]], pair[[paste0(prefix, "_tss")]],
                               pair[[paste0(prefix, "_self_any")]], pair[[paste0(prefix, "_self_end")]], pair[[paste0(prefix, "_hairpin")]])
    stats::setNames(
      data.frame(c("Forward (left)", "Reverse (right)"), rbind(side("left"), side("right")), stringsAsFactors = FALSE),
      c("primer", "sequence 5'-3'", "length", "Tm", "GC %", "converted C", "CpG", "genomic", "TSS-relative",
        "self any (Tm)", "self end (Tm)", "hairpin (Tm)"))
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$bsp_sequence_view <- shiny::renderUI(sequence_view_html(design(), selected_rank()))

  output$download_bsp_pairs <- shiny::downloadHandler(
    filename = function() sprintf("%s_bsp_pairs.tsv", design_label(design())),
    content = function(file) write_bsp_pairs(design(), file)
  )
  output$download_bsp_pairs_csv <- shiny::downloadHandler(
    filename = function() sprintf("%s_bsp_pairs.csv", design_label(design())),
    content = function(file) write_bsp_pairs(design(), file, sep = ",")
  )
  output$download_order_sheet <- shiny::downloadHandler(
    filename = function() sprintf("%s_order.tsv", pair_names(design(), design_pairs(design(), selected_rank()))),
    content = function(file) write_order_sheet(design(), file, ranks = selected_rank())
  )
  output$download_amplicons <- shiny::downloadHandler(
    filename = function() sprintf("%s_amplicon.fa", pair_names(design(), design_pairs(design(), selected_rank()))),
    content = function(file) write_amplicon_references(design(), file, ranks = selected_rank())
  )

  shiny::observeEvent(input$check_specificity, {
    value <- design()
    if (is.null(value) || is.null(value$pairs)) {
      message_text("Design primer pairs before checking specificity.")
      return()
    }
    tryCatch({
      genome <- optional_path(input$bismark_genome)
      if (is.null(genome)) stop("Set the Bismark genome folder under 'External tools'.", call. = FALSE)
      bowtie2 <- optional_path(input$bowtie2_path)
      shiny::withProgress(message = "Aligning primers to the converted genome", value = 0.3, {
        checked <- check_bsp_specificity(value, genome, top_n = input$spec_top_n, max_mismatch = input$spec_max_mismatch,
                                         max_product = input$spec_max_product, bowtie2 = if (is.null(bowtie2)) bowtie2_path() else bowtie2)
        design(checked)
      })
      done <- checked$pairs[!is.na(checked$pairs$specific), , drop = FALSE]
      message_text(sprintf("Specificity checked for %s pair%s: %s specific.", nrow(done), if (nrow(done) == 1L) "" else "s", sum(done$specific)))
    }, error = function(e) message_text(conditionMessage(e)))
  })
}
