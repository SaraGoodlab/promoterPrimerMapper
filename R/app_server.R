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
  message_type <- shiny::reactiveVal("info")
  extraction_snapshot <- shiny::reactiveVal(NULL)
  design_snapshot <- shiny::reactiveVal(NULL)
  design_failed <- shiny::reactiveVal(FALSE)
  gene_snapshot <- shiny::reactiveVal(NULL)
  selected_pair <- shiny::reactiveVal(NULL)
  selection_visible <- shiny::reactiveVal(TRUE)
  map_bounds <- shiny::reactiveVal(NULL)
  specificity_failed <- shiny::reactiveVal(FALSE)

  # Check metadata outside reactive observers: even a no-op reactivePoll timer
  # produces Shiny busy events. reactiveVal only invalidates on an actual change.
  reference_state <- shiny::reactiveVal(shiny::isolate(watched_file_state(input)))
  shiny::observe({ reference_state(watched_file_state(input)) })
  cancel_poll <- NULL
  schedule_poll <- function() {
    # Honor a session-provided clock (e.g. testServer's virtual clock) without
    # leaving a recurring timer in the loop used to resolve test output promises.
    if (is.function(session$.scheduleTask)) session$.scheduleTask(3000, poll_files)
    else later::later(poll_files, delay = 3)
  }
  poll_files <- function() {
    if (!is.null(cancel_poll)) { cancel_poll(); cancel_poll <<- NULL }
    if (session$isClosed()) return(invisible(NULL))
    shiny::withReactiveDomain(session, shiny::isolate(reference_state(watched_file_state(input))))
    cancel_poll <<- schedule_poll()
    invisible(NULL)
  }
  cancel_poll <- schedule_poll()
  session$onSessionEnded(function() { if (!is.null(cancel_poll)) cancel_poll() })
  last_controls <- shiny::reactiveVal(NULL)

  # Only explicit computations show a working state; background polling stays silent.
  with_operation_progress <- function(expr, message, value = 0, ...) {
    session$sendCustomMessage("ppm-operation", list(active = TRUE, label = message))
    on.exit(session$sendCustomMessage("ppm-operation", list(active = FALSE, label = "")), add = TRUE)
    shiny::withProgress(expr = expr, message = message, value = value, ...)
  }

  set_message <- function(text, type = "info") {
    message_type(type)
    message_text(text)
  }
  gene_inputs <- shiny::reactive({ reference_state(); list(annotation = reference_identity(input$gtf_path), gene = trimws(input_or(input$gene_query, ""))) })
  current_extraction <- shiny::reactive({ reference_state(); extraction_inputs(input, tss_records()) })
  extraction_changed <- shiny::reactive(!is.null(extraction_snapshot()) && !identical(current_extraction(), extraction_snapshot()))
  current_design <- shiny::reactive(design_inputs(input, target_mode(),
    tryCatch(validate_bsp_params(bsp_params()), error = function(e) bsp_params())))
  design_changed <- shiny::reactive(!is.null(design_snapshot()) && !identical(current_design(), design_snapshot()))

  output$has_extraction <- shiny::reactive(!is.null(extracted()))
  output$has_pairs <- shiny::reactive(has_primer_pairs(design()))
  shiny::outputOptions(output, "has_extraction", suspendWhenHidden = FALSE)
  shiny::outputOptions(output, "has_pairs", suspendWhenHidden = FALSE)

  reference_problems <- shiny::reactive({
    reference_state()
    Filter(Negate(is.null), list(file_problem(input$gtf_path, "Annotation"), file_problem(input$fasta_path, "FASTA"),
                                file_problem(input$fasta_index_path, "FASTA index")))
  })
  primer3_problem <- shiny::reactive({ reference_state(); tool_problem(input$primer3_path, primer3_path, "Primer3") })
  specificity_inputs <- shiny::reactive({
    reference_state()
    list(genome = reference_identity(input$bismark_genome), bowtie2 = trimws(input_or(input$bowtie2_path, "")),
         top_n = as.numeric(input_or(input$spec_top_n, 5)), max_mismatch = as.numeric(input_or(input$spec_max_mismatch, 2)),
         max_product = as.numeric(input_or(input$spec_max_product, 3000)))
  })
  specificity_problem <- shiny::reactive({
    reference_state()
    tryCatch({
      tool <- tool_problem(input$bowtie2_path, bowtie2_path, "Bowtie2")
      if (!is.null(tool)) return(tool)
      if (is.null(optional_path(input$bismark_genome))) return("Enter the Bismark genome folder to check specificity.")
      bismark_indexes(trimws(input$bismark_genome))
      p <- specificity_inputs()
      assert_scalar_integer(p$top_n, "Pairs per target", 1L)
      assert_scalar_integer(p$max_mismatch, "Max mismatches", 0L)
      assert_scalar_integer(p$max_product, "Max product", 50L)
      NULL
    }, error = function(e) conditionMessage(e))
  })
  output$reference_readiness <- shiny::renderUI({
    issues <- reference_problems()
    shiny::div(class = "ppm-help", role = "status", if (length(issues)) shiny::tags$ul(lapply(issues, shiny::tags$li)) else "Reference files are readable.")
  })
  output$primer3_readiness <- shiny::renderUI(shiny::p(class = "ppm-help", input_or(primer3_problem(), "Primer3 is ready.")))
  output$specificity_readiness <- shiny::renderUI(shiny::p(class = "ppm-help", input_or(specificity_problem(), "Bowtie2 and Bismark indexes are ready.")))
  shiny::observe({
    ready_gene <- is.null(file_problem(input$gtf_path, "Annotation")) && nzchar(trimws(input_or(input$gene_query, "")))
    ready_extract <- length(reference_problems()) == 0L && identical(gene_inputs(), gene_snapshot()) &&
      length(input$tss_choice) == 1L && nzchar(input$tss_choice)
    planned <- if (is.null(extracted())) NULL else planned_targets()
    controls <- list(find_gene = ready_gene, extract = ready_extract,
      use_selection = is.numeric(selection()),
      design_primers = !is.null(extracted()) && !extraction_changed() && is.null(parameter_error()) &&
        is.data.frame(planned) && is.null(primer3_problem()),
      check_specificity = has_primer_pairs(design()) && !extraction_changed() && !design_changed() && !design_failed() && is.null(specificity_problem()))
    if (!identical(controls, shiny::isolate(last_controls()))) {
      last_controls(controls)
      session$sendCustomMessage("ppm-controls", controls)
    }
  })
  shiny::observeEvent(input$fasta_path, {
    old <- shiny::isolate(input$fasta_index_path)
    if (is.null(old) || !nzchar(old) || identical(old, last_suggested_index())) {
      next_index <- paste0(trimws(input$fasta_path), ".fai")
      last_suggested_index(next_index)
      shiny::updateTextInput(session, "fasta_index_path", value = next_index)
    }
  }, ignoreInit = TRUE)
  last_suggested_index <- shiny::reactiveVal(NULL)

  output$design_context <- shiny::renderUI({
    value <- extracted()
    if (is.null(value)) return(empty_state("Start with a gene promoter", "Choose reference files and find a gene."))
    shiny::div(class = "ppm-context",
      shiny::strong(design_label(list(extracted = value))),
      sprintf(" | %s:%s-%s (%s strand) | TSS %s | %s bp", value$scaffold, value$region_start, value$region_end,
              value$strand, value$tss, value$sequence_length),
      if (!is.null(design_snapshot())) shiny::span(class = "ppm-badge", paste("Completed design:", design_snapshot()$mode)))
  })
  output$result_state <- shiny::renderUI({
    if (is.null(extracted())) return(NULL)
    if (extraction_changed()) return(shiny::div(class = "alert alert-warning", role = "status", "Reference, gene, TSS or extraction settings changed. The map and exports still belong to the previous extraction. Find the gene and extract again."))
    if (!is.null(design()) && (design_changed() || design_failed())) return(shiny::div(class = "alert alert-warning", role = "status",
      if (design_failed()) "The last design attempt failed. Previous results are retained." else "Design settings changed. Results and exports belong to the previous design.",
      " Choose Design BSP primers to update results."))
    NULL
  })

  output$status <- shiny::renderUI({
    value <- message_text()
    if (is.null(value)) return(NULL)
    shiny::div(class = paste("alert", paste0("alert-", message_type())), role = "status", `aria-live` = "polite", value)
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
      with_operation_progress(message = "Finding gene and TSSs", value = 0, {
        shiny::incProgress(0.35, detail = "Loading GTF annotation")
        gtf <- cached_annotation(trimws(input$gtf_path))
        shiny::incProgress(0.35, detail = "Calculating transcript TSSs")
        records <- collapse_tss(get_transcript_tss(gtf, gene = input$gene_query))
        shiny::incProgress(0.2, detail = "Preparing TSS choices")
        tss_records(records)
        gene_snapshot(gene_inputs())
        choices <- stats::setNames(seq_len(nrow(records)), sprintf("%s:%s (%s, %s transcripts)", records$scaffold, records$tss, records$strand, records$n_transcripts))
        shiny::updateSelectInput(session, "tss_choice", choices = choices, selected = if (nrow(records) == 1L) "1" else character(0))
        shiny::incProgress(0.1)
        set_message(sprintf("Found %s distinct TSS%s.", nrow(records), if (nrow(records) == 1L) "" else "s"), "success")
      })
    }, error = function(e) set_message(conditionMessage(e), "danger"))
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
  # Invalid input is shown explicitly rather than silently substituting defaults.
  parameter_error <- shiny::reactive(tryCatch({ validate_bsp_params(bsp_params()); NULL }, error = function(e) conditionMessage(e)))
  display_params <- shiny::reactive({ shiny::req(is.null(parameter_error())); validate_bsp_params(bsp_params()) })
  output$parameter_validation <- shiny::renderUI({
    error <- parameter_error()
    if (is.null(error)) return(NULL)
    shiny::div(class = "alert alert-danger", role = "alert", error)
  })
  output$design_readiness <- shiny::renderUI({
    if (extraction_changed()) return(shiny::p(class = "ppm-help", "Extract the updated promoter before designing again."))
    if (!is.null(primer3_problem())) return(shiny::p(class = "ppm-help", primer3_problem(), " Open Reference and gene, then Design tool settings."))
    if (!is.null(parameter_error())) return(shiny::p(class = "ppm-help", "Correct the highlighted parameters to enable design."))
    NULL
  })
  shiny::observeEvent(parameter_error(), {
    error <- parameter_error()
    if (is.null(error)) return()
    spec <- bsp_param_spec()
    groups <- unique(vapply(spec[vapply(spec, function(s) grepl(s$label, error, fixed = TRUE), logical(1))], `[[`, character(1), "group"))
    if (!length(groups)) groups <- if (grepl("Tm|Primer length", error)) "Primer size & Tm" else if (grepl("GC", error)) "Composition" else "Amplicon"
    bslib::accordion_panel_open("bsp_param_groups", values = groups, session = session)
    session$sendCustomMessage("ppm-validation", list())
  }, ignoreNULL = TRUE)
  target_mode <- shiny::reactive(if (is.null(input$bsp_target_mode)) "auto" else input$bsp_target_mode)

  # Design the current extraction; returns the status text.
  run_design <- function(value, mode) {
    tryCatch({
      params <- validate_bsp_params(bsp_params(), value$sequence_length)
      target <- if (identical(mode, "auto")) "auto" else mode_targets(value, mode, params, input$bsp_island, input$range_from, input$range_to)
      primer3 <- optional_path(input$primer3_path)
      strand <- if (is.null(input$bsp_strand)) "both" else input$bsp_strand
      result <- with_operation_progress(message = "Designing BSP primers", value = 0, {
        design_bsp_primers(value, strand, target, params, primer3 = if (is.null(primer3)) primer3_path() else primer3,
                           progress = function(step, steps, detail) shiny::setProgress((step - 1) / steps, detail = detail))
      })
      design_snapshot(design_inputs(input, mode, params))
      design_failed(FALSE)
      specificity_failed(FALSE)
      result$gui_provenance <- list(extraction = extraction_snapshot(), design = design_snapshot(), completed = as.character(Sys.time()))
      message_type(if (has_primer_pairs(result)) "success" else "warning")
      design(result)
      rank <- if (has_primer_pairs(result)) result$pairs$rank[1L] else NULL
      selected_pair(rank)
      shiny::updateSelectInput(session, "bsp_pair", choices = pair_choices(result$pairs), selected = rank)
      shiny::updateSelectInput(session, "pair_target_filter", choices = c("All targets" = "all", stats::setNames(result$targets$target_id, result$targets$target_id)), selected = "all")
      design_message(result)
    }, error = function(e) {
      design_failed(TRUE)
      message_type("danger")
      paste("Primer design failed:", conditionMessage(e))
    })
  }

  shiny::observeEvent(input$extract, {
    shiny::req(tss_records(), input$tss_choice, input$upstream, input$downstream)
    if (!identical(gene_inputs(), gene_snapshot())) {
      set_message("Reference or gene changed. Choose Find gene before extracting again.", "warning")
      return()
    }
    value <- tryCatch({
      with_operation_progress(message = "Extracting reference sequence", value = 0, {
        shiny::incProgress(0.2, detail = "Opening FASTA and index")
        record <- tss_records()[as.integer(input$tss_choice), , drop = FALSE]
        reference <- fasta_input()
        shiny::incProgress(0.3, detail = "Reading genomic interval")
        value <- extract_tss_region(record, reference, input$upstream, input$downstream)
        selected_tss(record)
        value
      })
    }, error = function(e) {
      set_message(conditionMessage(e), "danger")
      NULL
    })
    if (is.null(value)) return()
    extracted(value)
    extraction_snapshot(current_extraction())
    map_bounds(NULL)
    selection_visible(TRUE)
    design(NULL)
    design_snapshot(NULL)
    design_failed(FALSE)
    shiny::updateSelectInput(session, "bsp_pair", choices = character(0))
    islands <- cpg_island_scan(value)
    nearest <- nearest_island(cpg_island_scan(value), tss_sense_position(value))
    shiny::updateSelectInput(session, "bsp_island", choices = stats::setNames(as.character(seq_len(nrow(islands))),
                             sprintf("CGI%s: %s (%s bp)", seq_len(nrow(islands)), format_tss_interval(islands$tss_start, islands$tss_end), islands$length)),
                             selected = if (is.na(nearest)) character(0) else as.character(nearest))
    shiny::updateRadioButtons(session, "bsp_target_mode", selected = "auto")
    extracted_text <- sprintf("Extracted %s bp from %s:%s-%s.", value$sequence_length, value$scaffold, value$region_start, value$region_end)
    message_text(paste(extracted_text, run_design(value, "auto")))
    bslib::toggle_sidebar("setup_sidebar", open = FALSE, session = session)
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
    if (!is.null(parameter_error())) return(parameter_error())
    tryCatch(mode_targets(value, target_mode(), display_params(), input$bsp_island, input$range_from, input$range_to),
             error = function(e) conditionMessage(e))
  })

  shiny::observeEvent(input$cpg_brush, {
    value <- extracted()
    brush <- input$cpg_brush
    if (is.null(value) || is.null(brush$xmin) || is.null(brush$xmax)) return()
    selection_visible(TRUE)
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
    set_message(sprintf("Selected CpG island CGI%s as the target; press 'Design BSP primers' to design it.", k))
  })

  shiny::observeEvent(input$use_selection, {
    interval <- selection()
    if (!is.numeric(interval)) {
      set_message(if (is.character(interval)) interval else "Extract a sequence and select a range first.", "warning")
      return()
    }
    shiny::updateRadioButtons(session, "bsp_target_mode", selected = "range")
    selection_visible(TRUE)
    tss <- tss_sense_position(extracted())
    set_message(sprintf("Target set to the selected range %s; press 'Design BSP primers' to design it.",
                         format_tss_interval(tss_relative_position(interval[1L], tss), tss_relative_position(interval[2L], tss))))
  })

  shiny::observeEvent(list(input$range_from, input$range_to), selection_visible(TRUE), ignoreInit = TRUE)
  shiny::observeEvent(input$clear_selection, {
    selection_visible(FALSE)
    shiny::updateRadioButtons(session, "bsp_target_mode", selected = "auto")
    session$sendCustomMessage("ppm-clear-brush", list())
    set_message("Selection cleared. Target mode is All CpG islands; previous pairs keep their completed-design context.")
  })
  shiny::observeEvent(input$zoom_reset, map_bounds(NULL))
  shiny::observeEvent(input$zoom_target, {
    targets <- planned_targets()
    if (is.data.frame(targets) && nrow(targets)) map_bounds(c(min(targets$start), max(targets$end)))
  })
  shiny::observeEvent(input$zoom_pair, {
    if (!has_primer_pairs(design())) return()
    pair <- design_pairs(design(), selected_rank())
    template <- design()$templates[[pair$strand]]
    span <- template_to_sense(template, pair$template_left_start, pair$template_right_end)
    map_bounds(c(span$start, span$end))
  })

  output$map_context <- shiny::renderUI({
    if (is.null(extracted())) return(NULL)
    shiny::p(class = "ppm-help", if (is.null(design())) "Map shading: proposed targets." else
      "Map shading and pair: completed design. Proposed targets are listed below; rerun design to update shading.")
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
    rank <- selected_pair()
    targets <- NULL
    if (is.null(current)) {
      planned <- planned_targets()
      if (is.data.frame(planned)) targets <- planned
    }
    interval <- selection()
    island <- if (identical(target_mode(), "island")) suppressWarnings(as.integer(input$bsp_island)) else NULL
    plot_cpg_map(value, current, if (length(rank) == 1L && !is.na(rank)) rank else NULL,
                 targets = targets, selection = if (selection_visible() && is.numeric(interval)) interval else NULL, island = island,
                 view = map_bounds())
  })

  output$cursor_position <- shiny::renderText({
    value <- extracted()
    hover <- input$cpg_hover
    if (is.null(value)) return("")
    instruction <- "drag to select a range; click an island to target it. Target controls also work without the map."
    if (is.null(hover$x)) return(paste("Hover for coordinates;", instruction))
    paste(describe_position(value, hover$x), "|", instruction)
  })

  output$selection_table <- shiny::renderUI({
    interval <- selection()
    if (is.null(interval)) return(NULL)
    if (!selection_visible()) return(shiny::p(class = "ppm-help", "No range highlighted. Type a range or drag on the map."))
    if (is.character(interval)) return(shiny::div(class = "alert alert-danger", role = "alert", interval))
    region <- describe_region(extracted(), interval[1L], interval[2L])
    region$gc_percent <- sprintf("%.1f", region$gc_percent)
    shiny::tagList(metric_list(list("TSS-relative" = region$tss_relative, "Length (bp)" = region$length, "CpGs" = region$cpg, "GC (%)" = region$gc_percent)),
      shiny::p(class = "ppm-help", paste(region$genomic, "| Sequence", region$sequence, "| Islands:", region$cpg_islands)))
  })

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
  output$target_validation <- shiny::renderUI({
    if (is.null(extracted())) return(NULL)
    planned <- planned_targets()
    if (is.character(planned)) shiny::div(class = "alert alert-danger", role = "alert", planned) else NULL
  })

  output$sequence_header <- shiny::renderText({
    value <- extracted()
    if (is.null(value)) return(NULL)
    paste0(">", fasta_header(value))
  })
  output$sequence_text <- shiny::renderText({ if (!is.null(extracted())) extracted()$sequence })

  output$download_fasta <- shiny::downloadHandler(
    filename = function() "requested_promoter.fa",
    content = function(file) { shiny::req(extracted()); write_sequence_fasta(extracted(), "requested", file) }
  )
  output$download_complementary_fasta <- shiny::downloadHandler(
    filename = function() "complementary_promoter.fa",
    content = function(file) { shiny::req(extracted()); write_sequence_fasta(extracted(), "complementary", file) }
  )
  output$download_bisulfite_requested_fasta <- shiny::downloadHandler(
    filename = function() "bisulfite_requested_promoter.fa",
    content = function(file) { shiny::req(extracted()); write_sequence_fasta(extracted(), "bisulfite_requested", file) }
  )
  output$download_bisulfite_complementary_fasta <- shiny::downloadHandler(
    filename = function() "bisulfite_complementary_promoter.fa",
    content = function(file) { shiny::req(extracted()); write_sequence_fasta(extracted(), "bisulfite_complementary", file) }
  )

  # BSP design ---------------------------------------------------------------

  shiny::observeEvent(input$bsp_reset, {
    defaults <- bsp_primer_params()
    for (name in names(defaults)) shiny::updateNumericInput(session, bsp_param_input_id(name), value = defaults[[name]])
  })

  shiny::observeEvent(input$design_primers, {
    if (is.null(extracted())) {
      set_message("Extract a sequence before designing primers.", "warning")
      return()
    }
    if (extraction_changed()) {
      set_message("Extraction settings changed. Extract the promoter again before designing primers.", "warning")
      return()
    }
    message_text(run_design(extracted(), target_mode()))
  })

  selected_rank <- shiny::reactive({
    value <- design()
    rank <- selected_pair()
    shiny::req(value, value$pairs, length(rank) == 1L, !is.na(rank), rank %in% value$pairs$rank)
    rank
  })
  choose_pair <- function(rank) {
    rank <- suppressWarnings(as.integer(rank))
    if (!has_primer_pairs(design()) || length(rank) != 1L || is.na(rank) || !rank %in% design()$pairs$rank) return()
    selected_pair(rank)
  }
  shiny::observeEvent(input$bsp_pair, choose_pair(input$bsp_pair))
  shiny::observeEvent(input$pair_row, {
    choose_pair(input$pair_row)
    shiny::updateSelectInput(session, "bsp_pair", selected = selected_pair())
  })
  output$pair_state <- shiny::renderUI({
    value <- design()
    if (is.null(value)) return(empty_state("No design yet", "Choose a target and Design BSP primers."))
    if (!has_primer_pairs(value)) return(empty_state("No pairs passed the constraints", design_message(value)))
    NULL
  })

  output$bsp_targets <- shiny::renderTable({
    value <- design()
    if (is.null(value) || !nrow(value$targets)) return(NULL)
    t <- value$targets
    data.frame(target = t$target_id, kind = gsub("_", "-", t$kind), `TSS-relative` = format_tss_interval(t$tss_start, t$tss_end),
               `length (bp)` = t$length, CpGs = t$n_cpg, pairs = t$n_pairs, check.names = FALSE)
  }, striped = TRUE, bordered = TRUE, spacing = "s", caption = "Designed targets", caption.placement = "top")

  output$bsp_table <- shiny::renderUI({
    value <- design()
    if (!has_primer_pairs(value)) return(NULL)
    pairs <- filter_pairs(value$pairs, input_or(input$pair_target_filter, "all"), input_or(input$pair_strand_filter, "all"), input_or(input$pair_sort, "rank"))
    pair_table_ui(pairs, selected_pair())
  })

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

  output$bsp_pair_detail <- shiny::renderUI({
    if (!has_primer_pairs(design())) return(empty_state("No pair selected", "Design primers to review oligos and their amplicon."))
    rank <- selected_rank()
    pair <- design_pairs(design(), rank)
    pair_detail_ui(pair)
  })

  output$bsp_sequence_view <- shiny::renderUI({
    if (!has_primer_pairs(design())) return(empty_state("No pair selected", "Design primers to review oligos and their amplicon."))
    rank <- selected_rank()
    sequence_view_html(design(), rank, full = isTRUE(input$sequence_full))
  })

  output$export_context <- shiny::renderUI({
    value <- extracted()
    if (is.null(value)) return(NULL)
    text <- sprintf("Promoter exports: %s, %s:%s-%s (%s), TSS %s.", design_label(list(extracted = value)),
                    value$scaffold, value$region_start, value$region_end, value$strand, value$tss)
    if (has_primer_pairs(design())) text <- paste(text, sprintf("Selected-pair exports: %s.", pair_names(design(), design_pairs(design(), selected_rank()))))
    shiny::p(class = "ppm-help", role = "status", text,
      if (extraction_changed() || design_changed() || design_failed()) " Exports retain the previous completed result, not the edited settings.")
  })
  output$specificity_state <- shiny::renderUI({
    value <- design()
    if (!has_primer_pairs(value)) return(empty_state("Design pairs before checking specificity"))
    checked <- value$specificity
    if (is.null(checked)) return(shiny::p("Not checked. Specificity is optional and requires a prepared genome."))
    settings <- checked$settings
    count <- sum(!is.na(value$pairs$specific))
    saved <- value$gui_specificity
    changed <- !is.null(saved) && !identical(saved, specificity_inputs())
    shiny::div(shiny::p(sprintf("Completed check: %s pairs; up to %s per target, %s mismatches per hit, max product %s bp.",
      count, settings$top_n, settings$max_mismatch, settings$max_product)),
      if (specificity_failed()) shiny::p(class = "alert alert-warning", "The last check failed. Previous check results are retained."),
      if (changed) shiny::p(class = "alert alert-warning", "Specificity settings changed. Displayed results use the previous criteria; check again to update them."))
  })

  output$download_bsp_pairs_csv <- shiny::downloadHandler(
    filename = function() { shiny::req(has_primer_pairs(design())); sprintf("%s_primers.csv", design_label(design())) },
    contentType = "text/csv",
    content = function(file) { shiny::req(has_primer_pairs(design())); write_order_sheet(design(), file) }
  )
  output$download_order_sheet <- shiny::downloadHandler(
    filename = function() { rank <- selected_rank(); sprintf("%s_primers.csv", pair_names(design(), design_pairs(design(), rank))) },
    contentType = "text/csv",
    content = function(file) { rank <- selected_rank(); write_order_sheet(design(), file, ranks = rank) }
  )
  output$download_amplicons <- shiny::downloadHandler(
    filename = function() { rank <- selected_rank(); sprintf("%s_amplicon.fa", pair_names(design(), design_pairs(design(), rank))) },
    content = function(file) { rank <- selected_rank(); write_amplicon_references(design(), file, ranks = rank) }
  )

  shiny::observeEvent(input$check_specificity, {
    value <- design()
    if (is.null(value) || is.null(value$pairs)) {
      set_message("Design primer pairs before checking specificity.", "warning")
      return()
    }
    if (extraction_changed() || design_changed() || design_failed()) {
      set_message("Update the extraction and design before checking specificity.", "warning")
      return()
    }
    tryCatch({
      genome <- optional_path(input$bismark_genome)
      if (is.null(genome)) stop("Set the Bismark genome folder beside the specificity controls.", call. = FALSE)
      bowtie2 <- optional_path(input$bowtie2_path)
      with_operation_progress(message = "Aligning primers to the converted genome", value = 0.3, {
        checked <- check_bsp_specificity(value, genome, top_n = input$spec_top_n, max_mismatch = input$spec_max_mismatch,
                                         max_product = input$spec_max_product, bowtie2 = if (is.null(bowtie2)) bowtie2_path() else bowtie2)
        checked$gui_specificity <- specificity_inputs()
        design(checked)
        specificity_failed(FALSE)
      })
      done <- checked$pairs[!is.na(checked$pairs$specific), , drop = FALSE]
      set_message(sprintf("Specificity checked for %s pair%s: %s specific within checked criteria.", nrow(done), if (nrow(done) == 1L) "" else "s", sum(done$specific)), "success")
    }, error = function(e) { specificity_failed(TRUE); set_message(conditionMessage(e), "danger") })
  })
}
