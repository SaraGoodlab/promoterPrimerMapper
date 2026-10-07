app_empty_design <- function(extracted, strand, target, params, ...) {
  targets <- bsp_targets(extracted, target, params)
  targets$n_pairs <- rep(0L, nrow(targets))
  structure(list(extracted = extracted, targets = targets, params = params,
                 pairs = NULL, diagnostics = NULL, fallback = FALSE), class = "bsp_design")
}

app_reference_inputs <- function() {
  fasta <- fixture_reference("fa")
  list(gtf_path = fixture_reference("gtf"),
       fasta_path = fasta, fasta_index_path = paste0(fasta, ".fai"), gene_query = "BSPPLUS",
       upstream = 1000, downstream = 500, tss_choice = "1")
}

# Deliberately non-consecutive ranks expose row-index/identity mixups in the GUI.
app_pair_design <- function(extracted, strand, target, params, ...) {
  design <- app_empty_design(extracted, strand, target, params)
  design$templates <- list(sense = bisulfite_template(extracted, "sense"), antisense = bisulfite_template(extracted, "antisense"))
  rows <- lapply(seq_len(3L), function(i) {
    template <- design$templates[[if (i == 2L) "antisense" else "sense"]]
    raw <- data.frame(left_start = 300L + i, left_end = 321L + i, right_start = 630L + i, right_end = 651L + i,
      left_sequence = substr(template$converted, 300L + i, 321L + i),
      right_sequence = reverse_complement_text(substr(template$converted, 630L + i, 651L + i)),
      product_size = 352L, product_tm = 70 + i, cpg_in_amplicon = i + 3L, left_tm = 56, right_tm = 57, tm_diff = 1,
      left_gc_percent = 30, right_gc_percent = 30, left_length = 22L, right_length = 22L,
      left_converted_c = 4L, right_converted_c = 4L, left_cpg_in_primer = 0L, right_cpg_in_primer = 0L,
      pair_penalty = 4 - i, left_penalty = 1, right_penalty = 1, compl_any = 0, compl_end = 0,
      left_self_any = 0, left_self_end = 0, left_hairpin = 0, right_self_any = 0, right_self_end = 0, right_hairpin = 0)
    annotate_pairs(extracted, template, raw)
  })
  design$pairs <- cbind(rank = c(1L, 3L, 7L), target_id = c("CGI1", "CGI2", "CGI1"), do.call(rbind, rows))
  design$targets$n_pairs <- as.integer(table(factor(design$pairs$target_id, levels = design$targets$target_id)))
  design
}

test_that("initial pair outputs are guarded and explain the next step", {
  shiny::testServer(app_server, {
    expect_match(as.character(output$bsp_pair_detail$html), "Design primers to review")
    expect_match(as.character(output$bsp_sequence_view$html), "Design primers to review")
    expect_match(as.character(output$design_context$html), "Choose reference files and find a gene.", fixed = TRUE)
  })
})

test_that("GUI primer downloads use the CSV schema and correct pair scope", {
  shiny::testServer(app_server, {
    extracted(fixture_promoter("BSPMINUS"))
    design(app_pair_design(extracted(), "both", "auto", bsp_primer_params()))
    selected_pair(3L)
    session$setInputs(pair_target_filter = "CGI1", pair_strand_filter = "sense")
    all <- utils::read.csv(output$download_bsp_pairs_csv, comment.char = "#", check.names = FALSE)
    selected <- utils::read.csv(output$download_order_sheet, comment.char = "#", check.names = FALSE)
    expect_equal(nrow(all), 6L)
    expect_equal(nrow(selected), 2L)
    expect_equal(ncol(all), 14L)
    expect_identical(names(selected), names(all))
    expect_equal(selected$Primer, paste0("BSPMINUS_CGI2_BSP3_", c("F", "R")))
    expect_equal(selected[["Product Tm"]], c(72, 72))
    expect_true(startsWith(readLines(output$download_bsp_pairs_csv)[1L], "# "))
    expect_true(startsWith(readLines(output$download_order_sheet)[1L], "# "))
  })
})

test_that("invalid parameters are reported without a default target preview", {
  shiny::testServer(app_server, {
    extracted(fixture_promoter("BSPPLUS"))
    session$setInputs(bsp_tm_min = 70)
    expect_match(parameter_error(), "Tm values must satisfy")
    expect_match(output$target_preview, "Tm values must satisfy")
    expect_match(as.character(output$parameter_validation$html), "alert-danger")
  })
})

test_that("completed design snapshots track changes and failed reruns", {
  local_mocked_bindings(design_bsp_primers = app_empty_design, primer3_path = function() "unused")
  shiny::testServer(app_server, {
    do.call(session$setInputs, app_reference_inputs())
    session$setInputs(find_gene = 1, extract = 1)
    expect_false(extraction_changed())
    expect_false(design_changed())
    previous <- design()
    session$setInputs(bsp_tm_opt = 60)
    expect_true(design_changed())
    expect_match(as.character(output$result_state$html), "previous design")
    session$setInputs(bsp_tm_min = 70, design_primers = 1)
    expect_identical(design(), previous)
    expect_true(design_failed())
    expect_match(as.character(output$result_state$html), "attempt failed")
    expect_equal(message_type(), "danger")
    session$setInputs(bsp_tm_min = 52, bsp_tm_opt = 58, design_primers = 2)
    expect_false(design_changed())
    expect_false(design_failed())
  })
})

test_that("changing gene or reference cannot silently reuse an old extraction", {
  local_mocked_bindings(design_bsp_primers = app_empty_design, primer3_path = function() "unused")
  shiny::testServer(app_server, {
    do.call(session$setInputs, app_reference_inputs())
    session$setInputs(find_gene = 1, extract = 1)
    previous <- extracted()
    session$setInputs(gene_query = "BSPMINUS", design_primers = 1)
    expect_true(extraction_changed())
    expect_identical(extracted(), previous)
    expect_match(message_text(), "Extract the promoter again")
    session$setInputs(extract = 2)
    expect_match(message_text(), "Find gene")
    session$setInputs(find_gene = 2, extract = 3)
    expect_equal(extracted()$gene_id, "BSPMINUS")
    expect_false(extraction_changed())
    session$setInputs(fasta_path = "/missing/reference.fa")
    expect_true(extraction_changed())
  })
})

test_that("sorting and filtering preserve stable pair identity across every view", {
  local_mocked_bindings(design_bsp_primers = app_pair_design, primer3_path = function() "unused")
  shiny::testServer(app_server, {
    do.call(session$setInputs, app_reference_inputs())
    session$setInputs(find_gene = 1, extract = 1)
    session$setInputs(pair_sort = "score", pair_row = 7)
    expect_equal(selected_rank(), 7L)
    table <- as.character(output$bsp_table$html)
    expect_true(regexpr('data-pair-rank="7"', table, fixed = TRUE) < regexpr('data-pair-rank="1"', table, fixed = TRUE))
    expect_match(as.character(output$bsp_pair_detail$html), "Pair 7")
    expect_match(as.character(output$export_context$html), "CGI1_BSP7")
    session$setInputs(pair_target_filter = "CGI2")
    expect_equal(selected_rank(), 7L)
    expect_match(as.character(output$bsp_table$html), "outside the current filters")
    session$setInputs(pair_row = 3)
    expect_equal(selected_rank(), 3L)
    expect_match(as.character(output$bsp_pair_detail$html), "antisense")
    expect_match(as.character(output$export_context$html), "CGI2_BSP3")
    session$setInputs(pair_row = 2)
    expect_equal(selected_rank(), 3L)
    session$setInputs(bsp_pair = "1")
    expect_equal(selected_rank(), 1L)
  })
})

test_that("specificity reruns retain selection and distinguish previous criteria and failures", {
  check <- function(design, bismark_genome, top_n, max_mismatch, max_product, ...) {
    design$pairs$specific <- c(TRUE, FALSE, NA)
    design$pairs$on_target <- c(TRUE, TRUE, NA)
    design$pairs$n_offtarget_products <- c(0L, 2L, NA)
    design$specificity <- list(settings = list(top_n = top_n, max_mismatch = max_mismatch, max_product = max_product))
    design
  }
  local_mocked_bindings(design_bsp_primers = app_pair_design, primer3_path = function() "unused", check_bsp_specificity = check)
  shiny::testServer(app_server, {
    do.call(session$setInputs, app_reference_inputs())
    session$setInputs(find_gene = 1, extract = 1, pair_row = 7,
      bismark_genome = tempdir(), spec_top_n = 1, spec_max_mismatch = 2, spec_max_product = 3000)
    session$setInputs(check_specificity = 1)
    expect_equal(selected_rank(), 7L)
    expect_equal(specificity_labels(design()$pairs), c("Specific within checked criteria", "Off-target products found", "Not checked"))
    expect_match(as.character(output$specificity_state$html), "2 pairs")
    session$setInputs(spec_max_mismatch = 1)
    expect_match(as.character(output$specificity_state$html), "previous criteria")
    session$setInputs(bismark_genome = "", check_specificity = 2)
    expect_equal(selected_rank(), 7L)
    expect_match(as.character(output$specificity_state$html), "last check failed")
  })
})

test_that("sequence cropping and map zoom preserve strand coordinates", {
  for (gene in c("BSPPLUS", "BSPMINUS")) {
    design <- app_pair_design(fixture_promoter(gene), "both", "auto", bsp_primer_params())
    for (rank in c(1L, 3L)) {
      cropped <- as.character(sequence_view_html(design, rank, full = FALSE))
      full <- as.character(sequence_view_html(design, rank, full = TRUE))
      expect_lt(nchar(cropped), nchar(full))
      expect_match(cropped, 'class="primer"', fixed = TRUE)
      expect_match(cropped, design$pairs$assay_strand[match(rank, design$pairs$rank)], fixed = TRUE)
    }
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot_cpg_map(design$extracted, view = c(300L, 650L)))
  expect_equal(graphics::par("usr")[1:2], c(250, 700) + c(-1, 1) * .04 * 450)
})

test_that("exports record completed sources even after inputs change", {
  local_mocked_bindings(design_bsp_primers = app_pair_design, primer3_path = function() "unused")
  shiny::testServer(app_server, {
    reference <- app_reference_inputs()
    do.call(session$setInputs, reference)
    session$setInputs(find_gene = 1, extract = 1)
    session$setInputs(fasta_path = "/another/reference.fa", bsp_tm_opt = 60)
    comments <- bsp_param_comments(design())
    expect_true(paste("reference:", normalizePath(reference$fasta_path)) %in% comments)
    expect_false(any(grepl("/another/reference", comments, fixed = TRUE)))
    expect_match(comments[3L], "tm_opt=58", fixed = TRUE)
    expect_match(as.character(output$export_context$html), "previous completed result")
  })
})

test_that("FASTA index suggestions never replace a custom index", {
  session <- shiny::MockShinySession$new()
  sent <- list()
  session$sendInputMessage <- function(id, message) sent[[id]] <<- message
  shiny::testServer(app_server, session = session, {
    session$setInputs(fasta_path = "", fasta_index_path = "")
    session$setInputs(fasta_path = "/first.fa")
    expect_equal(sent$fasta_index_path$value, "/first.fa.fai")
    session$setInputs(fasta_index_path = "/custom.fai", fasta_path = "/second.fa")
    expect_equal(sent$fasta_index_path$value, "/first.fa.fai")
  })
})

test_that("clearing selection and zooming never change the completed design", {
  local_mocked_bindings(design_bsp_primers = app_pair_design, primer3_path = function() "unused")
  shiny::testServer(app_server, {
    do.call(session$setInputs, app_reference_inputs())
    session$setInputs(find_gene = 1, extract = 1)
    completed <- design()
    session$setInputs(clear_selection = 1)
    expect_false(selection_visible())
    expect_match(as.character(output$selection_table$html), "No range highlighted")
    session$setInputs(zoom_pair = 1)
    expect_length(map_bounds(), 2L)
    session$setInputs(zoom_reset = 1)
    expect_null(map_bounds())
    expect_identical(design(), completed)
  })
})

test_that("unchanged polling leaves outputs and control messages alone", {
  session <- shiny::MockShinySession$new()
  messages <- list()
  session$sendCustomMessage <- function(type, message) messages[[length(messages) + 1L]] <<- list(type = type, value = message)
  shiny::testServer(app_server, session = session, {
    do.call(session$setInputs, app_reference_inputs())
    updates <- 0L
    consumer <- shiny::observe({
      reference_state()
      current_extraction()
      reference_problems()
      specificity_problem()
      updates <<- updates + 1L
    })
    session$flushReact()
    initial_updates <- updates
    initial_messages <- length(messages)
    for (i in seq_len(4L)) session$elapse(3100)
    expect_equal(updates, initial_updates)
    expect_length(messages, initial_messages)
    expect_false(any(vapply(messages, function(m) identical(m$type, "ppm-operation"), logical(1))))
    consumer$destroy()
  })
})

test_that("on-disk reference changes are detected without rerunning a design", {
  designs <- 0L
  stub <- function(...) { designs <<- designs + 1L; app_pair_design(...) }
  local_mocked_bindings(design_bsp_primers = stub, primer3_path = function() "unused")
  reference <- app_reference_inputs()
  fasta <- tempfile(fileext = ".fa")
  file.copy(reference$fasta_path, fasta)
  file.copy(reference$fasta_index_path, paste0(fasta, ".fai"))
  reference$fasta_path <- fasta
  reference$fasta_index_path <- paste0(fasta, ".fai")
  shiny::testServer(app_server, {
    do.call(session$setInputs, reference)
    session$setInputs(find_gene = 1, extract = 1)
    expect_false(extraction_changed())
    previous <- design()
    write("", file = reference$fasta_path, append = TRUE)
    session$elapse(3100)
    expect_true(extraction_changed())
    expect_identical(design(), previous)
    expect_equal(designs, 1L)
    expect_match(as.character(output$result_state$html), "previous extraction")
    file.remove(reference$fasta_index_path)
    session$elapse(3100)
    expect_match(paste(unlist(reference_problems()), collapse = " "), "FASTA index file is missing")
    expect_false(last_controls()$extract)
  })
})

test_that("index changes are detected even when the genome folder is unchanged", {
  folder <- tempfile("bismark-poll")
  ct <- file.path(folder, "Bisulfite_Genome", "CT_conversion", "BS_CT.1.bt2")
  ga <- file.path(folder, "Bisulfite_Genome", "GA_conversion", "BS_GA.1.bt2")
  dir.create(dirname(ct), recursive = TRUE)
  dir.create(dirname(ga), recursive = TRUE)
  writeLines("fixture", ct)
  writeLines("fixture", ga)
  local_mocked_bindings(bowtie2_path = function() unname(Sys.which("Rscript")))
  shiny::testServer(app_server, {
    session$setInputs(bismark_genome = folder)
    expect_null(specificity_problem())
    folder_before <- reference_identity(folder)
    file.remove(ct)
    expect_identical(reference_identity(folder), folder_before)
    session$elapse(3100)
    expect_match(specificity_problem(), "No Bowtie2 bisulfite indexes")
  })
})

test_that("operation messages stop on both success and computation errors", {
  session <- shiny::MockShinySession$new()
  messages <- list()
  session$sendCustomMessage <- function(type, message) {
    if (type == "ppm-operation") messages[[length(messages) + 1L]] <<- message
  }
  local_mocked_bindings(design_bsp_primers = app_empty_design, primer3_path = function() "unused")
  shiny::testServer(app_server, session = session, {
    do.call(session$setInputs, app_reference_inputs())
    session$setInputs(find_gene = 1, extract = 1)
    expect_true(any(vapply(messages, function(m) isTRUE(m$active), logical(1))))
    expect_false(tail(messages, 1L)[[1L]]$active)
    session$setInputs(primer3_path = "/missing/program", design_primers = 1)
    # Simulate a tool failure inside the computation to exercise on.exit cleanup.
    expect_error(with_operation_progress(message = "Testing failure", { stop("tool failed") }), "tool failed")
    expect_false(tail(messages, 1L)[[1L]]$active)
  })
})

test_that("tool permission changes update readiness without touching file contents", {
  skip_on_os("windows")
  program <- tempfile("primer3-poll")
  writeLines(c("#!/bin/sh", "exit 0"), program)
  Sys.chmod(program, "0755")
  shiny::testServer(app_server, {
    session$setInputs(primer3_path = program)
    expect_null(primer3_problem())
    before <- reference_identity(program)
    Sys.chmod(program, "0600")
    expect_identical(reference_identity(program), before)
    session$elapse(3100)
    expect_match(primer3_problem(), "not executable")
  })
})

test_that("file polling stops when a session closes", {
  calls <- 0L
  original <- watched_file_state
  local_mocked_bindings(watched_file_state = function(input) { calls <<- calls + 1L; original(input) })
  shiny::testServer(app_server, {
    do.call(session$setInputs, app_reference_inputs())
    session$elapse(3100)
    session$close()
    before <- calls
    session$elapse(10000)
    expect_equal(calls, before)
  })
})
