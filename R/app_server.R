app_server <- function(input, output, session) {
  annotation <- shiny::reactiveVal(NULL)
  tss_records <- shiny::reactiveVal(NULL)
  selected_tss <- shiny::reactiveVal(NULL)
  extracted <- shiny::reactiveVal(NULL)
  results <- shiny::reactiveVal(NULL)
  message_text <- shiny::reactiveVal(NULL)

  output$status <- shiny::renderUI({
    value <- message_text()
    if (is.null(value)) return(NULL)
    shiny::div(class = "alert alert-info", value)
  })

  shiny::observeEvent(input$find_gene, {
    shiny::req(input$gtf_path, input$gene_query)
    tryCatch({
      shiny::withProgress(message = "Finding gene and TSSs", value = 0, {
        shiny::incProgress(0.35, detail = "Loading GTF annotation")
        gtf <- load_annotation(trimws(input$gtf_path))
        shiny::incProgress(0.35, detail = "Calculating transcript TSSs")
        records <- collapse_tss(get_transcript_tss(gtf, gene = input$gene_query))
        shiny::incProgress(0.2, detail = "Preparing TSS choices")
        annotation(gtf)
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

  shiny::observeEvent(input$extract, {
    shiny::req(tss_records(), input$tss_choice, input$upstream, input$downstream)
    tryCatch({
      shiny::withProgress(message = "Extracting reference sequence", value = 0, {
        shiny::incProgress(0.2, detail = "Opening FASTA and index")
        record <- tss_records()[as.integer(input$tss_choice), , drop = FALSE]
        reference <- fasta_input()
        shiny::incProgress(0.2, detail = "Reading genomic interval")
        value <- extract_tss_region(record, reference, input$upstream, input$downstream)
        shiny::incProgress(0.5, detail = "Orienting sequence")
        selected_tss(record)
        extracted(value)
        results(NULL)
        shiny::incProgress(0.1)
        message_text(sprintf("Extracted %s bp from %s:%s-%s.", value$sequence_length, value$scaffold, value$region_start, value$region_end))
      })
    }, error = function(e) message_text(conditionMessage(e)))
  })

  primer_data <- shiny::reactive({
    if (identical(input$primer_mode, "paste")) return(parse_primer_output(input$primer_text))
    data.frame(
      primer_type = c("Left primer", "Right primer"),
      start = c(input$left_start, input$right_start), size = c(input$left_size, input$right_size),
      stringsAsFactors = FALSE
    )
  })

  shiny::observeEvent(input$convert, {
    shiny::req(extracted())
    tryCatch({
      shiny::withProgress(message = "Converting primer coordinates", value = 0, {
        shiny::incProgress(0.3, detail = "Parsing primer Start and Size")
        primers <- primer_data()
        shiny::incProgress(0.5, detail = "Mapping intervals by strand")
        converted <- map_primer_coordinates(primers, extracted(), input$start_convention)
        results(converted)
        shiny::incProgress(0.2)
        message_text("Primer coordinates converted from extracted-sequence positions.")
      })
    }, error = function(e) message_text(conditionMessage(e)))
  })

  output$sequence_header <- shiny::renderText({
    value <- extracted()
    if (is.null(value)) return(NULL)
    paste0(">", fasta_header(value))
  })
  output$sequence_text <- shiny::renderText({ if (!is.null(extracted())) extracted()$sequence })
  output$primer_table <- shiny::renderTable({
    value <- results()
    if (is.null(value)) return(NULL)
    value[, c("primer_type", "relative_start", "relative_end", "size", "scaffold", "genomic_start", "genomic_end", "binding_strand"), drop = FALSE]
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  output$amplicon_interval <- shiny::renderText({
    value <- results()
    if (is.null(value) || nrow(value) < 2L || length(unique(value$scaffold)) != 1L) return(NULL)
    sprintf("Amplicon span: %s:%s-%s (%s bp)", value$scaffold[1L], min(value$genomic_start), max(value$genomic_end), max(value$genomic_end) - min(value$genomic_start) + 1L)
  })
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
  output$download_results <- shiny::downloadHandler(
    filename = function() "primer_coordinates.tsv",
    content = function(file) write_primer_results(results(), file)
  )
  output$download_results_csv <- shiny::downloadHandler(
    filename = function() "primer_coordinates.csv",
    content = function(file) write_primer_results(results(), file, sep = ",")
  )
}
