bsp_param_input_id <- function(name) paste0("bsp_", name)

bsp_params_ui <- function(params = bsp_primer_params()) {
  spec <- bsp_param_spec()
  groups <- unique(vapply(spec, `[[`, character(1), "group"))
  panels <- lapply(groups, function(group) {
    names <- names(spec)[vapply(spec, function(s) identical(s$group, group), logical(1))]
    inputs <- lapply(names, function(name) {
      s <- spec[[name]]
      help_id <- paste0(bsp_param_input_id(name), "_help")
      control <- shiny::numericInput(bsp_param_input_id(name), s$label, value = params[[name]], min = s$min, max = s$max, step = s$step)
      if (!is.null(s$help)) control <- htmltools::tagQuery(control)$find("input")$addAttrs(`aria-describedby` = help_id)$allTags()
      shiny::div(control, if (!is.null(s$help)) shiny::tags$small(id = help_id, class = "ppm-help", s$help))
    })
    bslib::accordion_panel(group, shiny::div(class = "ppm-parameter-grid", inputs), value = group)
  })
  shiny::tagList(
    do.call(bslib::accordion, c(list(id = "bsp_param_groups", open = FALSE), panels)),
    shiny::actionButton("bsp_reset", "Reset parameters", class = "btn-outline-secondary btn-sm mt-2")
  )
}

app_action <- function(id, label, class = "btn-outline-secondary", ...) {
  shiny::actionButton(id, label, class = paste(class, "ppm-action"), disabled = "disabled", ...)
}

scroll_table <- function(id, label) {
  shiny::div(class = "ppm-table-scroll", role = "region", `aria-label` = label, tabindex = "0", shiny::tableOutput(id))
}

bsp_view_css <- function() {
  colors <- app_colors()
  # Keep overrides after bslib's component styles, which load after Bootstrap.
  htmltools::attachDependencies(shiny::tagList(
    shiny::tags$style(shiny::HTML(paste0(":root {",
      paste(sprintf("--ppm-%s: %s;", names(colors), colors), collapse = " "), "}"))),
    shiny::includeCSS(system.file("app", "gui.css", package = "promoterPrimerMapper"))), app_assets())
}

app_ui <- function() {
  colors <- app_colors()
  bslib::page_sidebar(
    title = shiny::div(class = "ppm-header", shiny::span("Promoter Primer Mapper"), shiny::span(class = "ppm-subtitle", "Bisulfite PCR design")),
    window_title = "Promoter Primer Mapper", fillable = FALSE,
    theme = app_theme(colors), lang = "en",
    sidebar = bslib::sidebar(id = "setup_sidebar", title = "Reference and gene", width = 320,
      open = list(desktop = "open", mobile = "closed"),
      shiny::h5("1. Reference files"),
      shiny::textInput("gtf_path", "Annotation path (GTF/GFF3)", placeholder = "/path/to/annotation.gtf"),
      shiny::textInput("fasta_path", "Reference FASTA path", placeholder = "/path/to/genome.fa"),
      shiny::textInput("fasta_index_path", "FASTA index path (.fai)", placeholder = "/path/to/genome.fa.fai"),
      shiny::helpText("Files are read by the machine running the app. Container paths must be bind-mounted. An existing .fai index is required."),
      shiny::uiOutput("reference_readiness"),
      shiny::h5("2. Gene and TSS"),
      shiny::textInput("gene_query", "Gene ID or gene name"),
      app_action("find_gene", "Find gene"),
      shiny::selectInput("tss_choice", "TSS", choices = character(0)),
      shiny::tags$details(class = "ppm-details", shiny::tags$summary("Transcript details"), scroll_table("tss_table", "Transcript TSS records")),
      shiny::h5("3. Extract promoter"),
      shiny::numericInput("upstream", "Upstream (bp)", value = 1500, min = 0, step = 1),
      shiny::numericInput("downstream", "Downstream (bp)", value = 500, min = 0, step = 1),
      app_action("extract", "Extract and design"),
      shiny::helpText("Extracts the promoter and designs primers for all CpG islands, with non-island fallback."),
      shiny::tags$details(class = "ppm-details", shiny::tags$summary("Design tool settings"),
        shiny::textInput("primer3_path", "primer3_core path", placeholder = "blank = search PATH"),
        shiny::uiOutput("primer3_readiness"))),
    bsp_view_css(),
    shiny::uiOutput("design_context"), shiny::uiOutput("result_state"), shiny::uiOutput("status"),
    shiny::div(id = "operation_status", role = "status", `aria-live` = "polite", class = "ppm-help"),
    shiny::conditionalPanel("output.has_extraction",
      bslib::card(
        bslib::card_header("1. Promoter and design target"),
        shiny::uiOutput("region_summary"),
        shiny::plotOutput("cpg_map", height = "340px",
          brush = shiny::brushOpts("cpg_brush", direction = "x", resetOnNew = FALSE, fill = "transparent", stroke = colors[["primary"]],
                                   opacity = 1, delayType = "debounce", delay = 200),
          click = "cpg_click", hover = shiny::hoverOpts("cpg_hover", delay = 80, delayType = "throttle")),
        shiny::div(class = "ppm-help ppm-map-help", shiny::textOutput("cursor_position")),
        shiny::div(class = "ppm-actions",
          shiny::actionButton("zoom_target", "Zoom to target", class = "btn-outline-secondary btn-sm"),
          shiny::actionButton("zoom_pair", "Zoom to selected pair", class = "btn-outline-secondary btn-sm"),
          shiny::actionButton("zoom_reset", "Show full promoter", class = "btn-outline-secondary btn-sm")),
        shiny::uiOutput("map_context"),
        shiny::div(class = "ppm-target-grid",
          shiny::div(
            shiny::radioButtons("bsp_target_mode", "Target", c("All CpG islands" = "auto", "Selected island" = "island", "Selected range" = "range"), inline = TRUE),
            shiny::conditionalPanel("input.bsp_target_mode == 'island'", shiny::selectInput("bsp_island", "CpG island", choices = character(0))),
            shiny::conditionalPanel("input.bsp_target_mode == 'range'",
              shiny::div(class = "ppm-filter-grid",
                shiny::numericInput("range_from", "Range from (TSS-relative)", value = -50, step = 1),
                shiny::numericInput("range_to", "Range to (TSS-relative)", value = 50, step = 1)),
              shiny::div(class = "ppm-actions", app_action("use_selection", "Use range as target", "btn-outline-secondary btn-sm"),
                         shiny::actionButton("clear_selection", "Clear selection", class = "btn-outline-secondary btn-sm")),
              shiny::uiOutput("selection_table"),
              shiny::helpText("The TSS is +1; the base upstream is -1. There is no position 0."))),
          shiny::div(
            shiny::radioButtons("bsp_strand", "Template strand", c("Sense" = "sense", "Antisense" = "antisense", "Both" = "both"), selected = "both", inline = TRUE),
            shiny::uiOutput("target_validation"), scroll_table("target_preview", "Proposed design targets"),
            shiny::helpText("Every amplicon covers its whole target. All-island mode uses non-island CpG targets if no island yields a pair."),
            shiny::tags$details(id = "advanced_parameters", class = "ppm-details", shiny::tags$summary("Advanced design parameters"), bsp_params_ui()),
            shiny::uiOutput("parameter_validation"),
            shiny::uiOutput("design_readiness"),
            shiny::div(class = "ppm-actions", app_action("design_primers", "Design BSP primers", "btn-primary")))),
        shiny::tags$details(class = "ppm-details", shiny::tags$summary("CpG island statistics"), scroll_table("cpg_islands", "CpG island statistics"))),
      bslib::card(
        bslib::card_header("2. Compare primer pairs"),
        shiny::uiOutput("pair_state"),
        scroll_table("bsp_targets", "Completed design targets"),
        shiny::conditionalPanel("output.has_pairs",
          shiny::div(class = "ppm-review-grid",
            shiny::div(
              shiny::div(class = "ppm-filter-grid",
                shiny::selectInput("pair_target_filter", "Target filter", c("All targets" = "all")),
                shiny::selectInput("pair_strand_filter", "Strand filter", c("Both strands" = "all", "Sense" = "sense", "Antisense" = "antisense")),
                shiny::selectInput("pair_sort", "Sort by", c("Pair rank" = "rank", "Score" = "score", "Product size" = "product_size", "CpG count" = "cpg_in_amplicon"))),
              shiny::uiOutput("bsp_table"),
              shiny::helpText("Score combines the Primer3 pair penalty and Tm difference; lower is better. It is not a quality probability.")),
            shiny::div(
              shiny::selectInput("bsp_pair", "Selected pair", choices = character(0)),
              shiny::uiOutput("bsp_pair_detail"),
              shiny::tags$details(class = "ppm-details", open = "open", shiny::tags$summary("Amplicon sequence"),
                shiny::checkboxInput("sequence_full", "Show full promoter", FALSE), shiny::uiOutput("bsp_sequence_view"))))),
        shiny::uiOutput("bsp_hints"),
        shiny::tags$details(class = "ppm-details", shiny::tags$summary("Why targets have no pair"), scroll_table("bsp_diagnostics", "Candidate rejection diagnostics"))),
      bslib::card(
        bslib::card_header("3. Check specificity"),
        shiny::uiOutput("specificity_state"),
        shiny::textInput("bismark_genome", "Bismark genome folder", placeholder = "/path/to/genome_folder"),
        shiny::tags$details(class = "ppm-details", shiny::tags$summary("Alignment tool and criteria"),
          shiny::textInput("bowtie2_path", "bowtie2 path", placeholder = "blank = search PATH"),
          shiny::div(class = "ppm-parameter-grid",
            shiny::numericInput("spec_top_n", "Pairs per target", 5, min = 1, step = 1),
            shiny::numericInput("spec_max_mismatch", "Max mismatches per hit", 2, min = 0, max = 5, step = 1),
            shiny::numericInput("spec_max_product", "Max product (bp)", 3000, min = 50, step = 50))),
        shiny::uiOutput("specificity_readiness"), app_action("check_specificity", "Check specificity"),
        shiny::helpText("Optional: requires Bowtie2 indexes prepared by Bismark. Checks the best pairs per target on both converted strands, under the displayed criteria.")),
      bslib::card(
        bslib::card_header("4. Export results"),
        shiny::uiOutput("export_context"),
        shiny::conditionalPanel("output.has_pairs",
          shiny::h5("All pairs"),
          shiny::div(class = "ppm-actions", shiny::downloadButton("download_bsp_pairs_csv", "Primers CSV")),
          shiny::h5("Selected pair"),
          shiny::div(class = "ppm-actions", shiny::downloadButton("download_order_sheet", "Selected primers CSV"), shiny::downloadButton("download_amplicons", "Amplicon references")),
          shiny::helpText("Amplicon FASTA contains genomic, converted-unmethylated and converted-methylated references.")),
        shiny::tags$details(class = "ppm-details", shiny::tags$summary("Extracted promoter and FASTA variants"),
          shiny::textOutput("sequence_header"), shiny::verbatimTextOutput("sequence_text"),
          shiny::div(class = "ppm-actions", shiny::downloadButton("download_fasta", "Sense FASTA"), shiny::downloadButton("download_complementary_fasta", "Antisense FASTA"),
            shiny::downloadButton("download_bisulfite_requested_fasta", "Converted sense"), shiny::downloadButton("download_bisulfite_complementary_fasta", "Converted antisense")),
          shiny::helpText("Converted promoter exports replace every C with T."))))
  )
}
