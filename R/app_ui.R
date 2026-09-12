bsp_param_input_id <- function(name) paste0("bsp_", name)

# Shiny inputs for every bsp_primer_params() field, grouped as in
# bsp_param_spec(). Defaults come from bsp_primer_params() so the GUI and the R
# API cannot drift apart.
bsp_params_ui <- function(params = bsp_primer_params()) {
  spec <- bsp_param_spec()
  groups <- unique(vapply(spec, `[[`, character(1), "group"))
  panels <- lapply(groups, function(group) {
    names <- names(spec)[vapply(spec, function(s) identical(s$group, group), logical(1))]
    inputs <- lapply(names, function(name) {
      s <- spec[[name]]
      shiny::column(4,
        shiny::numericInput(bsp_param_input_id(name), s$label, value = params[[name]], min = s$min, max = s$max, step = s$step),
        if (!is.null(s$help)) shiny::tags$small(class = "text-muted d-block mb-2", s$help))
    })
    bslib::accordion_panel(group, do.call(shiny::fluidRow, inputs))
  })
  shiny::tagList(
    do.call(bslib::accordion, c(list(id = "bsp_param_groups", open = FALSE), panels)),
    shiny::actionButton("bsp_reset", "Reset parameters to defaults", class = "btn-outline-secondary btn-sm mt-2")
  )
}

bsp_view_css <- function() {
  shiny::tags$style(shiny::HTML(paste(
    ".bsp-seq { font-size: 0.85em; line-height: 1.5; white-space: pre; overflow-x: auto; }",
    ".bsp-seq .amp, .bsp-legend .amp { background: #eaf2fb; }",
    ".bsp-seq .primer, .bsp-legend .primer { background: #ffe08a; text-decoration: underline; }",
    ".bsp-seq .cpg, .bsp-legend .cpg { color: #c0392b; font-weight: bold; }",
    ".bsp-seq .tss, .bsp-legend .tss { outline: 1px solid #2c3e50; }",
    ".bsp-legend span { padding: 0 4px; margin-right: 8px; }",
    ".card { margin-bottom: 1rem; }"
  )))
}

app_ui <- function() {
  shiny::fluidPage(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    bsp_view_css(),
    shiny::titlePanel("Promoter Primer Mapper"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width = 3,
        shiny::h4("1. Reference files"),
        shiny::textInput("gtf_path", "GTF annotation path", placeholder = "/path/to/annotation.gtf"),
        shiny::textInput("fasta_path", "Reference FASTA path", placeholder = "/path/to/genome.fa"),
        shiny::textInput("fasta_index_path", "FASTA index path (.fai)", placeholder = "/path/to/genome.fa.fai"),
        shiny::helpText("Files are opened directly from this computer. An existing .fai index is required."),
        shiny::tags$details(
          shiny::tags$summary("External tools"),
          shiny::textInput("primer3_path", "primer3_core path", placeholder = "blank = search PATH"),
          shiny::textInput("bowtie2_path", "bowtie2 path", placeholder = "blank = search PATH"),
          shiny::textInput("bismark_genome", "Bismark genome folder", placeholder = "/path/to/genome_folder"),
          shiny::helpText("The Bismark folder must contain Bisulfite_Genome/ from 'bismark_genome_preparation --bowtie2'. Only needed for the specificity check.")
        ),
        shiny::hr(),
        shiny::h4("2. Select gene and TSS"),
        shiny::textInput("gene_query", "Gene ID or gene name"),
        shiny::actionButton("find_gene", "Find gene", class = "btn-primary"),
        shiny::selectInput("tss_choice", "TSS", choices = character(0)),
        shiny::hr(),
        shiny::h4("3. Extract sequence"),
        shiny::numericInput("upstream", "Upstream (bp)", value = 1500, min = 0, step = 1),
        shiny::numericInput("downstream", "Downstream (bp)", value = 500, min = 0, step = 1),
        shiny::actionButton("extract", "Extract sequence", class = "btn-primary"),
        shiny::helpText("Extracting a sequence designs primers for all of its CpG islands straight away.")
      ),
      shiny::mainPanel(width = 9,
        shiny::uiOutput("status"),
        bslib::card(
          bslib::card_header("Promoter CpG map"),
          shiny::uiOutput("region_summary"),
          shiny::plotOutput("cpg_map", height = "280px",
            brush = shiny::brushOpts("cpg_brush", direction = "x", resetOnNew = TRUE, fill = "#2980b9", stroke = "#2980b9",
                                     opacity = 0.15, delayType = "debounce", delay = 5000),
            click = "cpg_click",
            hover = shiny::hoverOpts("cpg_hover", delay = 80, delayType = "throttle")),
          shiny::div(class = "small text-muted", shiny::textOutput("cursor_position", inline = TRUE)),
          shiny::helpText("Drag across the map to select a range; click a CpG island bar (green) to design primers for that island."),
          shiny::fluidRow(
            shiny::column(4,
              shiny::h5("Selected range"),
              shiny::fluidRow(
                shiny::column(6, shiny::numericInput("range_from", "From", value = -50, step = 1)),
                shiny::column(6, shiny::numericInput("range_to", "To", value = 50, step = 1))),
              shiny::actionButton("use_selection", "Use selection as BSP target", class = "btn-outline-primary btn-sm")),
            shiny::column(8, shiny::tableOutput("selection_table"))
          ),
          shiny::helpText("TSS-relative positions: the TSS is +1 and the base upstream is -1 (there is no 0)."),
          shiny::h5("CpG islands"),
          shiny::tableOutput("cpg_islands")
        ),
        bslib::card(
          bslib::card_header("BSP primer design"),
          shiny::fluidRow(
            shiny::column(4, shiny::radioButtons("bsp_strand", "Template strand",
              c("Sense" = "sense", "Antisense" = "antisense", "Both" = "both"), selected = "both", inline = TRUE)),
            shiny::column(8, shiny::radioButtons("bsp_target_mode", "Target",
              c("All CpG islands" = "auto", "Selected CpG island" = "island", "Selected range" = "range"), inline = TRUE))
          ),
          shiny::conditionalPanel("input.bsp_target_mode == 'island'",
            shiny::selectInput("bsp_island", "CpG island", choices = character(0), width = "100%")),
          shiny::helpText("Every amplicon contains its whole target. 'All CpG islands' designs each island's window and, when no island yields a pair, the densest non-island CpG regions."),
          shiny::tableOutput("target_preview"),
          bsp_params_ui(),
          shiny::div(class = "mt-3", shiny::actionButton("design_primers", "Design BSP primers", class = "btn-success"))
        ),
        bslib::card(
          bslib::card_header("Primer pairs"),
          shiny::tableOutput("bsp_targets"),
          shiny::tableOutput("bsp_table"),
          shiny::uiOutput("bsp_hints"),
          shiny::tableOutput("bsp_diagnostics"),
          shiny::selectInput("bsp_pair", "Selected pair", choices = character(0), width = "100%"),
          shiny::tableOutput("bsp_pair_detail"),
          shiny::uiOutput("bsp_sequence_view"),
          shiny::div(
            shiny::downloadButton("download_bsp_pairs", "Download pairs (TSV)"),
            shiny::downloadButton("download_bsp_pairs_csv", "Download pairs (CSV)"),
            shiny::downloadButton("download_order_sheet", "Download order sheet (selected pair)"),
            shiny::downloadButton("download_amplicons", "Download amplicon references (selected pair)")),
          shiny::h5("Specificity", class = "mt-3"),
          shiny::fluidRow(
            shiny::column(4, shiny::numericInput("spec_top_n", "Pairs to check per target", value = 5, min = 1, step = 1)),
            shiny::column(4, shiny::numericInput("spec_max_mismatch", "Max mismatches per hit", value = 2, min = 0, max = 5, step = 1)),
            shiny::column(4, shiny::numericInput("spec_max_product", "Max off-target product (bp)", value = 3000, min = 50, step = 50))
          ),
          shiny::div(shiny::actionButton("check_specificity", "Check specificity (Bowtie2)", class = "btn-primary")),
          shiny::helpText("Aligns the best pairs of each target to the bisulfite-converted genome and predicts every product they could amplify.")
        ),
        bslib::card(
          bslib::card_header("Extracted sequence"),
          shiny::tags$details(
            shiny::tags$summary("Show sequence and FASTA downloads"),
            shiny::tableOutput("tss_table"),
            shiny::textOutput("sequence_header"),
            shiny::verbatimTextOutput("sequence_text"),
            shiny::downloadButton("download_fasta", "Download requested FASTA"),
            shiny::downloadButton("download_complementary_fasta", "Download complementary FASTA"),
            shiny::downloadButton("download_bisulfite_requested_fasta", "Download bisulfite requested FASTA"),
            shiny::downloadButton("download_bisulfite_complementary_fasta", "Download bisulfite complementary FASTA"))
        )
      )
    )
  )
}
