app_ui <- function() {
  shiny::fluidPage(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    shiny::titlePanel("Promoter Primer Mapper"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::h4("1. Reference files"),
        shiny::textInput("gtf_path", "GTF annotation path", placeholder = "/path/to/annotation.gtf"),
        shiny::textInput("fasta_path", "Reference FASTA path", placeholder = "/path/to/genome.fa"),
        shiny::textInput("fasta_index_path", "FASTA index path (.fai)", placeholder = "/path/to/genome.fa.fai"),
        shiny::helpText("Files are opened directly from this computer. An existing .fai index is required."),
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
        shiny::hr(),
        shiny::h4("4. Primer coordinates"),
        shiny::radioButtons("primer_mode", "Input mode", c("Paste table" = "paste", "Manual entry" = "manual"), inline = TRUE),
        shiny::conditionalPanel("input.primer_mode == 'paste'", shiny::textAreaInput("primer_text", "Primer Start/Size table", rows = 5, placeholder = "Left primer    1074    22\nRight primer   1299    28")),
        shiny::conditionalPanel("input.primer_mode == 'manual'",
          shiny::numericInput("left_start", "Left primer start", value = NA, min = 0),
          shiny::numericInput("left_size", "Left primer size", value = NA, min = 1),
          shiny::numericInput("right_start", "Right primer start", value = NA, min = 0),
          shiny::numericInput("right_size", "Right primer size", value = NA, min = 1)),
        shiny::selectInput("start_convention", "Primer Start convention", c("1-based" = "one-based", "0-based" = "zero-based")),
        shiny::actionButton("convert", "Convert coordinates", class = "btn-success")
      ),
      shiny::mainPanel(
        shiny::uiOutput("status"),
        shiny::h4("Selected TSS"),
        shiny::tableOutput("tss_table"),
        shiny::h4("Extracted sequence"),
        shiny::textOutput("sequence_header"),
        shiny::verbatimTextOutput("sequence_text"),
        shiny::downloadButton("download_fasta", "Download requested FASTA"),
        shiny::downloadButton("download_complementary_fasta", "Download complementary FASTA"),
        shiny::downloadButton("download_bisulfite_requested_fasta", "Download bisulfite requested FASTA"),
        shiny::downloadButton("download_bisulfite_complementary_fasta", "Download bisulfite complementary FASTA"),
        shiny::h4("Primer coordinates"),
        shiny::tableOutput("primer_table"),
        shiny::textOutput("amplicon_interval"),
        shiny::downloadButton("download_results_csv", "Download CSV"),
        shiny::downloadButton("download_results", "Download TSV")
      )
    )
  )
}
