test_that("every BSP parameter has a GUI input defaulting to bsp_primer_params()", {
  html <- as.character(bsp_params_ui())
  defaults <- bsp_primer_params()
  for (name in names(defaults)) {
    expect_match(html, sprintf("id=\"%s\"", bsp_param_input_id(name)), fixed = TRUE)
    expect_match(html, sprintf("id=\"%s\"[^>]*value=\"%s\"", bsp_param_input_id(name), defaults[[name]]))
  }
})

test_that("GUI parameter inputs round-trip into Primer3 tags and reset to defaults", {
  shiny::testServer(app_server, {
    defaults <- bsp_primer_params()
    inputs <- stats::setNames(defaults, vapply(names(defaults), bsp_param_input_id, character(1)))
    do.call(session$setInputs, inputs)
    expect_identical(bsp_params(), defaults)
    session$setInputs(bsp_tm_opt = 60, bsp_product_size_max = 400)
    expect_equal(primer3_global_tags(bsp_params())$PRIMER_OPT_TM, 60)
    expect_equal(primer3_global_tags(bsp_params())$PRIMER_PRODUCT_SIZE_RANGE, "200-400")
    session$setInputs(bsp_tm_min = 70, design_primers = 1)
    expect_match(message_text(), "Extract a sequence")
  })
})

test_that("reset restores every parameter input", {
  session <- shiny::MockShinySession$new()
  sent <- list()
  session$sendInputMessage <- function(inputId, message) sent[[inputId]] <<- message$value
  shiny::testServer(app_server, session = session, {
    session$setInputs(bsp_tm_opt = 61, bsp_reset = 1)
  })
  defaults <- bsp_primer_params()
  expect_setequal(names(sent), vapply(names(defaults), bsp_param_input_id, character(1)))
  expect_equal(vapply(names(defaults), function(n) as.numeric(sent[[bsp_param_input_id(n)]]), numeric(1)), unlist(defaults))
})

test_that("the import tab is gone and sequence and design share one screen", {
  html <- as.character(app_ui())
  expect_false(grepl("primer_text|id=\"convert\"|Import primers", html))
  expect_false(grepl("demo_gene|load_demo|Load demo", html))
  expect_false(grepl("tabbable|nav-tabs", html))
  for (id in c("cpg_map", "range_from", "range_to", "use_selection", "bsp_target_mode", "bsp_island", "design_primers", "bsp_table")) {
    expect_match(html, sprintf("id=\"%s\"", id), fixed = TRUE)
  }
})

test_that("brushing the map sets the TSS-relative range and clicks pick islands", {
  session <- shiny::MockShinySession$new()
  sent <- list()
  session$sendInputMessage <- function(inputId, message) sent[[inputId]] <<- message
  shiny::testServer(app_server, session = session, {
    extracted(fixture_promoter("BSPPLUS"))
    tss <- tss_sense_position(extracted())
    session$setInputs(cpg_brush = list(xmin = tss - 320.4, xmax = tss + 44.6))
    expect_equal(as.numeric(sent$range_from$value), -320)
    expect_equal(as.numeric(sent$range_to$value), 46)
    session$setInputs(range_from = -320, range_to = 46)
    expect_equal(selection(), c(tss - 320L, tss + 45L))
    expect_match(output$selection_table, "chrTest1:2681-3046 (+)", fixed = TRUE)
    session$setInputs(range_from = 0)
    expect_match(selection(), "0 is not a valid")
    session$setInputs(range_from = -320, use_selection = 1)
    expect_equal(sent$bsp_target_mode$value, "range")
    expect_match(message_text(), "selected range -320..+46", fixed = TRUE)
    island <- cpg_island_scan(extracted())[1L, ]
    session$setInputs(cpg_click = list(x = island$start + 10, y = 0.9))
    expect_null(sent$bsp_island)
    session$setInputs(cpg_click = list(x = island$start + 10, y = 0))
    expect_equal(sent$bsp_island$value, "1")
    expect_equal(sent$bsp_target_mode$value, "island")
    session$setInputs(cpg_hover = list(x = tss, y = 0.5))
    expect_match(output$cursor_position, "^Cursor: \\+1 ")
    session$setInputs(bsp_target_mode = "island", bsp_island = "2")
    expect_match(output$target_preview, "CGI2")
    session$setInputs(bsp_target_mode = "range", range_from = -600, range_to = 50)
    expect_match(output$target_preview, "can be at most 406 bp")
  })
})

test_that("extraction designs all CpG islands and the other target modes design on demand", {
  skip_if_no_primer3()
  session <- shiny::MockShinySession$new()
  sent <- list()
  session$sendInputMessage <- function(inputId, message) sent[[inputId]] <<- message
  shiny::testServer(app_server, session = session, {
    session$setInputs(gtf_path = fixture_reference("gtf"), gene_query = "BSPPLUS", find_gene = 1)
    fasta <- fixture_reference("fa")
    session$setInputs(fasta_path = fasta, fasta_index_path = paste0(fasta, ".fai"), tss_choice = "1",
                      upstream = 1000, downstream = 500, extract = 1)
    expect_match(message_text(), "^Extracted 1501 bp .* Designed [0-9]+ BSP primer pairs for CGI1 .*CGI2")
    expect_equal(design()$targets$target_id, c("CGI1", "CGI2"))
    expect_equal(sent$bsp_target_mode$value, "auto")
    expect_equal(sent$bsp_island$value, "2")
    expect_match(sent$bsp_island$options, "CGI1: -544..-211", fixed = TRUE)
    expect_match(output$bsp_targets, "CGI2")
    session$setInputs(bsp_strand = "sense", bsp_target_mode = "range", range_from = -50, range_to = 50,
                      bsp_num_return = 3, design_primers = 1)
    expect_match(message_text(), "Designed 3 BSP primer pairs for R1 (range, -50..+50)", fixed = TRUE)
    expect_equal(design()$params$num_return, 3L)
    session$setInputs(bsp_pair = "1")
    expect_match(as.character(output$bsp_sequence_view$html), "class=\"primer\"")
    expect_match(output$bsp_table, design()$pairs$left_sequence[1L])
    session$setInputs(bsp_target_mode = "island", bsp_island = "1", design_primers = 2)
    expect_equal(unique(design()$pairs$target_id), "CGI1")
    session$setInputs(bsp_tm_min = 70, design_primers = 3)
    expect_match(message_text(), "^Primer design failed: Tm values must satisfy")
  })
})
