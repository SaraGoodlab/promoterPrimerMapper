# Shared GUI and map tokens. Plum marks actions and selections; biological
# features keep distinct slate, teal, and amber colors.
app_colors <- function() {
  c(background = "#F7F7F8", surface = "#FFFFFF", border = "#E3E3E8",
    ink = "#1C1B22", muted = "#5C5B66", control = "#C9C8D0",
    primary = "#6B3F5A", `primary-hover` = "#57324A", `primary-soft` = "#F3ECF0",
    island = "#3F7F6E", cpg = "#44435A", primer = "#D98E04",
    warning = "#B45309", danger = "#B42318", line = "#C9C8D0")
}

app_theme <- function(colors = app_colors()) {
  bslib::bs_theme(version = 5,
    bg = colors[["background"]], fg = colors[["ink"]],
    base_font = bslib::font_collection("IBM Plex Sans", "Inter", "system-ui", "sans-serif"),
    code_font = bslib::font_collection("IBM Plex Mono", "JetBrains Mono", "monospace"),
    primary = colors[["primary"]], secondary = colors[["muted"]],
    success = colors[["muted"]], info = colors[["primary"]],
    warning = colors[["warning"]], danger = colors[["danger"]],
    "font-size-base" = "1rem", "font-weight-base" = "400", "font-weight-bold" = "600",
    "headings-font-weight" = "600", "card-bg" = colors[["surface"]], "card-cap-bg" = colors[["background"]],
    "input-bg" = colors[["surface"]], "input-border-color" = colors[["control"]],
    "input-placeholder-color" = colors[["muted"]], "border-color" = colors[["border"]],
    "body-secondary-color" = colors[["muted"]], "link-color" = colors[["ink"]],
    "link-hover-color" = colors[["ink"]], "input-btn-focus-box-shadow" = "none")
}

app_assets <- function() {
  htmltools::htmlDependency("ppm-gui", "0.2.0",
    src = c(file = system.file("app", package = "promoterPrimerMapper")),
    stylesheet = "fonts.css", script = "gui.js", all_files = TRUE)
}
