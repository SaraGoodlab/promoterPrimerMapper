# Promoter Primer Mapper

An R package with a local Shiny GUI for selecting transcript-derived transcription start sites (TSSs), extracting promoter sequences, and mapping primer positions back to genomic coordinates.

## Install and run

Install the required packages in R:

```r
install.packages(c("BiocManager", "shiny", "bslib", "testthat"))
BiocManager::install(c(
  "rtracklayer",
  "GenomicRanges",
  "IRanges",
  "Biostrings",
  "Rsamtools"
))
```

Install the package from the project directory:

```bash
R CMD INSTALL .
```

Then launch the GUI in R:

```r
promoterPrimerMapper::run_app()
```

## GitHub release installation

For a published GitHub Release, install the Bioconductor dependencies first, then download the release asset named `promoterPrimerMapper_0.1.0.tar.gz` and install it locally:

```bash
R CMD INSTALL promoterPrimerMapper_0.1.0.tar.gz
```

Alternatively, after replacing `<OWNER>/<REPOSITORY>` with the GitHub repository path:

```r
install.packages("remotes")
remotes::install_url(
  "https://github.com/<OWNER>/<REPOSITORY>/releases/download/v0.1.0/promoterPrimerMapper_0.1.0.tar.gz",
  dependencies = FALSE
)
promoterPrimerMapper::run_app()
```

Maintainers can build and validate the release asset locally with:

```bash
R CMD build .
R CMD check promoterPrimerMapper_0.1.0.tar.gz
```

Create the GitHub Release as tag `v0.1.0` and attach the resulting `.tar.gz` file.
This repository also includes a GitHub Actions workflow at `.github/workflows/release.yml`; after pushing tag `v0.1.0`, it builds the tarball and creates the release automatically:

```bash
git add .
git commit -m "Prepare v0.1.0 release"
git tag v0.1.0
git push origin main --tags
```

The GUI requires paths to a GTF file, FASTA file, and pre-generated FASTA `.fai` index.

Primer coordinates are mapped without sequence alignment. Mapping uses the extracted sequence coordinates, genomic extraction interval, strand, and primer `Start` and `Size`.

After extraction, four FASTA outputs are available:

* transcription-oriented sequence
* reverse-complement sequence
* bisulfite-converted transcription-oriented sequence
* bisulfite-converted reverse-complement sequence

Bisulfite conversion replaces all `C` with `T`.

## Coordinate conventions

Genomic and extracted-sequence coordinates are 1-based and inclusive.

Primer `Start` values are assumed to be 1-based by default. For tools reporting zero-based starts, select `0-based` in the GUI or use:

```r
map_primer_coordinates(..., start_convention = "zero-based")
```

For `+` strand regions, extracted sequence position 1 corresponds to `region_start`. For `-` strand regions, the sequence is reverse-complemented and position 1 corresponds to `region_end`. Primer intervals are reported in ascending genomic order.

The TSS is defined as the transcript start on `+` and transcript end on `-`. Transcripts sharing the same gene, scaffold, strand, and TSS are collapsed into a single selectable TSS while retaining all transcript IDs.

## Tests

Run from the package root:

```r
devtools::test()
```

The synthetic GTF and FASTA under `inst/extdata/` are used by the automated integration test.
