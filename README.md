# Promoter Primer Mapper

An R package with a local Shiny GUI for designing **bisulfite-sequencing PCR (BSP)** primers for gene promoters. It selects transcript-derived transcription start sites (TSSs) and extracts a transcription-oriented promoter sequence. An interactive CpG map shows the promoter's CpG sites and CpG islands. The app designs ranked BSP primer pairs on the bisulfite-converted template for every CpG island, or for a range or island you pick on the map. It checks the pairs for off-target products in the converted genome, and exports an oligo order sheet together with amplicon references for analysing the sequencing reads.

## Run the app

Choose one of the following paths:

| Where | Direct R installation | Apptainer/Singularity container |
|---|---|---|
| Your computer | Use when R and the required tools can be installed locally. | Use when you want a self-contained, reproducible environment. |
| HPC | Use the cluster's R and software modules or Conda environment. | Usually the simplest HPC option when Apptainer/Singularity is available. |

The app listens on `127.0.0.1` by default. Keep that setting on shared systems and access an HPC job through an SSH tunnel.

### Local computer: direct R installation

Install R and the package dependencies. Bioconductor packages must be installed with an R version supported by the current Bioconductor release.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "rtracklayer", "GenomicRanges", "IRanges", "Biostrings", "Rsamtools",
  "shiny", "bslib"
), ask = FALSE, update = FALSE)
```

Install the command-line tools used by the app:

| Tool | Required for | Example installation |
|---|---|---|
| [Primer3](https://primer3.org) `primer3_core` (>= 2.5.0) | Primer design | `conda install -c bioconda primer3` |
| [Bowtie2](https://github.com/BenLangmead/bowtie2) and [Bismark](https://github.com/FelixKrueger/Bismark) | Optional specificity check | `conda install -c bioconda bowtie2 bismark` |

On Debian/Ubuntu, `sudo apt install primer3 bowtie2` supplies the first two tools; install Bismark with Conda or your system package manager. Confirm that the commands are available before starting the app:

```bash
command -v primer3_core
command -v bowtie2
command -v bismark_genome_preparation
```

Clone and install the package, then launch the GUI:

```bash
git clone https://github.com/SaraGoodlab/promoterPrimerMapper.git
cd promoterPrimerMapper
R CMD INSTALL .
R --vanilla -e 'promoterPrimerMapper::run_app()'
```

Open <http://127.0.0.1:3838> if a browser does not open automatically. You may instead install a published release asset with `R CMD INSTALL promoterPrimerMapper_<version>.tar.gz`.

The app looks for `primer3_core` and `bowtie2` on `PATH`. If necessary, supply their full paths in the GUI under **External tools**, or before launching R:

```r
options(
  promoterPrimerMapper.primer3 = "/path/to/primer3_core",
  promoterPrimerMapper.bowtie2 = "/path/to/bowtie2"
)
```

### HPC: direct R installation

Run the app in an interactive compute allocation, not on a login node. The module names and scheduler commands below are examples; use your cluster's documented equivalents.

```bash
# On the login node: obtain an interactive allocation.
salloc --time=02:00:00 --mem=8G --cpus-per-task=2

# On the allocated compute node: load locally provided software.
module load R primer3 bowtie2 bismark

# Install R packages where you have write access, once per R version.
export R_LIBS_USER="$HOME/R/library"
mkdir -p "$R_LIBS_USER"
git clone https://github.com/SaraGoodlab/promoterPrimerMapper.git
cd promoterPrimerMapper
R CMD INSTALL --library="$R_LIBS_USER" .

# Keep the app private to the compute node.
R --vanilla -e 'promoterPrimerMapper::run_app(host = "127.0.0.1", port = 3838, launch.browser = FALSE)'
```

From your own computer, forward the port through the login node. Replace `node123` with the allocated node name from `squeue`, `srun hostname`, or the scheduler's equivalent.

```bash
ssh -N -J user@login.cluster -L 3838:127.0.0.1:3838 user@node123
```

Then open <http://127.0.0.1:3838> locally. If the compute node is not reachable through SSH, follow your cluster's web-proxy or port-forwarding instructions.

### Local computer: container

The Apptainer/Singularity image includes R, this package, Primer3, Bowtie2, and Bismark. Build it once from the repository root; building requires root or Apptainer's `--fakeroot` support.

```bash
apptainer build --fakeroot promoterPrimerMapper.sif container/promoterPrimerMapper.def
```

Use `singularity` in place of `apptainer` on systems that provide that command. The build runs the package tests; add `--notest` to skip them, or run them later with `apptainer test promoterPrimerMapper.sif`.

Bind the directory that holds the annotation, reference FASTA, `.fai`, and any Bismark genome folder. The bind preserves the same path inside the container, so paste that host path into the GUI.

```bash
apptainer run --bind /path/to/genomes promoterPrimerMapper.sif
# Use another local port if 3838 is already occupied.
apptainer run --bind /path/to/genomes promoterPrimerMapper.sif --port 4000
```

Open <http://127.0.0.1:3838> (or port `4000`) in your browser. Leave **External tools** empty: the image already places `primer3_core` and `bowtie2` on `PATH`.

Prepare a Bismark index with the same image when you need specificity checking:

```bash
apptainer exec --bind /path/to/genomes promoterPrimerMapper.sif \
  bismark_genome_preparation --bowtie2 /path/to/genomes/hg38
```

### HPC: container

Copy a previously built `.sif` to the cluster if the cluster prohibits image builds. Then request an interactive allocation and run the image on its compute node:

```bash
# On the login node
salloc --time=02:00:00 --mem=8G --cpus-per-task=2

# On the allocated compute node
apptainer run --bind /project,/scratch promoterPrimerMapper.sif
```

All paths entered in the GUI must be available inside the container. Bind every directory that contains reference data or output files, for example `--bind /project/mygroup,/scratch/$USER`. The container already contains the external tools, so no modules or host R libraries are needed.

From your computer, create the same SSH tunnel used for a direct installation:

```bash
ssh -N -J user@login.cluster -L 3838:127.0.0.1:3838 user@node123
```

Open <http://127.0.0.1:3838> in your local browser. Do not use `--host 0.0.0.0` on a shared network: the app can read every file available to your account.

## Workflow in the GUI

The sidebar takes you from reference files to an extracted promoter. Everything else is on one screen.

1. **Reference files.** Give paths to a GTF/GFF annotation, a reference FASTA and its `.fai` index. A GFF3 file works too: `mRNA` records are linked to their gene through `Parent`. The annotation is cached between gene searches.
2. **Gene and TSS.** Search by gene ID or gene name and pick a TSS.
3. **Extract sequence.** Set the bases to extract upstream and downstream of the TSS. Extraction immediately designs primers for **all CpG islands** (see [How targets are chosen](#how-targets-are-chosen)).

**Promoter CpG map.** The map shows every CpG (red), the predicted CpG islands (green bars, labelled `CGI1`, `CGI2`, …) and the TSS. The bottom axis gives TSS-relative positions and the top axis genomic coordinates.
- **Hover** to read the TSS-relative, sequence and genomic position and the base under the cursor.
- **Drag** across the map to select a range. The range fills the *From*/*To* boxes, which you can also type into, and is described below the map: TSS-relative, sequence and genomic interval, length, CpG count, GC % and overlapping islands. **Use selection as BSP target** makes it the design target.
- **Click** an island bar to make that island the target.
- The island table lists every predicted island with its interval, CpG count, GC %, CpG observed/expected ratio and the window that will be designed.

**BSP primer design.**
- **Template strand:** sense, antisense or both.
- **Target:** *All CpG islands* (default), *Selected CpG island* (choose from the list of all islands), or *Selected range*. The targets about to be designed are listed and shaded on the map.
- **Parameters.** Every design parameter (below) is in the collapsible panels, pre-filled with the defaults from `bsp_primer_params()`. **Reset parameters to defaults** restores them.
- **Design BSP primers** re-runs the design with the current target and parameters.

**Primer pairs.**
- **Tables.** Results show the pairs found per target and a ranked pair table.
- **Targets without a pair** get candidate counts for each filtering stage, the main reasons candidates were rejected, and a hint on what to change. For example, *allow 1–2 CpGs per primer* where CpG-dense islands leave no CpG-free primer site.
- **Selected pair.** Selecting a pair shows per-primer statistics, draws the amplicon on the map, and shows the converted template with the primers, CpGs and TSS highlighted.
- **Downloads:**
  - the pair table (TSV/CSV), headed by the design region, the targets and every parameter used;
  - an order sheet for the selected pair, with oligos named like `GENE_CGI1_BSP3_F`;
  - amplicon reference FASTA (genomic, converted-unmethylated and converted-methylated) for BiQ Analyzer, QUMA or Bismark.
- **Check specificity** aligns the best pairs of each target to a Bismark-prepared genome and reports off-target products (see below).

**Extracted sequence** (collapsed): the TSS record, the sequence and the four FASTA downloads.

## Using the R API

```r
library(promoterPrimerMapper)
annotation <- load_annotation("annotation.gtf")
tss <- collapse_tss(get_transcript_tss(annotation, gene = "MLH1"))
promoter <- extract_tss_region(tss[1, ], "genome.fa", upstream = 1000, downstream = 500)

cpg_island_scan(promoter)                # predicted CpG islands
design <- design_bsp_primers(promoter, strand = "both",
                             params = list(product_size_max = 400, tm_opt = 60))
design                                   # targets and ranked pairs

# Other targets: one island, a TSS-relative range, or the non-island CpG regions
design_bsp_primers(promoter, target = island_targets(promoter, islands = 2))
design_bsp_primers(promoter, target = c(-300, 100))
design_bsp_primers(promoter, target = "non_island")

design <- check_bsp_specificity(design, "/path/to/bismark_genome", top_n = 5)
write_bsp_pairs(design, "pairs.tsv")
write_order_sheet(design, "order.tsv", ranks = 1)
write_amplicon_references(design, "amplicon.fa", ranks = 1)
plot_cpg_map(promoter, design, rank = 1)
```

`design$targets` lists the targets with the number of pairs found for each. If no pair passes, `design$pairs` is `NULL`. To see which constraint to relax:
- `design$diagnostics` shows how many candidate primers survived each stage for every target and strand, with the main rejection reasons per side.
- `design$rejections` holds the full counts.
- Printing the design adds a hint naming the main problem, e.g. CpGs under the primer, or sequence too AT-rich after conversion.

## How targets are chosen

Every amplicon contains its whole target, and the primers lie outside it, so all CpGs in the target are read.

- **CpG islands.** Islands are predicted with Gardiner-Garden/Frommer-style criteria: 100 bp windows with GC ≥ 50 % and CpG observed/expected ≥ 0.6, merged, at least 200 bp long. Each island gives one target (`CGI1`, `CGI2`, …).
  - An island no longer than the *target width* (mid product size minus two maximum-length primers; 265 bp with the defaults) is targeted whole.
  - A longer island gets a window of that width, centred on the TSS when the island contains it and on the island centre otherwise. To cover another part of a long island, select it as a range.
- **Non-island fallback.** In *All CpG islands* mode, if no island yields a primer pair, or the sequence has no island, the densest CpG clusters outside islands are designed instead (`NC1`, `NC2`, …).
  - Windows of the target width start at each non-island CpG and stop before the next island.
  - Each window is trimmed to its first and last CpG and must hold at least `min_cpg_in_amplicon` CpGs.
  - The windows with the most CpGs are taken, without overlap and nearest the TSS on ties, up to `fallback_targets`.
- **Selected range.** Any TSS-relative range up to `product_size_max − 2 × primer_size_min` bp (406 bp with the defaults). Narrower ranges leave more room for primers and give more pairs.

## How primers are designed

1. **Template.** The extracted sequence (sense) or its reverse complement (antisense) is converted in silico. Every non-CpG `C` becomes `T`. CpG cytosines are methylation-dependent and are tracked separately.
2. **Candidate enumeration.** Every primer footprint of allowed length on either side of the target is tested against the bisulfite rules:
   - no CpG cytosine under the primer, unless `max_cpg_in_primer` > 0; allowed CpGs must be in the 5′ half and are written as degenerate `Y` (forward) or `R` (reverse);
   - at least `min_converted_c` converted cytosines under the primer, so only fully converted DNA is amplified;
   - a converted cytosine within `three_prime_window` bases of the 3′ end;
   - no run longer than `max_poly_non_t` of a single A, C or G base.
3. **Primer3 scoring.** Surviving oligos are checked with `primer3_core` (`check_primers`) for Tm, GC, homopolymers, self-dimers and hairpins. The best `candidates_per_side` per side are kept.
4. **Pairing.** Pairs must fit the product size range and the Tm difference limit, and enclose at least `min_cpg_in_amplicon` CpGs. Each primer site must occur only once in its amplicon, so primers in tandem repeats are dropped. A shortlist is re-checked with Primer3 for cross-dimers.
5. **Ranking.** Pairs are ranked by Primer3 pair penalty plus Tm difference. No more than `max_primer_reuse` reported pairs may share a primer site; primers whose 3′ ends lie within half the minimum primer length count as one site.

| Parameter | Default | Why |
|---|---|---|
| `product_size_min`, `product_size_max` | 200–450 bp | Bisulfite treatment fragments DNA; Sanger/cloning works best with amplicons under ~500 bp |
| `primer_size_min/opt/max` | 22 / 26 / 30 nt | Converted DNA has a three-letter alphabet, so primers need extra length for specificity |
| `tm_min/opt/max` | 52 / 58 / 64 °C | Primer3 nearest-neighbour Tm |
| `max_tm_diff` | 3 °C | Balanced annealing of the pair |
| `gc_min`, `gc_max` | 15–55 % | Converted primers are AT-rich |
| `max_poly_x` | 8 | T runs are unavoidable after conversion |
| `max_poly_non_t` | 5 | Limits other homopolymers |
| `min_converted_c` | 4 | Discriminates converted from unconverted DNA |
| `three_prime_window` | 3 nt | A converted C near the 3′ end blocks extension on unconverted template |
| `max_cpg_in_primer` | 0 | Primer binding must not depend on methylation |
| `min_cpg_in_amplicon` | 3 | Enough CpGs to read out |
| `candidates_per_side` | 200 | Search breadth |
| `num_return` | 10 | Pairs reported per target |
| `max_primer_reuse` | 2 | Diversity of reported pairs |
| `fallback_targets` | 3 | Non-island CpG regions designed when no island yields a pair |

## Specificity check

Prepare the converted genome once with Bismark:

```bash
bismark_genome_preparation --bowtie2 /path/to/genome_folder   # folder containing the FASTA
```

Then set **Bismark genome folder** in the GUI, or call `check_bsp_specificity(design, "/path/to/genome_folder")`.

- **Alignment.** Each primer is aligned with Bowtie2 to both converted strands (Bismark's CT and GA indexes), with no gaps and up to `max_mismatch` mismatches.
- **Priming hits.** A hit counts only if none of its mismatches falls in the 3′-terminal three bases.
- **Predicted products.** Any forward hit followed by a reverse hit on the same converted strand within `max_product` bp is a predicted product. Products formed by a single primer are included.
- **Pairs checked.** The best `top_n` pairs of each target are checked.
- **Output.** The pair table gains `n_hits_left`, `n_hits_right`, `n_products`, `on_target`, `n_offtarget_products` and `specific`.

## Checking published primers from R

The GUI designs primers itself. To check primers designed elsewhere (MethPrimer, Primer3, a publication), map their positions on the extracted sequence back to the genome from R:

```r
primers <- parse_primer_output("
Left primer   1074  22  59.36  54.55  4  TTATTTAGTTTTTGGAGGTTAG
Right primer  1299  28  54.70  50.00 10  AAAACTAAAACCTACAAAAACTAAACTA
Product size: 226")
map_primer_coordinates(primers, promoter, sequence_variant = "bisulfite_requested")
```

- **Right-primer convention.** Primer3 and MethPrimer report a right primer's `Start` as its **5′ end**, which is its highest position. The right primer above therefore covers positions 1272–1299, giving the reported 226 bp product. Use `right_start = "leftmost"` only for tools that report the lowest position.
- **Product-size check.** When a `Product size` line is present, the mapped amplicon must match it.
- **Sequence check.** When primer sequences are present, each is checked against the converted template at its mapped position. The `sequence_match` column reports `converted`, `converted (methylated CpG)`, `unconverted` or `mismatch`.
- **Sequence variant.** `sequence_variant` names the sequence the primers were designed on (requested, complementary, or either bisulfite variant). Positions on the complementary sequences are flipped back before mapping.

## Coordinates and strands

- Genomic and extracted-sequence coordinates are 1-based and inclusive. Imported primer `Start` values are 1-based by default; use `start_convention = "zero-based"` for tools that report zero-based starts.
- The extracted sequence is in transcription orientation:
  - For `+` strand genes, position 1 is `region_start`.
  - For `-` strand genes, the sequence is reverse-complemented and position 1 is `region_end`.
  - Primer intervals are reported in ascending genomic order.
- **TSS-relative positions** follow the promoter convention: the TSS is `+1`, the base upstream is `-1`, and there is no position 0.
- **TSS definition.** The TSS is the transcript start on `+` and the transcript end on `-`. Transcripts sharing the same gene, scaffold, strand and TSS are collapsed into one selectable TSS, keeping all transcript IDs.
- **Bisulfite strands.** The two strands are no longer complementary after conversion, so each assay targets one of them:
  - The *sense* template is the gene's own strand; *antisense* is the opposite strand.
  - Results also report the Bismark label of the genomic strand: `OT` (original top, `+`) or `OB` (original bottom, `-`).
- **FASTA conversion.** The four FASTA downloads convert every `C` to `T`. `bisulfite_convert(x, cpg = "methylated")` keeps CpG cytosines.

## Tests

Run from the package root:

```r
install.packages("testthat")
testthat::test_local()
```

The tests use two reference fixtures:

- `inst/extdata/synthetic.*`: a minimal GTF/FASTA used for the TSS and coordinate tests;
- `tests/testthat/fixtures/bsp_reference.*`: two synthetic promoters with CpG islands, generated by `data-raw/make_bsp_fixture.R`, used only by the design tests and not installed with the tool.

Tests that call Primer3, Bowtie2 or Bismark are skipped when those tools are not installed. CI installs all three.
