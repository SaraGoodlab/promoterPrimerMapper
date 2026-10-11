# GUI layout and design review

Reviewed: October 14, 2026

Baseline: `eb58b3a` — Clarify local and HPC run instructions

Implementation: applied in three passes on October 14, 2026. The findings below describe the baseline; see the implementation record at the end for the resulting behavior and validation.

## Scope and main recommendation

This review covers the current Shiny interface, its generated HTML, and selected server behaviors. Sources reviewed: [R/app_ui.R](R/app_ui.R), [R/app_server.R](R/app_server.R), [R/bsp_view.R](R/bsp_view.R), [R/bsp_design.R](R/bsp_design.R), and [tests/testthat/test-app.R](tests/testthat/test-app.R). Behavioral findings were checked with `shiny::testServer()` against the working source. No browser screenshots or measurements of rendered overflow were collected; viewport, contrast, and keyboard recommendations require browser validation.

Keep the promoter map central, bring target selection and pair review closer together, and make the relationship between current inputs and completed results explicit. Improve result provenance and empty states before changing the theme. The current interface already has useful scientific content; the main opportunity is to reduce scrolling, simplify comparison, and prevent users from misreading previous results as a new design.

## Current layout

```text
Promoter Primer Mapper
┌──────────────────────┬─────────────────────────────────────────────┐
│ Setup sidebar: 25%   │ Main content: 75%, vertically stacked       │
│                      │                                             │
│ 1. Reference files   │ Status message                              │
│ Annotation / FASTA   │                                             │
│ FASTA index          │ Promoter CpG map                            │
│ External tools ▸     │ Summary, 280 px map, hover text              │
│                      │ Range inputs + seven-column selection table │
│ 2. Gene and TSS      │ Eight-column CpG island table               │
│ Find gene / TSS      │                                             │
│                      │ BSP primer design                           │
│ 3. Extraction        │ Strand, target mode, target preview          │
│ Upstream/downstream  │ Five closed parameter groups; 21 parameters │
│ Extract sequence     │ Reset / Design BSP primers                  │
│                      │                                             │
│                      │ Primer pairs                                │
│                      │ Target summary, 12–14-column pair table      │
│                      │ Hints, diagnostics, selected-pair dropdown   │
│                      │ 12-column primer detail table, sequence      │
│                      │ Four downloads, specificity controls         │
│                      │                                             │
│                      │ Extracted sequence ▸ + four FASTA downloads  │
└──────────────────────┴─────────────────────────────────────────────┘
```

The setup is numbered, but target selection, design, specificity, and export are not presented as equally clear stages. After setup, the sidebar continues to occupy one quarter of the desktop grid while the most demanding comparison views share the remaining width.

### What to preserve

- The continuous promoter-centered workspace, with access to both map and results.
- TSS-relative and genomic coordinates, including correct minus-strand orientation and the no-zero convention.
- Explicit target windows, strand choice, and scientific parameter defaults shared with the R API.
- Collapsible advanced parameters and extracted-sequence content.
- Existing progress feedback and actionable hints for targets without pairs.
- Separate pair-table, oligo-order, and amplicon-reference exports.

## Priorities

P0 means address before relying on the interface for another design; P1 means address in the first usability pass; P2 means subsequent refinement. Effort estimates are relative, not delivery estimates.

| Priority | Improvement | Why it matters | Effort |
|---|---|---|---|
| P0 | Show completed-design provenance and changed-input state | Previous results can remain visible after a new gene search or failed redesign. | Medium |
| P1 | Replace initial pair-output errors with guided empty states | A fresh session currently raises errors in pair detail and sequence outputs. | Small |
| P1 | Bring target controls and selected-pair review beside the map | The user must repeatedly move through several large sections. | Medium |
| P1 | Reduce default table width and connect row selection to details | Long oligos and technical columns obstruct pair comparison. | Medium |
| P1 | Make range selection responsive and describe clicks accurately | Brush updates wait five seconds; island-click help implies immediate design. | Small |
| P1 | Keep validation and status close to the affected control | Errors share an informational banner; invalid settings silently use defaults in previews. | Medium |
| P1 | Label specificity outcomes as checked, unchecked, or inconclusive | A Boolean or blank value does not explain what was assessed. | Medium |
| P1 | Clarify local, HPC, and container file paths | “This computer” can be mistaken for the browser user's computer. | Small |
| P2 | Add setup readiness indicators | New users currently need external instructions to get started. | Small–medium |
| P2 | Focus sequence review on the chosen amplicon | The current selected-pair view renders the entire extracted template. | Medium |
| P2 | Improve responsive layout, keyboard access, and visual consistency | Dense grids and plot-only interactions need accessible alternatives. | Medium |
| P2 | Group exports by scope | Eight download buttons use long, inconsistent labels. | Small |

## Detailed suggestions

### 1. Make completed results distinguishable from current settings

**Evidence:** `app_server.R:99` updates gene/TSS choices without clearing or labeling the existing extraction and design. `run_design()` at line 134 replaces the design only on success. `display_params()` at line 130 substitutes defaults when parameter validation fails. Once a design exists, the map uses its targets rather than the current planned targets (`app_server.R:238`; `bsp_view.R:25`).

Behavioral probes confirmed that a successful search for `BSPMINUS` leaves a previous `BSPPLUS` extraction and design in state, and an invalid-Tm redesign leaves the previous design intact. Preserving useful results is reasonable, but the interface must identify them as belonging to the previous inputs.

**Suggestions:**

- Add a persistent context strip: gene, selected/extracted TSS, scaffold, genomic strand, extracted interval, and completed-design target mode.
- Store a snapshot of all inputs that determine a completed design: reference identity, TSS, extraction bounds, template strand, targets, and parameter values.
- When those inputs change, show “Settings changed — results below are from the previous design” and a clear action to extract or redesign as appropriate.
- A failed redesign should retain a visibly labeled previous result or clear the result; it must not appear to have succeeded.
- Keep proposed targets and completed-design targets visually distinct. Label each layer and avoid showing a new target preview as though existing pairs cover it.
- If previous results remain exportable, label exports with their completed-design context and show that context beside the export controls.

**Acceptance:** Searching for another gene, changing a target, editing parameters, and submitting an invalid redesign all leave an unambiguous indication of which inputs produced displayed and downloadable results.

### 2. Give each stage an intentional empty and readiness state

**Evidence:** All four main cards are declared immediately (`app_ui.R:70`). Only the map summary gives an initial instruction. The pair-detail and sequence outputs call helpers that reject `NULL` designs before the selected-rank requirement can suppress rendering (`app_server.R:374`, line 386; `bsp_design.R:737`). Both initial outputs raised “design must be a bsp_design object from design_bsp_primers()” in `testServer()`.

**Suggestions:**

- Before extraction, show one short welcome panel: “Choose reference files and find a gene.”
- Gate pair detail, sequence review, and downloads on a valid design with pairs before calling their helpers.
- Replace missing-data areas with stage-specific guidance: “Extract a promoter to choose a target,” “Design primers to compare pairs,” and “Select a pair to review oligos.”
- Disable unavailable actions and place the reason beside the action. Preserve server-side checks as well.
- Distinguish missing prerequisites, a valid design with zero pairs, and a computation failure.
- Guide users through supplying their reference files, finding a gene, and extracting its promoter.

**Acceptance:** A fresh session and a zero-pair design contain useful guidance and no raw helper errors. Reference setup is clear without opening the README.

### 3. Rebalance the page around selection and comparison

**Evidence:** The current 3/9 Bootstrap grid and four vertically stacked cards are defined in `app_ui.R:42`. Target editing lives below the map; selected-pair details live below the full result table and diagnostics.

**Suggestions:**

- Use a collapsible setup sidebar; provide a compact reference/gene summary after extraction. Keep setup readily reachable when changing genes.
- Put target mode, selected island/range, template strand, and the design action in a compact control area beside or immediately below the map.
- Keep the result list and selected-pair inspector adjacent on wider screens. Put diagnostics below these views or inside an expandable section when a target fails.
- Maintain a continuous workspace rather than requiring navigation away from the map to inspect a pair.
- On narrower screens, stack map, target controls, results, and details in task order. Avoid sticky panels that consume most of the viewport at high zoom.

Suggested desktop arrangement:

```text
Header: app name | gene / TSS / region | design status
┌──────────────────┬─────────────────────────────────────────────────┐
│ Collapsible      │ Promoter map                                   │
│ setup            │ Genomic + TSS-relative axes, target/pair legend │
│                  ├───────────────────────────┬─────────────────────┤
│ Reference files  │ Target + strand controls  │ Target summary      │
│ Gene / TSS       │ Parameters ▸   [Design]   │ Validation / status │
│ Extraction       ├───────────────────────────┴─────────────────────┤
│                  │ Pair results + selected-pair inspector          │
│ Tool settings ▸  │ Filters/list        | oligos, metrics, products │
│                  ├─────────────────────────────────────────────────┤
│                  │ Specificity status + check action               │
│                  │ Exports     Diagnostics ▸     Full sequence ▸   │
└──────────────────┴─────────────────────────────────────────────────┘
```

**Acceptance:** At a typical laptop viewport, selecting a target and reviewing a pair require substantially less scrolling; setup can be collapsed without losing gene/TSS context.

### 4. Make tables support comparison instead of displaying every field

**Evidence:** `bsp_table` exposes 12 columns before specificity and up to 14 afterward (`app_server.R:350`), including both complete oligo sequences. The selection table has seven columns, the island table eight, diagnostics eleven, and primer detail twelve. These use ordinary `renderTable()` outputs without explicit page-level table overflow styling or row selection.

**Suggestions:**

- Default pair columns: Pair, Target, Template strand, Product (bp), CpGs, Forward Tm, Reverse Tm, Score, Specificity.
- Move full sequences, assay labels, coordinates, and thermodynamic details into the selected-pair inspector. Retain all scientific fields in exports.
- Add target/strand filters, sortable metrics, and pagination or a bounded scroll area for larger designs.
- Make a row selection update the same selected-rank state used by the map, inspector, and selected-pair exports. Keep a keyboard-operable selection control.
- Show the selected row clearly. Use stable pair rank as identity even when sorting or filtering changes row position.
- Use monospace for sequences, aligned numeric columns, consistent decimal precision, and units in display headings. Explain how the score is constructed rather than implying a calibrated quality probability.
- Replace the two-row, twelve-column primer detail table with two stacked oligo cards and a compact metrics grid.
- Convert the single-row selection table into a short definition list or metric strip. Put detailed island statistics in a collapsible table or detail view.

**Acceptance:** Primary pair comparison is readable at laptop width; long sequences and diagnostics do not force the whole page to scroll horizontally.

### 5. Tighten map interactions and target feedback

**Evidence:** The brush debounce is `5000` ms (`app_ui.R:75`). Help says clicking an island designs primers, but `app_server.R:204` only selects the island and requests a later design action. Numeric From/To controls and the island dropdown already offer useful alternatives to plot interaction.

**Suggestions:**

- Reduce brush delay to approximately 150–300 ms, then verify that selection feedback and redraw do not disrupt dragging.
- Change help to “Drag to select a range. Click an island to select it, then choose Design primers.”
- Keep the selected range, active target, and selected pair as separate concepts with distinct labels and visible boundaries.
- Place “Use range as target” beside the range inputs, with an immediate target summary after activation.
- Add a clear-selection action. Show extracted bounds and the no-zero rule near the numeric range controls.
- Make island selection available from the island list/table as well as the existing dropdown.
- Keep target, CpG, primer, amplicon, and TSS styles consistent between map and sequence review. Use boundaries, labels, or patterns as well as color.
- Offer zoom-to-target/amplicon and reset-to-promoter controls if wide extraction windows make loci difficult to inspect. Treat this as a later enhancement, not a prerequisite for the first layout pass.

**Acceptance:** Selection feedback appears promptly, help matches actual behavior, and every target-selection task can be completed without dragging or hovering.

### 6. Clarify computation, validation, and specificity states

**Evidence:** The shared status output always uses `alert-info` (`app_server.R:83`), including caught failures. Parameter previews silently fall back to defaults when invalid. Extraction immediately runs automatic design (`app_server.R:176`). Specificity is configured at the bottom of the pair card while its required genome folder is hidden in sidebar tool settings.

**Suggestions:**

- Use distinct states for success, warning, failure, and in-progress work, with status text announced to assistive technology.
- Place validation near the field or group responsible, and open an accordion containing the invalid field. Do not silently display a default target preview for invalid parameters.
- Rename the extraction action to “Extract and design” if automatic design remains the default, or offer a clearly labeled automatic-design toggle.
- Disable duplicate actions while a computation runs and show the current stage using the existing progress callbacks. Report elapsed time where useful; only add cancellation if the execution path supports it.
- Put specificity readiness and its genome-folder input or settings shortcut beside the check action.
- Display “Not checked,” “Specific within checked criteria,” “Off-target products found,” and “On-target product not found” instead of exposing only `TRUE`/`FALSE`/blank.
- Distinguish unchecked pairs from a failed check. Show how many pairs per target were checked and the mismatch/product limits used for completed checks; flag changes to those settings until rerun.
- Preserve the selected rank when the specificity columns update and show outcomes in both the table and pair inspector.

**Acceptance:** Users can identify the affected input, know whether an operation is running or failed, and understand the scope of each specificity result without interpreting raw Boolean fields.

### 7. Make setup accurate for local, HPC, and container use

**Evidence:** The file help says “Files are opened directly from this computer” (`app_ui.R:48`), although the app may execute on a remote compute node or in a container. The annotation label only mentions GTF even though GFF/GFF3 are supported. All paths are manually entered.

**Suggestions:**

- Use “Annotation path (GTF/GFF3)” and “Files are read by the machine running the app. In a container, these paths must be bind-mounted.”
- Add lightweight readiness indicators for annotation, FASTA/index, Primer3, and optional specificity tools/index. Readiness checks should not re-import a large annotation on each keystroke.
- Suggest the `.fai` path from the FASTA path while allowing a different existing index.
- Give targeted messages for missing files, unreadable files, missing indexes, and tools unavailable on PATH.
- Keep direct paths as the standard HPC workflow. If a file picker is added, identify whether it browses server files or uploads browser-local files; large genome uploads should not be implied by a generic Browse label.
- Offer concise TSS option labels with transcript counts, and show full transcript IDs in an adjacent expandable detail view.

**Acceptance:** A user accessing the app through an HPC tunnel knows where files must exist and can determine readiness without trial-and-error design attempts.

### 8. Focus sequence and export review on the chosen pair

**Evidence:** `sequence_view_html()` renders the complete template in 60-base lines (`bsp_view.R:166`), even when only one amplicon is selected. Four pair exports and four extracted-sequence exports are separate buttons with long labels (`app_ui.R:117`, line 138).

**Suggestions:**

- Default sequence review to the selected amplicon with short flanks; provide “Show full promoter.”
- Label template orientation, assay strand, and the methylated-CpG display assumption prominently. Preserve the genomic and TSS-relative meaning when viewing an antisense template.
- Give forward and reverse oligos copy actions with a clear copied state and 5′→3′ labels.
- Group exports by scope: “All pairs” with TSV/CSV choice; “Selected pair” with order sheet and amplicon references; “Extracted promoter” with variant choice.
- Keep gene, target, rank, and design provenance visible at export time. Explain the amplicon-reference variants in brief help.

**Acceptance:** A selected pair can be reviewed and exported without scrolling through the complete promoter or interpreting eight separate long button labels.

### 9. Apply a consistent visual and accessibility system

**Evidence:** The interface uses Flatly, nested Bootstrap grids, a fixed-height plot, several small-text help elements, and custom scientific colors (`app_ui.R:25`, line 39). Numerical controls use fixed three-column group layouts. The map is a base-graphics image; source review does not establish actual contrast or keyboard accessibility.

**Suggestions:**

- Use a small spacing scale, consistent card headings, and one primary action per active stage. Keep secondary actions such as reset and export visually subordinate.
- Make gene/region context and pair metrics easier to scan through hierarchy rather than extra decoration.
- Arrange parameter fields in three columns only when space permits, two on medium widths, and one on narrow widths; keep related min/optimal/max values together.
- Retain numeric range entry and explicit island selection as keyboard alternatives to the plot. Provide a textual summary of the current map selection and pair.
- Verify contrast for muted help, green island labels, red CpGs, and yellow primer highlights. Ensure status and selection remain distinguishable without color.
- Ensure visible focus indicators, associated labels/help, accessible table headings, and appropriate live-status announcements.
- Verify that page content reflows at 200% zoom; horizontal scrolling should be limited to deliberate sequence/table viewports.

## Suggested implementation order

1. **Result integrity and empty states:** add completed-input snapshots, changed-input banners, output guards, explicit readiness, and local validation. Make clicks/help consistent and reduce brush delay.
2. **Layout and comparison:** collapse completed setup, move target controls next to the map, introduce compact pair comparison with a synchronized inspector, and collapse diagnostics/full sequence.
3. **Finishing pass:** group exports, add setup checks, refine specificity presentation, and validate responsive behavior and accessibility in a browser.

Prefer small changes to the existing components first. A theme replacement or a plot-library migration is not necessary to achieve the highest-priority improvements. Any later table or layout component should be evaluated for the actual Shiny/bslib versions supported by the package.

## Verification plan for a future implementation

| Scenario | Expected outcome |
|---|---|
| Fresh launch, missing files or Primer3 | Guided empty states, actionable readiness feedback, no raw pair-helper error. |
| Bundled plus- and minus-strand examples | Context, axes, target selection, oligo orientation, and pair details agree. |
| New gene/TSS/reference after a completed design | Previous results are cleared or explicitly identified; exports preserve their original provenance. |
| Parameter/target edits, including failed redesign | Pending changes are visible; previous results cannot be mistaken for a successful rerun. |
| Valid design with zero pairs or fallback targets | No-pair/fallback explanation is visible and diagnostics are accessible. |
| Brush, typed range, dropdown/table island selection | Consistent targets and prompt feedback; all operations have keyboard alternatives. |
| Sorting/filtering pair rows | Selected rank, map, inspector, and exports refer to the same stable pair. |
| Partial, failed, and completed specificity checks | Unchecked pairs and check outcomes are distinct; settings and checked scope are shown. |
| Desktop 1440×900, laptop 1280×800, tablet 768×1024, narrow 390×844 | Controls remain reachable; layout stacks coherently; wide content scrolls within its own container. |
| Keyboard-only operation and 200% browser zoom | Focus is visible, all essential tasks remain possible, and content reflows. |
| HPC/container session accessed through a tunnel | File guidance describes the execution host and bound paths accurately. |

Extend the existing interaction tests where behavior changes, especially initial outputs, result provenance, validation, and stable pair selection. Update layout-specific tests intentionally: the current suite explicitly expects a single screen without navigation tabs. Browser checks remain necessary for actual width, focus behavior, and visual contrast; server tests cannot establish those properties.

## Implementation record

### Pass 1: result integrity and empty states

- Added extraction and completed-design input snapshots, including reference paths, sizes and modification times. Lightweight file checks also detect on-disk reference changes.
- Added a persistent gene/TSS/interval context strip and warnings for changed inputs or failed redesigns. Previous results remain explicitly labeled and exportable.
- Blocked reuse of an old extraction for a changed reference/gene and removed silent parameter-default substitution from previews.
- Guarded empty pair-detail, sequence and download paths; added guided initial and zero-pair states.
- Added severity-specific, live status messages and parameter validation that opens the affected parameter group.

### Pass 2: workspace and comparison

- Replaced the fixed 25% sidebar with a collapsible, resizable setup panel that closes after successful extraction.
- Moved target, range, strand and design controls into a responsive grid directly below the promoter map.
- Introduced a compact, bounded pair table with target/strand filters, sort controls, stable-rank row selection and an adjacent oligo inspector on wide screens.
- Added oligo copy controls, compact metrics and expandable thermodynamic details. Full scientific data remain in exports.
- Reduced brush debounce to 200 ms, corrected island-click help, and added clear-selection and map-zoom controls.
- Made sequence review default to the selected amplicon with 20 bp flanks, with a full-promoter toggle. Coordinate conversions are preserved on both gene and template strands.
- Made coordinate tick density and legend layout adapt to plot width.

### Pass 3: readiness, specificity, export and accessibility

- Added reference/tool readiness, automatic `.fai` suggestions that preserve a custom index, and execution-host/container file guidance.
- Disabled actions without prerequisites and during browser-visible computation; retained server-side checks.
- Moved specificity into its own stage with genome/tool settings, completed criteria, explicit outcome labels and warnings for edited criteria or failed checks. Pair selection survives specificity updates.
- Grouped exports by scope and displayed their saved-result context. GUI pair-table exports include completed source paths/time and specificity criteria when available.
- Added responsive parameter/result grids, table-scoped overflow, accessible row buttons, keyboard alternatives, focus styling and live copy/status feedback.
- Adjusted status and secondary-action colors for contrast; the rendered success message measured approximately 6.93:1 text/background contrast.
- Updated the README for the revised GUI. Declared the bslib/htmltools versions used by the new layout and excluded this review document from R source builds.

### Validation performed

- Added regression coverage for initial/zero-pair states, invalid previews, changed reference/gene/design settings, failed reruns, stable rank selection across sorting/filtering, specificity outcomes and reruns, sequence cropping on both strands, map zoom/clear-selection, custom index preservation and saved-result export provenance.
- The complete test suite passed in the existing container with real Primer3, Bowtie2 and Bismark. The host also passed its runnable GUI regressions, skipping the real-Primer3 integration test because that executable is absent on the host.
- A built-source `R CMD check --no-manual` completed with `Status: OK` on the host.
- Firefox checks exercised empty launch, reference setup, real design, row selection, invalid-parameter feedback, clear selection, pair zoom, full-sequence expansion and keyboard interaction. Invalid parameters opened the relevant accordion group.
- Checked desktop/laptop widths of 1440 and 1280 pixels, a 768-pixel tablet view and the actual application in a 390-pixel frame. Each had no page-level horizontal overflow or visible Shiny output errors in the exercised states.
- Verified that activating a pair button with Enter updates the selected rank and inspector. Native 200% Firefox zoom produced a 640-pixel viewport at laptop width, with no page-level horizontal overflow or visible Shiny output errors.

### Follow-up: quiet polling and conditional range controls

- Replaced the periodic reactive invalidation with a cancellable background file check. File fingerprints are published only when paths, sizes, modification times, readability, executability, or Bismark index parts change. Checks stop when their session closes.
- Send action-readiness messages only when their enabled states change. Working feedback and button locking are driven by explicit computations rather than generic Shiny busy events.
- Show the From/To inputs and associated range controls only for **Selected range**. Switching target modes preserves their values.
- Regression checks cover unchanged polling, on-disk reference/index changes, executable permission changes, operation cleanup and closed sessions. A real-browser observation across 13 idle seconds recorded zero busy events, button toggles, output refreshes or working messages. Actual input changes continued to update warnings and readiness.
- Browser checks confirmed that the range fields hide in all-island and selected-island modes, show in selected-range mode, and retain edited bounds when switching modes.

### Earlier iteration: Antique color design (superseded below)

- Centralized the exact `rcartocolor::carto_pal(11, "Antique")` swatches in `R/app_theme.R`, with a `#F7F3E3` page background and warm card/sidebar surfaces. The swatches are stored locally, so no additional runtime package is required.
- Applied mauve primary actions, slate selection/focus accents, sage islands, terracotta CpGs, and gold primers consistently across Bootstrap controls, statuses, map layers, and sequence highlights.
- Used darker palette shades for small text, including CpGs overlapping primer highlights. Representative rendered text combinations measured at least 4.98:1 contrast in Firefox, accounting for translucent backgrounds.
- Existing GUI and plot tests passed. Browser checks covered empty launch, a real Primer3 design, selected-pair review, invalid-parameter alerts, and a 768-pixel view without page-level horizontal overflow or Shiny output errors.

### Follow-up: neutral surfaces, plum accent, and technical typography

- Replaced the cream/Antique theme with `#F7F7F8` page and status backgrounds, white panels, `#E3E3E8` borders, dark text, and muted helper text. The white header has a 2px plum bottom border; section headers are neutral and use weight 600.
- Reserved `#6B3F5A` plum for the primary design action, selected states, focus, and information markers; its hover color is `#57324A`. Routine buttons have neutral outlines. Status banners use a 3px left border, with amber/red reserved for warnings/errors. The completed-design chip uses `#F3ECF0` with plum text.
- Bundled IBM Plex Sans and Mono at weights 400/600 with the upstream license, so fonts load locally without a CDN. Added tabular UI numerals and monospaced chart ticks. Chart labels are at least 11px, the map helper is 13px, and body text is 16px.
- Changed CpGs to thin, 70%-opaque slate marks, islands to teal, and primers to amber. Target windows share 9%-opaque plum shading with dashed outlines; range and amplicon outlines are solid. The active island has an accent tint and plum outline. Larger chart labels have additional vertical room and the pair caption stays inside the plotted interval.
- Consolidated map instructions, enlarged radio targets to 44px with clickable labels and spacing, and used 2px plum keyboard focus rings. Applied the same feature colors to the sequence view, preserving readable text where highlights overlap.
- Firefox verified local loading of all four font files, a real Primer3 design, selected-pair changes, label clicks, conditional range fields, the exact primary hover color, keyboard focus, and neutral status/error presentation. Representative rendered text combinations measured at least 6.24:1 contrast; 1440px, 768px, and 390px views had no page-level horizontal overflow or Shiny output errors. The narrow map has extra space for its three-row legend.

### Follow-up: primer CSV downloads

- Both all-primers and selected-pair downloads now use a 14-column CSV with one row per oligo. Every metadata line begins with `#`; the column-name row and data are ordinary CSV with correct quoting.
- Metadata identifies the design, parameters, included pairs, template strands, coordinate convention, and units. Start/Stop are 1-based inclusive genomic bounds; forward/reverse names identify each pair's two primers.
- Retained Primer3 product Tm in the design and requested its calculation with a nonbinding Tm bound. Exported self-complementarity fields use the existing Primer3 thermodynamic values, with units explained in the header. CpGs in Product counts complete original CpG dinucleotides across the full amplicon, including primer-binding sites; it excludes a CpG extending past the product boundary.
- The detailed pair table remains available through the R API's `write_bsp_pairs()`. GUI primer downloads and `write_order_sheet()` default to the requested CSV schema.
- Regression checks cover metadata line prefixes, exact columns, multiline metadata, quoted fields, metrics, pair scope/order, legacy missing values, and actual GUI download handlers. The full container suite passed, and a real minus-strand design produced an independently parsed 14-column CSV with populated product Tm values.

Additional enhancements such as uploading large references, migrating the plot library, named biological parameter presets, and process cancellation remain outside these implemented passes. The interface retains server-side paths, the existing base-graphics plot and the established scientific defaults.
