test_that("Boulder-IO records are written and parsed", {
  expect_equal(boulder_record(list(SEQUENCE_ID = "a", PRIMER_PRODUCT_SIZE_RANGE = "100-200")),
               c("SEQUENCE_ID=a", "PRIMER_PRODUCT_SIZE_RANGE=100-200", "="))
  records <- parse_boulder(readLines(test_path("fixtures", "primer3_output.txt")))
  expect_length(records, 2L)
  expect_equal(records[[1L]][["PRIMER_PAIR_0_PRODUCT_SIZE"]], "130")
  expect_equal(records[[2L]][["PRIMER_LEFT_NUM_RETURNED"]], "0")
})

test_that("Primer3 right-primer positions are converted from the 5' end", {
  record <- parse_boulder(readLines(test_path("fixtures", "primer3_output.txt")))[[1L]]
  left <- primer3_oligo(record, "LEFT")
  right <- primer3_oligo(record, "RIGHT")
  expect_equal(c(left$start, left$end), c(1L, 26L))
  expect_equal(c(right$start, right$end), c(109L, 130L))
  expect_equal(right$end - left$start + 1L, as.integer(record[["PRIMER_PAIR_0_PRODUCT_SIZE"]]))
  expect_equal(right$tm, 59.765)
  template <- record[["SEQUENCE_TEMPLATE"]]
  expect_equal(reverse_complement_text(substr(template, right$start, right$end)), right$sequence)
})

test_that("Primer3 tags follow the BSP parameters", {
  params <- utils::modifyList(bsp_primer_params(), list(tm_opt = 60, product_size_max = 400L))
  tags <- primer3_global_tags(params)
  expect_equal(tags$PRIMER_OPT_TM, 60)
  expect_equal(tags$PRIMER_PRODUCT_SIZE_RANGE, "200-400")
  expect_equal(tags$PRIMER_FIRST_BASE_INDEX, 1L)
  expect_equal(tags$PRIMER_MAX_NS_ACCEPTED, 0L)
})

test_that("a missing primer3_core gives install hints", {
  old <- options(promoterPrimerMapper.primer3 = file.path(tempdir(), "no-such-primer3"))
  on.exit(options(old), add = TRUE)
  expect_error(primer3_path(), "was not found")
})

test_that("primer3_core scores oligos and pairs", {
  skip_if_no_primer3()
  params <- bsp_primer_params()
  scored <- primer3_check_oligos(c("GGGTATGTTGTTAGTTAGTGTGTTAG", "GGGGGGGGGGGGGGGGGGGGGG", "ATATATTAATATATATTTATAT"), "LEFT", params)
  expect_equal(scored$candidate, 1L)
  expect_true(scored$tm > 50 && scored$tm < 65)
  expect_equal(attr(scored, "rejected"), c("GC out of range" = 2L))
  record <- parse_boulder(readLines(test_path("fixtures", "primer3_output.txt")))[[1L]]
  pairs <- primer3_check_pairs(record[["SEQUENCE_TEMPLATE"]], record[["SEQUENCE_PRIMER"]], record[["SEQUENCE_PRIMER_REVCOMP"]],
                               utils::modifyList(params, list(product_size_min = 100L)))
  expect_equal(nrow(pairs), 1L)
  expect_true(is.finite(pairs$compl_any))
})

test_that("Primer3 explain strings become readable rejection reasons", {
  expect_equal(primer3_reject_reason("considered 1, GC content failed 1, ok 0"), "GC out of range")
  expect_equal(primer3_reject_reason("considered 1, low tm 1, ok 0"), "Tm too low")
  expect_equal(primer3_reject_reason("considered 1, high hairpin stability 1, ok 0"), "hairpin")
  expect_equal(primer3_reject_reason("considered 1, something new 1, ok 0"), "something new")
  expect_equal(primer3_reject_reason(NA_character_), "rejected by Primer3")
})

test_that("a primer3_core crash only drops the records that cause it", {
  skip_on_os("windows")
  fake <- tempfile("primer3")
  writeLines(c("#!/bin/sh",
               "input=$(cat)",
               "case \"$input\" in *CRASH*) kill -ABRT $$ ;; esac",
               "printf '%s\\n' \"$input\" | grep '^SEQUENCE_ID=' | while read -r id; do printf '%s\\nPRIMER_LEFT_NUM_RETURNED=0\\n=\\n' \"$id\"; done"), fake)
  Sys.chmod(fake, "0755")
  records <- list(list(SEQUENCE_ID = 1L, SEQUENCE_PRIMER = "ACGT"), list(SEQUENCE_ID = 2L, SEQUENCE_PRIMER = "CRASH"),
                  list(SEQUENCE_ID = 3L, SEQUENCE_PRIMER = "ACGT"))
  expect_error(run_primer3(records, fake), "exited with status")
  output <- run_primer3_tolerant(records, fake)
  expect_equal(vapply(output, primer3_crashed, logical(1)), c(FALSE, TRUE, FALSE))
  expect_equal(vapply(output, record_value, character(1), key = "SEQUENCE_ID"), c("1", "2", "3"))
  expect_error(run_primer3_tolerant(records[2L], fake), "exited with status")
})
