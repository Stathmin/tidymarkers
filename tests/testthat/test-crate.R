# J8: GWAS loci + dossier -> marker crate.  The soybean-named fixture is a small
# synthetic region, not Wm82: it exists to prove the crate contract and its
# determinism on a committed input.

test_that("a soybean locus crate is deterministic and carries full provenance", {
  skip_if_no_oligoarrayaux()
  created <- "2026-01-01T00:00:00Z"
  genome <- read_genome_fasta(scenario_file("j8_soy_genome.fa"))
  gwas <- scenario_file("j8_soy_gwas-loci.tsv")
  dossier <- scenario_file("j8_soy_dossier.tsv")
  variants <- scenario_file("j8_soy_variants.tsv")
  haplotypes <- scenario_file("j8_soy_haplotypes.csv")

  build <- function(dir) {
    marker_locus_crate(
      gwas, dossier, genome,
      variants = variants,
      haplotypes = haplotypes,
      dir = dir,
      created = created
    )
  }
  a <- build(withr::local_tempdir("crate-a-"))
  b <- build(withr::local_tempdir("crate-b-"))

  files <- function(dir) {
    out <- list.files(dir, recursive = TRUE)
    sort(out[!grepl("(^|/)solved\\.rds$", out)])
  }
  expect_identical(files(a), files(b))
  for (file in files(a)) {
    expect_identical(
      readBin(file.path(a, file), "raw", n = file.size(file.path(a, file))),
      readBin(file.path(b, file), "raw", n = file.size(file.path(b, file))),
      info = file
    )
  }
  expect_false(file.exists(file.path(a, "solved.rds")))
  expect_true(file.exists(file.path(a, "artefact", "solved.rds")))

  record <- jsonlite::fromJSON(file.path(a, "marker-crate.json"), simplifyVector = FALSE)
  expect_identical(record$schema, "tidymarkers/crate@1")
  expect_identical(record$contract$version, "tidymarkers/crate-contract@1")
  expect_identical(record$gwas$runs, list("GWAS-2026-11"))
  expect_identical(record$gwas$genome_build, "Wm82.a4.v1")
  expect_identical(record$dossier$sources, list("Merovingen"))
  expect_equal(record$dossier$confidence$min, 0.65)
  expect_equal(record$dossier$confidence$max, 0.80)
  expect_identical(record$selection$selected, list("soy_oc_locus"))

  candidates <- vapply(
    record$selection$candidates,
    function(x) paste(x$pos, x$status, sep = ":"),
    character(1)
  )
  expect_setequal(candidates, c("300:not_selected", "100:selected", "300:duplicate"))
  selected <- Filter(function(x) x$status == "selected", record$selection$candidates)[[1L]]
  expect_identical(selected$pos, 100L)
  expect_equal(selected$ld_r2, 1)
  expect_identical(selected$maf_source, "table")

  expect_identical(record$panel$targets, 1L)
  expect_identical(record$panel$designed, 1L)
  expect_identical(record$panel$one_pot$groups, 1L)
  expect_identical(record$specificity$builds, list("Wm82.a4.v1"))
  expect_identical(record$specificity$rows, 3L)
  expect_identical(record$specificity$index, "memory")

  specificity <- utils::read.csv(
    file.path(a, "specificity.csv"),
    colClasses = "character"
  )
  expect_equal(nrow(specificity), 3L)
  expect_true(all(specificity$locus == "soy_oc_locus"))
  expect_true(all(specificity$build == "Wm82.a4.v1"))

  artefact <- jsonlite::fromJSON(file.path(a, "artefact.json"), simplifyVector = FALSE)
  roles <- vapply(artefact$inputs, function(x) x$role, character(1))
  expect_true(all(c("gwas", "dossier", "candidates") %in% roles))
  input <- function(role) artefact$inputs[[which(roles == role)[[1L]]]]
  expect_identical(input("gwas")$source, "gwas:GWAS-2026-11")
  expect_identical(input("dossier")$source, "dossier:Merovingen")
  expect_match(input("gwas")$sha256, "^[0-9a-f]{64}$")
  expect_match(input("dossier")$sha256, "^[0-9a-f]{64}$")
  expect_identical(input("candidates")$name, "2 candidate(s)")

  meta <- jsonlite::fromJSON(
    file.path(a, "ro-crate-metadata.json"),
    simplifyVector = FALSE
  )
  root <- Filter(
    function(x) identical(x[["@id"]], "./"),
    meta[["@graph"]]
  )[[1L]]
  expect_identical(root$identifier, record$panel$artefact_id)
  parts <- vapply(root$hasPart, function(x) x[["@id"]], character(1))
  expect_true(all(c(
    "gwas-loci.csv", "dossier.csv", "selection.csv", "specificity.csv",
    "marker-crate.json"
  ) %in% parts))

  back <- marker_import_rocrate(a)
  expect_identical(back$id, record$panel$artefact_id)
})

test_that("a crate refuses inputs its contract cannot type", {
  skip_if_no_oligoarrayaux()
  genome <- read_genome_fasta(scenario_file("j8_soy_genome.fa"))
  gwas <- data.frame(
    locus = "soy_oc_locus", chrom = "Gm01", pos = 300L,
    ref = "A", alt = "G", population = "soybean_panel",
    genome_build = "Wm82.a4.v1", gwas_run = "GWAS-2026-11",
    maf = 0.4
  )
  dossier <- data.frame(
    locus = "soy_oc_locus", gene = "Glyma.01G000100",
    source = "Merovingen", source_version = "2026-11",
    retrieved_at = "2026-11-28", confidence = 0.8
  )
  crate <- function(gwas_tbl = gwas, dossier_tbl = dossier, ...) {
    marker_locus_crate(
      gwas_tbl, dossier_tbl, genome,
      dir = withr::local_tempdir(),
      created = "2026-01-01T00:00:00Z",
      ...
    )
  }

  expect_error(crate(gwas_tbl = gwas[, c("locus", "chrom", "pos")]), "gwas")
  expect_error(crate(gwas_tbl = rbind(gwas, gwas)), "repeats locus")
  expect_error(
    crate(dossier_tbl = transform(dossier, locus = "other_locus")),
    "no row for GWAS locus"
  )
  expect_error(
    crate(dossier_tbl = rbind(dossier, transform(dossier, locus = "other_locus"))),
    "absent from the GWAS table"
  )
  expect_error(crate(maf_min = 0.99), "No candidate")
  expect_error(
    crate(
      variants = data.frame(
        locus = "soy_oc_locus", chrom = "Gm01", pos = 100L,
        ref = "C", alt = "G", population = "soybean_panel"
      ),
      ld_min = 0.5
    ),
    "LD is unmeasured"
  )
})
