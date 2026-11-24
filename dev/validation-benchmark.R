# =============================================================================
#  dev/validation-benchmark.R  --  a task set with an established verdict
# =============================================================================
#  Validation, not generation: every task here has a verdict that was fixed
#  before the package ran it (`expected`).  The package's own checks produce an
#  `observed` value, and the pair is compared.  A task without an independently
#  established verdict is marked `gap`, and a task that could not run is
#  `not_run` with its reason — neither is ever invented or reported as a pass.
#
#  Families (25 tasks):
#    fixture   2   hand-written panel, expectations live in test-panel.R
#    soybean  13   recorded KASP assays whose operator outcome is "worked"
#    measured  4   cross-tool / cross-paper numbers quoted in VALIDATION.md
#    contract  6   vignette prose against live code (the ZX3 contradictions)
#
#  verdict: pass | fail | not_run | gap
#
#  Run:  Rscript dev/validation-benchmark.R [--check]
#        --check runs the set twice and errors if the task count or the verdict
#        distribution moves between runs.
#
#  Writes .tmp/validation-benchmark.csv (git-ignored):
#    task | family | expected | observed | verdict | where | note
# =============================================================================

# ---- plumbing ----------------------------------------------------------------

benchmark_root <- function() {
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script)) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))))
  }
  wd <- getwd()
  repeat {
    if (file.exists(file.path(wd, "DESCRIPTION")) && dir.exists(file.path(wd, "R"))) {
      return(wd)
    }
    parent <- dirname(wd)
    if (identical(parent, wd)) {
      return(getwd())
    }
    wd <- parent
  }
}

bench_extdata <- function(name) {
  shipped <- system.file("extdata", name, package = "tidymarkers")
  if (nzchar(shipped)) {
    return(shipped)
  }
  file.path(benchmark_root(), "inst", "extdata", name)
}

bench_path <- function(...) file.path(benchmark_root(), ...)

bench_rows <- function(task, family, expected, observed, verdict, where, note) {
  n <- length(task)
  recycle <- function(x) {
    if (length(x) == n) x else rep(x, length.out = n)
  }
  tibble::tibble(
    task = as.character(task),
    family = recycle(as.character(family)),
    expected = recycle(as.character(expected)),
    observed = recycle(as.character(observed)),
    verdict = recycle(as.character(verdict)),
    where = recycle(as.character(where)),
    note = recycle(as.character(note))
  )
}

# A family that cannot run reports every one of its tasks as `not_run` with the
# error, instead of throwing out of the whole benchmark.
bench_not_run <- function(task, family, expected, reason) {
  bench_rows(task, family, expected, NA_character_, "not_run", NA_character_, reason)
}

bench_gap <- function(task, family, expected, reason) {
  bench_rows(task, family, expected, NA_character_, "gap", NA_character_, reason)
}

bench_has_tool <- function(tool) nzchar(Sys.which(tool))

bench_verdict <- function(pass) if (isTRUE(pass)) "pass" else "fail"

# ---- family: fixture ---------------------------------------------------------

# The hand-written panel in tests/testthat/fixtures/panel has expectations that
# were fixed in the test file, not read off the package's output: one target
# that must design, one that must be rejected with a cause from every chemistry.
bench_fixture <- function() {
  family <- "fixture"
  task <- c("fixture:Chr1:300 A>C", "fixture:Chr2:150 A>C")
  expected <- c("designed", "rejected with a cause from every chemistry")
  dir <- bench_path("tests", "testthat", "fixtures", "panel")
  genome_p <- file.path(dir, "genome.fa")
  hap_p <- file.path(dir, "haplotypes.csv")
  if (!file.exists(genome_p) || !file.exists(hap_p)) {
    return(bench_not_run(task, family, expected, "fixture panel absent from this tree"))
  }

  res <- tryCatch(
    {
      panel <- marker_design_panel(
        genome = Biostrings::readDNAStringSet(genome_p),
        haplotypes = hap_p,
        chemistries = c("snp", "hrm", "caps"),
        targets = tibble::tibble(
          chrom = c("Chr1", "Chr2"),
          pos = c(300L, 150L),
          ref = c("A", "A"),
          alt = c("C", "C")
        ),
        dir = tempfile()
      )
      list(ok = TRUE, panel = panel)
    },
    error = function(e) list(ok = FALSE, message = conditionMessage(e))
  )
  if (!res$ok) {
    return(bench_not_run(task, family, expected, res$message))
  }

  panel <- res$panel
  t <- panel$targets
  dec <- tryCatch(artefact_decisions(panel$solved), error = function(e) NULL)
  where_of <- function(tg) {
    if (is.null(dec) || !"target" %in% names(dec)) {
      return(NA_character_)
    }
    d <- dec[dec$target == tg, , drop = FALSE]
    if (!nrow(d)) {
      return(NA_character_)
    }
    paste(sprintf("%s:%s", d$scheme, d$verdict), collapse = "; ")
  }

  i <- match("Chr1:300 A>C", t$target)
  good <- if (is.na(i)) NULL else t[i, , drop = FALSE]
  good_ok <- !is.null(good) && isTRUE(good$ok[[1L]]) &&
    identical(as.character(good$scheme[[1L]]), "snp") &&
    inherits(good$assay[[1L]], "marker_assay") &&
    is.finite(good$footprint_dropout[[1L]]) &&
    good$footprint_dropout[[1L]] >= 0 &&
    good$footprint_dropout[[1L]] <= KASP_LIMITS$footprint_dropout
  observed_1 <- if (good_ok) {
    "designed"
  } else if (is.null(good)) {
    "no row for this target"
  } else {
    sprintf(
      "ok=%s scheme=%s assay=%s footprint=%s",
      good$ok[[1L]], good$scheme[[1L]],
      paste(class(good$assay[[1L]]), collapse = "/"),
      format(good$footprint_dropout[[1L]])
    )
  }

  j <- match("Chr2:150 A>C", t$target)
  bad <- if (is.na(j)) NULL else t[j, , drop = FALSE]
  causes <- if (is.null(bad)) character(0) else bad$cause[[1L]]
  bad_ok <- !is.null(bad) &&
    !isTRUE(bad$ok[[1L]]) &&
    is.na(bad$scheme[[1L]]) &&
    is.null(bad$assay[[1L]]) &&
    length(causes) > 0L &&
    setequal(sub(":.*$", "", causes), c("snp", "hrm", "caps")) &&
    all(grepl("^[a-z]+: .+", causes))
  observed_2 <- if (bad_ok) {
    "rejected with a cause"
  } else if (is.null(bad)) {
    "no row for this target"
  } else {
    sprintf(
      "ok=%s scheme=%s assay=%s causes=[%s]",
      bad$ok[[1L]], bad$scheme[[1L]],
      paste(class(bad$assay[[1L]]), collapse = "/"),
      paste(causes, collapse = " | ")
    )
  }

  bench_rows(
    task, family, expected,
    observed = c(observed_1, observed_2),
    verdict = c(bench_verdict(good_ok), bench_verdict(bad_ok)),
    where = c(where_of("Chr1:300 A>C"), where_of("Chr2:150 A>C")),
    note = c(
      "one flanking pair must survive on a clean template",
      "no chemistry may design on the homopolymer block"
    )
  )
}

# ---- family: soybean ---------------------------------------------------------

# 13 KASP assays the operator recorded as working in the field.  Re-running the
# design is not possible here (the genome build is not in tree), so the task
# reads the recorded gate table: an assay whose recorded gates reject it is a
# real contradiction between the package's rules and the field.
bench_soybean <- function() {
  family <- "soybean"
  gates_p <- bench_extdata("soybean_worked_panel_kasp_gates.csv")
  assays_p <- bench_extdata("soybean_worked_panel_assays.csv")
  if (!file.exists(gates_p) || !file.exists(assays_p)) {
    return(bench_gap(
      "soybean:(all)", family, "worked",
      "recorded soybean panel absent; cannot establish a task set"
    ))
  }
  gates <- utils::read.csv(gates_p, stringsAsFactors = FALSE)
  assays <- utils::read.csv(assays_p, stringsAsFactors = FALSE)
  ids <- sort(unique(gates$assay_id))
  task <- paste0("soybean:", ids)
  expected <- rep("worked", length(ids))

  rows <- lapply(seq_along(ids), function(k) {
    id <- ids[[k]]
    g <- gates[gates$assay_id == id, , drop = FALSE]
    if (!id %in% assays$assay_id) {
      return(bench_gap(task[[k]], family, "worked", "assay missing from the assays table"))
    }
    n_fail <- sum(g$status == "fail", na.rm = TRUE)
    n_unknown <- sum(g$status == "unknown", na.rm = TRUE)
    observed <- if (n_fail > 0L) {
      "gates_reject"
    } else if (n_unknown > 0L) {
      "gates_partial"
    } else {
      "gates_accept"
    }
    fails <- g[g$status == "fail", , drop = FALSE]
    where <- if (nrow(fails)) {
      tab <- sort(table(fails$criterion), decreasing = TRUE)
      paste(sprintf("%s x%d", names(tab), as.integer(tab)), collapse = "; ")
    } else {
      NA_character_
    }
    bench_rows(
      task[[k]], family, "worked", observed,
      bench_verdict(identical(observed, "gates_accept")),
      where,
      sprintf(
        "%d unknown gate rows; %d gate rows; outcome recorded as %s",
        n_unknown, nrow(g),
        unique(assays$outcome[assays$assay_id == id])[[1L]]
      )
    )
  })
  do.call(rbind, rows)
}

# ---- family: measured --------------------------------------------------------

# The four numbers VALIDATION.md quotes as measured.  Each is recomputed from a
# committed artifact, so the quoted claim is checked rather than repeated.
bench_measured <- function() {
  rows <- list(
    bench_measured_mm2020(),
    bench_measured_primer3_tm(),
    bench_measured_primer3_design(),
    bench_measured_specificity()
  )
  do.call(rbind, rows)
}

bench_measured_mm2020 <- function() {
  task <- "measured:mm2020 ranking"
  family <- "measured"
  expected <- "r > 0.9, spearman > 0.85, -25 < mean(tm - measured) < -10, sd < 2"
  path <- bench_extdata("tm_mm2020_s1_measured.csv")
  if (!file.exists(path)) {
    return(bench_not_run(task, family, expected, "tm_mm2020_s1_measured.csv absent"))
  }
  res <- tryCatch(
    {
      d <- utils::read.csv(path, stringsAsFactors = FALSE)
      can <- d[d$mismatch_class == "canonical", , drop = FALSE]
      tm <- primer_tm(
        can$top_18,
        conditions = reaction_conditions(Na = 50, Mg = 0, dNTPs = 0, DMSO = 0)
      )
      bias <- tm - can$measured_tm_c
      list(
        ok = TRUE,
        n = length(tm),
        r = stats::cor(tm, can$measured_tm_c),
        rho = stats::cor(tm, can$measured_tm_c, method = "spearman"),
        mean = mean(bias),
        sd = stats::sd(bias)
      )
    },
    error = function(e) list(ok = FALSE, message = conditionMessage(e))
  )
  if (!isTRUE(res$ok)) {
    return(bench_not_run(task, family, expected, res$message))
  }
  observed <- sprintf(
    "n %d, r %.3f, rho %.3f, mean %.2f, sd %.2f",
    res$n, res$r, res$rho, res$mean, res$sd
  )
  pass <- res$n == 64L &&
    all(is.finite(c(res$r, res$rho, res$mean, res$sd))) &&
    res$r > 0.9 && res$rho > 0.85 &&
    res$mean > -25 && res$mean < -10 && res$sd < 2
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "test-tm-mm2020.R:88-113",
    "64 canonical duplexes; the level bias is the salt term, not corrected"
  )
}

bench_measured_primer3_tm <- function() {
  task <- "measured:primer3 Tm offset"
  family <- "measured"
  expected <- "n 34, mean +3.92, sd 0.25, max 4.66 (primer3 - pkg)"
  path <- bench_path("results", "primer3-thermo.csv")
  if (!file.exists(path)) {
    return(bench_not_run(task, family, expected, "results/primer3-thermo.csv absent"))
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  t <- d[d$measure == "tm_c", , drop = FALSE]
  if (!nrow(t)) {
    return(bench_not_run(task, family, expected, "no tm_c rows in primer3-thermo.csv"))
  }
  delta <- as.numeric(t$primer3) - as.numeric(t$pkg)
  observed <- sprintf(
    "n %d, mean %+.2f, sd %.2f, max %.2f", length(delta),
    mean(delta), stats::sd(delta), max(delta)
  )
  pass <- length(delta) == 34L &&
    abs(mean(delta) - 3.92) < 0.05 &&
    stats::sd(delta) < 0.35 &&
    abs(max(delta) - 4.66) < 0.1
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "results/primer3-thermo.csv",
    "cross-tool offset recorded from a local primer3 run; additive, not modelled"
  )
}

bench_measured_primer3_design <- function() {
  task <- "measured:primer3 design (one locus)"
  family <- "measured"
  expected <- "8/20 variants have a primer3 pair passing every gate"
  path <- bench_path("results", "primer3-design.csv")
  if (!file.exists(path)) {
    return(bench_not_run(task, family, expected, "results/primer3-design.csv absent"))
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  flag <- d$p3_any_pair_pass
  n <- sum(flag %in% TRUE)
  observed <- sprintf("%d/%d", n, nrow(d))
  pass <- nrow(d) == 20L && n == 8L
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "results/primer3-design.csv",
    paste(
      "only one of the three loci in the VALIDATION.md claim is recorded here;",
      "the other two need a live primer3_core and are not invented"
    )
  )
}

bench_measured_specificity <- function() {
  task <- "measured:specificity oracle"
  family <- "measured"
  expected <- "3'-seed 8/8, exact PCR 3/3, BLAST recall 1.000"
  path <- bench_path("results", "specificity-oracle.csv")
  if (!file.exists(path)) {
    return(bench_not_run(task, family, expected, "results/specificity-oracle.csv absent"))
  }
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  seed <- d[d$test == "seed", , drop = FALSE]
  pcr <- d[d$test == "pcr", , drop = FALSE]
  blast <- d[d$test == "blast_full_length", , drop = FALSE]
  seed_ok <- nrow(seed) == 8L && all(seed$equal %in% TRUE)
  pcr_ok <- nrow(pcr) == 3L && all(pcr$equal %in% TRUE)
  truth <- sum(as.numeric(blast$truth))
  found <- sum(as.numeric(blast$found))
  recall <- if (truth > 0) found / truth else NA_real_
  blast_ok <- nrow(blast) == 8L && isTRUE(abs(recall - 1) < 1e-9)
  observed <- sprintf(
    "seed %d/%d, pcr %d/%d, blast recall %.3f",
    sum(seed$equal %in% TRUE), nrow(seed),
    sum(pcr$equal %in% TRUE), nrow(pcr), recall
  )
  bench_rows(task, family, expected, observed,
    bench_verdict(seed_ok && pcr_ok && blast_ok),
    "results/specificity-oracle.csv",
    "BLAST rows carry equal = NA by design; recall is compared instead"
  )
}

# ---- family: contract --------------------------------------------------------

# The six prose claims the ZX3 audit flagged: the vignette says the code does X,
# so the benchmark asks the live code directly.  A matching answer is a pass;
# a divergence is a red task, not a silent skip.
bench_contract <- function() {
  rows <- list(
    bench_contract_raw_hits(),
    bench_contract_off_by_one(),
    bench_contract_crispr_flow(),
    bench_contract_duplex_na(),
    bench_contract_decision_graph(),
    bench_contract_ladder_na()
  )
  do.call(rbind, rows)
}

# A 20 nt spacer that is not its own reverse complement, and a genome holding it
# once.  `strrep` padding is T/G, which the spacer cannot match.
bench_contract_spacer <- function() "ACGTTGCAAGTCCTAGGATC"

bench_contract_raw_hits <- function() {
  task <- "contract:crispr_spacer_hits counts raw occurrences"
  family <- "contract"
  expected <- "1 occurrence -> 1; 2 occurrences -> 2"
  res <- tryCatch(
    {
      spacer <- bench_contract_spacer()
      one <- Biostrings::DNAStringSet(c(
        Chr1 = paste0(strrep("T", 10L), spacer, strrep("G", 10L))
      ))
      two <- Biostrings::DNAStringSet(c(
        Chr1 = paste0(strrep("T", 10L), spacer, strrep("G", 10L), spacer, strrep("T", 10L))
      ))
      c(
        once = crispr_spacer_hits(spacer, one),
        twice = crispr_spacer_hits(spacer, two)
      )
    },
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(bench_not_run(task, family, expected, conditionMessage(res)))
  }
  observed <- sprintf("once %d, twice %d", res[["once"]], res[["twice"]])
  pass <- identical(unname(res), c(1L, 2L))
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "R/marker_crispr.R:268",
    "counts forward and reverse-complement occurrences; the design subtracts 1"
  )
}

bench_contract_off_by_one <- function() {
  task <- "contract:off_target == hits - 1"
  family <- "contract"
  expected <- "design_crispr_marker()$off_target equals crispr_spacer_hits(guide) - 1"
  if (!bench_has_tool("hybrid-ss-min")) {
    return(bench_not_run(task, family, expected, "OligoArrayAux absent (hybrid-ss-min not on PATH)"))
  }
  res <- tryCatch(
    {
      w <- paste(sample(c("A", "C", "G", "T"), 600L, TRUE), collapse = "")
      substr(w, 290L, 303L) <- paste0("TTTA", "GCTAGCATGCATGC")
      substr(w, 300L, 300L) <- "T"
      genome <- Biostrings::DNAStringSet(c(Chr1 = w))
      snp <- tibble::tibble(chrom = "Chr1", pos = 300L, ref = "T", alt = "A")
      m <- suppressWarnings(design_crispr_marker(genome, snp, flank = 200L, cas = "cas13"))
      if (is.null(m)) NULL else list(off = m$off_target, hits = crispr_spacer_hits(m$guide, genome))
    },
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(bench_not_run(task, family, expected, conditionMessage(res)))
  }
  if (is.null(res)) {
    return(bench_not_run(task, family, expected, "no guide passed the Cas13 filters on this template"))
  }
  observed <- sprintf("off_target %s, hits %s", format(res$off), format(res$hits))
  pass <- identical(as.integer(res$off), as.integer(res$hits) - 1L)
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "R/marker_crispr.R:471",
    "the on-target occurrence is not an off-target"
  )
}

bench_contract_crispr_flow <- function() {
  task <- "contract:crispr flow uses a preamp pair, no RPA gate"
  family <- "contract"
  expected <- "no gate_rpa_structure in scheme_def_crispr(); preamp defaults to \"rpa\""
  res <- tryCatch(
    {
      list(
        gates = scheme_def_crispr()$gates,
        preamp = formals(design_crispr_marker)$preamp
      )
    },
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(bench_not_run(task, family, expected, conditionMessage(res)))
  }
  pass <- !"gate_rpa_structure" %in% res$gates && identical(res$preamp, "rpa")
  observed <- sprintf(
    "rpa gate present: %s; preamp default: %s",
    "gate_rpa_structure" %in% res$gates, deparse(res$preamp)
  )
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "R/marker_crispr.R:600",
    "the vignette prose already matches this"
  )
}

bench_contract_duplex_na <- function() {
  task <- "contract:design_duplex warns on undeclared homology"
  family <- "contract"
  expected <- "warning 'undeclared', gate_duplex_targets not_applicable, verdict fail"
  if (!bench_has_tool("hybrid-ss-min")) {
    return(bench_not_run(task, family, expected, "OligoArrayAux absent (hybrid-ss-min not on PATH)"))
  }
  genome_p <- bench_path("tests", "testthat", "fixtures", "panel", "genome.fa")
  if (!file.exists(genome_p)) {
    return(bench_not_run(task, family, expected, "fixture panel absent from this tree"))
  }
  res <- tryCatch(
    {
      genome <- Biostrings::readDNAStringSet(genome_p)
      loci <- tibble::tibble(
        chrom = c("Chr1", "Chr1"),
        pos = c(100L, 300L),
        ref = c("C", "A"),
        alt = c("G", "C")
      )
      warned <- FALSE
      d <- withCallingHandlers(
        design_duplex(genome, loci, readout = "co-dominant", homologous = NA),
        warning = function(w) {
          if (grepl("undeclared", conditionMessage(w))) warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
      compat <- d$compatibility
      status <- if (is.null(compat) || !nrow(compat)) {
        NA_character_
      } else {
        as.character(compat$status[compat$gate == "gate_duplex_targets"])
      }
      list(warned = warned, status = status, verdict = d$verdict)
    },
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(bench_not_run(task, family, expected, conditionMessage(res)))
  }
  pass <- isTRUE(res$warned) &&
    identical(res$status, "not_applicable") &&
    identical(res$verdict, "fail")
  observed <- sprintf(
    "warned %s, gate_duplex_targets %s, verdict %s",
    res$warned, res$status, res$verdict
  )
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "R/marker_duplex.R:400-405",
    "an undeclared topology cannot pass, so the pairing fails"
  )
}

bench_contract_decision_graph <- function() {
  task <- "contract:decision graph stamps gate status from the selection"
  family <- "contract"
  expected <- "one gate: no record ok, relaxed record relaxed, violated record violated"
  res <- tryCatch(
    {
      # Gate nodes are emitted by the selection segment, not by the criterion
      # segment: a scheme-only graph holds the criterion names.
      nodes <- tibble::as_tibble(
        marker_decision_graph(list(scheme = "snp", selection = list()))
      )$name
      gates <- criterion_spec(KASP_LIMITS)
      gates <- sub("^.*:gate:", "", gates$id[gates$role == "gate"])
      candidates <- nodes[grepl(":gate:", nodes)]
      bare <- sub("^.*:gate:", "", candidates)
      hit <- match(gates, bare)
      hit <- hit[!is.na(hit)]
      if (length(hit)) {
        node <- candidates[[hit[[1L]]]]
        bare <- sub("^.*:gate:", "", node)
        status <- function(sel) {
          g <- tibble::as_tibble(marker_decision_graph(list(scheme = "snp", selection = sel)))
          as.character(g$status[g$name == node])[1L]
        }
        got <- c(
          status(list()),
          status(list(relaxed = bare)),
          status(list(violated = bare))
        )
        list(ok = TRUE, pass = identical(got, c("ok", "relaxed", "violated")),
          observed = sprintf("%s -> %s", node, paste(got, collapse = "/")))
      } else {
        list(ok = FALSE, pass = FALSE, observed = "no known gate node on the snp graph")
      }
    },
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(bench_not_run(task, family, expected, conditionMessage(res)))
  }
  bench_rows(task, family, expected, res$observed, bench_verdict(res$pass),
    "tests/testthat/test-falsifiable.R:405-440",
    "the graph reads the selection records; it does not restate the rules"
  )
}

bench_contract_ladder_na <- function() {
  task <- "contract:select_with_ladder refuses an NA margin"
  family <- "contract"
  expected <- "NA COM Tm -> best NULL; 60 -> best non-NULL"
  res <- tryCatch(
    {
      pair_table <- function(com_tm) {
        tibble::tibble(
          discrimination = 0.9,
          com_homodimer = -1, com_hairpin = -0.5, com_tm = com_tm,
          amplicon = 80, n_var = 0L, as_n_3p = 0L, run = 0L, ssr = 0L,
          com_start = seq_along(com_tm)
        )
      }
      pick <- function(x) {
        d <- pair_table(x)
        g <- grade_assay_pairs(d, KASP_LIMITS, as_tm_mean = 58)
        select_with_ladder(
          g$data, g$gates, KASP_RELAX_PAIR, g$tiers, g$ties, g$gate_margin, g$gate_floor
        )$best
      }
      list(na_best = pick(NA_real_), ok_best = pick(60))
    },
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(bench_not_run(task, family, expected, conditionMessage(res)))
  }
  pass <- is.null(res$na_best) && !is.null(res$ok_best)
  observed <- sprintf(
    "NA -> %s; 60 -> %s",
    if (is.null(res$na_best)) "NULL" else "non-NULL",
    if (is.null(res$ok_best)) "NULL" else "non-NULL"
  )
  bench_rows(task, family, expected, observed, bench_verdict(pass),
    "R/marker_select.R:458-465",
    "an unmeasured margin is never a pass"
  )
}

# ---- driver ------------------------------------------------------------------

benchmark_tasks <- function() {
  rows <- list(
    bench_fixture(),
    bench_soybean(),
    bench_measured(),
    bench_contract()
  )
  out <- do.call(rbind, rows)
  out[, c("task", "family", "expected", "observed", "verdict", "where", "note")]
}

benchmark_write <- function(x, path = NULL) {
  if (is.null(path)) {
    path <- bench_path(".tmp", "validation-benchmark.csv")
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  path
}

benchmark_distribution <- function(x) {
  tab <- table(factor(x$verdict, levels = c("pass", "fail", "not_run", "gap")))
  as.integer(tab)
}

main <- function() {
  check <- "--check" %in% commandArgs(trailingOnly = TRUE)
  if (!"tidymarkers" %in% loadedNamespaces()) {
    suppressMessages(pkgload::load_all(benchmark_root(), quiet = TRUE))
  }
  x <- benchmark_tasks()
  path <- benchmark_write(x)
  n <- nrow(x)
  dist <- benchmark_distribution(x)
  names(dist) <- c("pass", "fail", "not_run", "gap")
  message(sprintf(
    "validation benchmark: %d task(s); %s", n,
    paste(sprintf("%s %d", names(dist), dist), collapse = ", ")
  ))
  message("wrote ", path)
  if (check) {
    y <- benchmark_tasks()
    if (nrow(y) != n || !identical(benchmark_distribution(y), unname(dist))) {
      stop("--check: the task count or verdict distribution moved between runs", call. = FALSE)
    }
    message("--check: reproduced")
  }
  invisible(x)
}

if (sys.nframe() == 0L) {
  main()
}
