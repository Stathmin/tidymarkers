#!/usr/bin/env Rscript
# ZW3.10 measurement: the Tm re-anchor from Schildkraut2010 to SantaLucia1998-2.
#
# Measurement only.  Nothing here changes TM_SALT_METHOD, a design window,
# KASP_LIMITS, tm_tolerances() or tm_range.  It runs the real designers on the
# test fixture genome and asks what the SL2 coordinate would do to the primers
# they emit.
#
# How the two corrections actually enter TmCalculator (read from
# `.tm_nn_chunk`, not assumed):
#   Schildkraut2010    tm <- 1000*dh/(ds + R*log(k)) - 273.15 + 16.6*log10(mon)
#   SantaLucia1998-2   ds <- ds + 0.368*(N - 1)*ln(mon), then the same 1/Tm form
# so SL2 is an *entropy* shift, not an additive Tm offset.  The identity check
# below therefore reconstructs dh/ds exactly from two `salt_method = "none"`
# probes (vary dnac to separate the two unknowns) and predicts both Tms.
#
# `mon` in TmCalculator is mol/L (`Mon/1000`), so ln(mon) is negative here.

pkgload::load_all(".", quiet = TRUE)

SCHILD <- "Schildkraut2010"
SL2 <- "SantaLucia1998-2"
R_GAS <- 1.987

# Every scheme shares SCHEME_BUFFER (empty): same Na/Mg/dNTPs/DMSO, and the
# scheme temperature does not enter the NN Tm.  State the conditions once.
cond <- scheme_conditions("pcr")
for (sc in c("snp", "dpcr", "lamp", "snapshot")) {
  other <- scheme_conditions(sc)
  stopifnot(
    isTRUE(all.equal(
      c(other$Na, other$Mg, other$dNTPs, other$DMSO),
      c(cond$Na, cond$Mg, cond$dNTPs, cond$DMSO)
    ))
  )
}

mon_mol <- (cond$Na + 120 * sqrt(cond$Mg - cond$dNTPs)) / 1000
ln_mon <- log(mon_mol)
corr_s <- 16.6 * log10(mon_mol) # Schildkraut2010, additive in Tm
corr_l <- function(n) 0.368 * (n - 1) * ln_mon # SantaLucia1998-2, additive in dS
# primer_tm() passes dnac_high = dnac_low = 25 nM, so k = (25 - 25/2) * 1e-9.
k_design <- 25 / 2 * 1e-9

get_tm <- function(seq, scheme, method, dnac = NULL) {
  primer_tm(
    seq,
    conditions = scheme_conditions(scheme),
    salt_method = method,
    dnac = dnac
  )
}

# Recover the NN enthalpy (A = 1000*dh) and entropy (S = ds) terms from two
# salt-free probes at different dnac.  A and S are salt-method independent.
nn_terms <- function(seq, scheme) {
  d1 <- 25
  d2 <- 2500
  k1 <- d1 / 2 * 1e-9
  k2 <- d2 / 2 * 1e-9
  t1 <- get_tm(seq, scheme, "none", dnac = d1) + 273.15
  t2 <- get_tm(seq, scheme, "none", dnac = d2) + 273.15
  A <- R_GAS * (log(k1) - log(k2)) / (1 / t1 - 1 / t2)
  S <- A / t1 - R_GAS * log(k1)
  list(A = A, S = S)
}

# Exact prediction of both salt coordinates from the reconstructed NN terms.
predict_pair <- function(seq, scheme) {
  nn <- nn_terms(seq, scheme)
  base <- nn$A / (nn$S + R_GAS * log(k_design)) - 273.15
  c(
    schild = base + corr_s,
    sl2 = nn$A / (nn$S + corr_l(nchar(seq)) + R_GAS * log(k_design)) - 273.15
  )
}

rows <- list()
lens <- list()
checks <- numeric(0)

add_len <- function(seq, delta) {
  lens[[length(lens) + 1L]] <<- tibble::tibble(
    len = nchar(seq), delta = delta
  )
}

push <- function(scheme, role, primer, lo, hi, tm_s, tm_l) {
  rows[[length(rows) + 1L]] <<- tibble::tibble(
    scheme, role, primer,
    len = nchar(primer),
    window_lo = lo, window_hi = hi,
    tm_schild = tm_s, tm_sl2 = tm_l, delta = tm_l - tm_s
  )
}

add1 <- function(scheme, role, primer, window) {
  primer <- as.character(primer)[1L]
  if (is.na(primer) || !nzchar(primer)) {
    return(invisible())
  }
  ts <- get_tm(primer, scheme, SCHILD)
  tl <- get_tm(primer, scheme, SL2)
  pred <- predict_pair(primer, scheme)
  checks <<- c(checks, abs(pred[["schild"]] - ts), abs(pred[["sl2"]] - tl))
  add_len(primer, tl - ts)
  push(scheme, role, primer, window[1L], window[2L], ts, tl)
}

# LAMP's inner gate is pmax(arm5 Tm, arm3 Tm), so the design-valid unit is the
# fused FIP/BIP and its Tm is the max over the two arms.  Measure both arms,
# check the identity on each, and record the gated pair for the window analysis.
add_lamp_inner <- function(marker, genome, window) {
  arm <- function(co) {
    as.character(Biostrings::subseq(genome[[1L]], co[1L], co[2L]))
  }
  for (nm in c("FIP", "BIP")) {
    co <- marker$lamp_arms[[nm]]
    if (is.null(co)) next
    seqs <- c(arm(co$arm5), arm(co$arm3))
    ts <- vapply(seqs, get_tm, numeric(1), scheme = "lamp", method = SCHILD)
    tl <- vapply(seqs, get_tm, numeric(1), scheme = "lamp", method = SL2)
    for (i in seq_along(seqs)) {
      pred <- predict_pair(seqs[i], "lamp")
      checks <<- c(
        checks,
        abs(pred[["schild"]] - ts[i]),
        abs(pred[["sl2"]] - tl[i])
      )
      add_len(seqs[i], tl[i] - ts[i])
    }
    push(
      "lamp", sprintf("LAMP inner (%s gated)", nm),
      seqs[which.max(ts)], window[1L], window[2L], max(ts), max(tl)
    )
  }
}

# windows, read from the same defaults the designers apply
kasp_as <- KASP_LIMITS$as_tm
kasp_com <- KASP_LIMITS$com_tm
dpcr_tm <- eval(formals(design_dpcr_marker)$tm_range)
lamp_inner <- eval(formals(design_lamp_marker)$inner_tm)
lamp_outer <- eval(formals(design_lamp_marker)$outer_tm)
snap_ext <- SNAPSHOT_EXT_TM
snap_amp <- KASP_LIMITS$com_tm

n_seeds <- 10L
positions <- c(260L, 300L, 340L)
built <- c(snp = 0L, dpcr = 0L, lamp = 0L, snapshot = 0L)
lamp_tm_resid <- numeric(0)

for (seed in seq_len(n_seeds)) {
  for (pos in positions) {
    fx <- snp_fixture(n = 600L, seed = seed, pos = pos)

    k <- suppressMessages(tryCatch(
      assemble_kasp_assay(fx$genome, fx$snp, flank = 150L),
      error = function(e) NULL
    ))
    if (!is.null(k) && !inherits(k, "marker_doomed") && !is.null(k$assay)) {
      built[["snp"]] <- built[["snp"]] + 1L
      add1("snp", "KASP AS", k$assay$as_ref[1L], kasp_as)
      add1("snp", "KASP AS", k$assay$as_alt[1L], kasp_as)
      add1("snp", "KASP COM", k$assay$com[1L], kasp_com)
    }

    d <- suppressMessages(tryCatch(
      design_dpcr_marker(fx$genome, fx$snp, flank = 150L),
      error = function(e) NULL
    ))
    if (!is.null(d) && !is.null(d$primers)) {
      built[["dpcr"]] <- built[["dpcr"]] + 1L
      add1("dpcr", "dPCR primer", d$primers$forward[1L], dpcr_tm)
      add1("dpcr", "dPCR primer", d$primers$reverse[1L], dpcr_tm)
    }

    l <- suppressMessages(tryCatch(
      design_lamp_marker(fx$genome, fx$snp, flank = 150L),
      error = function(e) NULL
    ))
    if (!is.null(l) && !is.null(l$lamp_arms)) {
      built[["lamp"]] <- built[["lamp"]] + 1L
      add1("lamp", "LAMP outer (F3/B3)", l$primers$F3[1L], lamp_outer)
      add1("lamp", "LAMP outer (F3/B3)", l$primers$B3[1L], lamp_outer)
      before <- length(rows)
      add_lamp_inner(l, fx$genome, lamp_inner)
      # sanity: the gated Schild Tm should reproduce the designer's own value
      if (length(rows) > before) {
        new <- do.call(rbind, rows[(before + 1L):length(rows)])
        for (i in seq_len(nrow(new))) {
          nm <- sub("^LAMP inner \\((FIP|BIP) gated\\)$", "\\1", new$role[i])
          if (nm %in% names(l$tm)) {
            lamp_tm_resid <- c(lamp_tm_resid, new$tm_schild[i] - l$tm[[nm]])
          }
        }
      }
    }

    s <- suppressMessages(tryCatch(
      design_snapshot_marker(fx$genome, fx$snp, flank = 250L),
      error = function(e) NULL
    ))
    if (!is.null(s) && !inherits(s, "marker_doomed") && !is.null(s$extension)) {
      built[["snapshot"]] <- built[["snapshot"]] + 1L
      add1("snapshot", "SNaPShot amp", s$amplification$forward[1L], snap_amp)
      add1("snapshot", "SNaPShot amp", s$amplification$reverse[1L], snap_amp)
      add1("snapshot", "SNaPShot ext", s$extension$primer[1L], snap_ext)
    }
  }
}

res <- do.call(rbind, rows)
len_tab <- do.call(rbind, lens)
stopifnot(nrow(res) > 0L, length(checks) > 0L)
max_check <- max(checks)
if (max_check >= 1e-6) {
  cat(sprintf("identity check FAILED: max |measured - predicted| = %.6g C\n", max_check))
  diag <- res |>
    dplyr::mutate(pred = predict_pair(primer, scheme)[["sl2"]]) |>
    dplyr::arrange(dplyr::desc(abs(tm_sl2 - pred))) |>
    utils::head(10L)
  print(as.data.frame(diag), digits = 6, row.names = FALSE)
  stop("TmCalculator's salt terms are not the assumed pair; inspect the rows above.")
}

dir.create(".tmp-sdd/roadmap", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, ".tmp-sdd/roadmap/task-4-primers.csv", row.names = FALSE)
utils::write.csv(len_tab, ".tmp-sdd/roadmap/task-4-length-deltas.csv", row.names = FALSE)

cat("\n## Conditions\n\n")
cat(sprintf(
  "Na %.1f mM, Mg %.1f mM, dNTPs %.1f mM, DMSO %.0f %%, dnac 25 nM each strand (50 nM total)\n",
  cond$Na, cond$Mg, cond$dNTPs, cond$DMSO
))
cat(sprintf(
  "monovalent equivalent mon = Na + 120*sqrt(Mg - dNTPs) = %.3f mM = %.6f mol/L\n",
  mon_mol * 1000, mon_mol
))
cat(sprintf(
  "Schildkraut2010 corr = 16.6*log10(mon) = %.4f C; SantaLucia1998-2 corr(N) = %.4f*(N-1) dS\n",
  corr_s, 0.368 * ln_mon
))
cat(sprintf(
  "identity check: max |measured - predicted| = %.3g C on %d Tm values\n\n",
  max_check, length(checks)
))
if (length(lamp_tm_resid) > 0L) {
  cat(sprintf(
    "LAMP sanity: gated Schild Tm vs designer l$tm: max |diff| = %.3g C on %d arms\n\n",
    max(abs(lamp_tm_resid)), length(lamp_tm_resid)
  ))
}

cat("## Fixtures\n\n")
cat(paste(
  sprintf("%s: %d/%d", names(built), built, n_seeds * length(positions)),
  collapse = "; "
), "\n\n")

cat("## Offset by population (C)\n\n")
agg <- res |>
  dplyr::group_by(scheme, role) |>
  dplyr::summarise(
    n = dplyr::n(),
    len_min = min(len), len_med = stats::median(len), len_max = max(len),
    delta_min = min(delta), delta_mean = mean(delta),
    delta_med = stats::median(delta), delta_max = max(delta),
    .groups = "drop"
  )
print(as.data.frame(agg), digits = 3, row.names = FALSE)

cat("\n## Blast radius per window\n\n")
win <- res |>
  dplyr::mutate(
    in_now = tm_schild >= window_lo & tm_schild <= window_hi,
    out_sl2 = tm_sl2 < window_lo | tm_sl2 > window_hi
  ) |>
  dplyr::group_by(scheme, role, window_lo, window_hi) |>
  dplyr::summarise(
    n = dplyr::n(),
    valid_now = sum(in_now),
    already_out = sum(!in_now),
    fall_out = sum(in_now & out_sl2),
    frac_fall_out = fall_out / max(valid_now, 1L),
    fall_low = sum(in_now & tm_sl2 < window_lo),
    fall_high = sum(in_now & tm_sl2 > window_hi),
    .groups = "drop"
  )
print(as.data.frame(win), digits = 3, row.names = FALSE)

cat("\n## Delta against length (measured, all schemes)\n\n")
lt <- len_tab |>
  dplyr::group_by(len) |>
  dplyr::summarise(
    n = dplyr::n(), delta_mean = mean(delta),
    delta_min = min(delta), delta_max = max(delta),
    delta_sd = stats::sd(delta), .groups = "drop"
  )
print(as.data.frame(lt), digits = 3, row.names = FALSE)

cat("\nprimer table: .tmp-sdd/roadmap/task-4-primers.csv\n")
cat("length table: .tmp-sdd/roadmap/task-4-length-deltas.csv\n")
