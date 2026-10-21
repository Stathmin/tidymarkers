# Minimal OLE2 + Word-97 piece-table text extractor.
#
# Used only by data-raw/dong4platform_poisson.R to read the .doc supplementary
# file of 10.1038/srep13174.  There is no antiword/catdoc/LibreOffice on this
# machine and pandoc cannot read Word 97, so the pieces are decoded by hand.
# The implementation was validated byte-for-byte against an independent Python
# reference decoder (same 18343 characters) before landing.
#
# Output: the full document text, pieces concatenated, with \r line breaks.

wd_text <- function(path) {
  raw <- readBin(path, "raw", file.size(path))
  if (!identical(raw[1:8], as.raw(c(0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1)))) {
    stop("not an OLE2 compound file: ", path)
  }

  u16 <- function(b, i) as.integer(b[i]) + 256L * as.integer(b[i + 1L])
  u32 <- function(b, i) {
    # unsigned 32-bit as a double; 0xFFFFFFFF stays finite (no integer overflow)
    as.numeric(u16(b, i)) + 65536 * as.numeric(u16(b, i + 2L))
  }

  # CFB header fields are 0-based byte offsets; u16()/u32() take 1-based indices
  secsz <- 2^u16(raw, 31L)
  minisz <- 2^u16(raw, 33L)
  nfat <- u32(raw, 45L)
  dirstart <- u32(raw, 49L)
  cutoff <- u32(raw, 57L)
  minifat <- u32(raw, 61L)
  sector <- function(i) raw[(513 + i * secsz):(512 + (i + 1) * secsz)]

  fat <- numeric(0)
  for (i in 0:(nfat - 1)) {
    s <- sector(u32(raw, 77 + 4 * i))
    n <- length(s) %/% 4
    fat <- c(fat, vapply(0:(n - 1), function(k) u32(s, 4 * k + 1), 0))
  }
  follow <- function(start, tab, limit) {
    out <- numeric(0)
    while (length(start) == 1 && is.finite(start) && start >= 0 &&
      start < limit && length(out) < 100000) {
      out <- c(out, start)
      # FAT is indexed by sector number: fat[sector]
      start <- if (start + 1 <= length(tab)) tab[start + 1] else NA_real_
    }
    out
  }
  limit <- length(raw) %/% secsz + 8

  dir_bytes <- unlist(lapply(follow(dirstart, fat, limit), sector), use.names = FALSE)
  n_ent <- length(dir_bytes) %/% 128L
  ents <- lapply(seq_len(n_ent), function(j) {
    k <- (j - 1L) * 128L + 1L
    e <- dir_bytes[k:(k + 127)]
    # directory entry: name length u16 @64, type @66, start @116, size @120
    nl <- if (length(e) >= 66) u16(e, 65L) else 0L
    if (nl < 2 || nl > 128) {
      return(NULL)
    }
    nb <- e[seq(1, nl - 2, 2)]
    hi <- e[seq(2, nl - 2, 2)]
    list(
      name = intToUtf8(as.integer(nb) + 256L * as.integer(hi)),
      type = as.integer(e[67]),
      start = u32(e, 117L), size = u32(e, 121L)
    )
  })
  ents <- Filter(Negate(is.null), ents)

  root <- Filter(function(e) e$type == 5, ents)[[1]]
  ministream <- unlist(lapply(follow(root$start, fat, limit), sector), use.names = FALSE)
  mfat <- numeric(0)
  for (i in follow(minifat, fat, limit)) {
    s <- sector(i)
    n <- length(s) %/% 4
    mfat <- c(mfat, vapply(0:(n - 1), function(k) u32(s, 4 * k + 1), 0))
  }
  stream <- function(name) {
    hit <- Filter(function(z) z$name == name && z$type == 2, ents)
    if (!length(hit)) stop("stream missing: ", name, "; have: ", paste(vapply(ents, function(z) z$name, ""), collapse = "|"))
    e <- hit[[1]]
    if (e$size < cutoff) {
      buf <- unlist(lapply(follow(e$start, mfat, limit), function(i) {
        ministream[(i * minisz + 1):((i + 1) * minisz)]
      }), use.names = FALSE)
    } else {
      buf <- unlist(lapply(follow(e$start, fat, limit), sector), use.names = FALSE)
    }
    buf[seq_len(e$size)]
  }

  wd <- stream("WordDocument")
  tbl <- stream("1Table")

  i32 <- function(b, i) {
    v <- sum(as.numeric(b[i:(i + 3)]) * c(1, 256, 65536, 16777216))
    if (v >= 2^31) v - 2^32 else v
  }

  # FIB (Word 97): csw(u16), rgw97, cslw(u16), rglw97, cbRgFcLcb(u16), rgFcLcb.
  # Index taken from the FIB layout; asserted to be 14 for this file.
  csw <- u16(wd, 33)
  if (csw != 14L) stop("unexpected FIB csw: ", csw)
  p <- 35L + csw * 2L
  p <- p + 2L + u16(wd, p) * 4L
  cb_rg <- u16(wd, p)
  p <- p + 2L
  rg <- vapply(seq_len(cb_rg * 2L), function(k) i32(wd, p + 4L * (k - 1L)), 0)
  fcClx <- rg[67]
  lcbClx <- rg[68]

  clx <- tbl[(fcClx + 1L):(fcClx + lcbClx)]
  i <- 1L
  repeat {
    if (as.integer(clx[i]) == 1L) {
      # Prc: clxt(1), cbGrpprl(u16), GrpPrl
      i <- i + 3L + u16(clx, i + 2L)
    } else if (as.integer(clx[i]) == 2L) {
      # Pcdt: clxt(1), lcb(u32), PlcPcd
      lcb <- i32(clx, i + 1L)
      plc <- clx[(i + 5L):(i + 4L + lcb)]
      break
    } else {
      stop("unexpected CLX record byte ", as.integer(clx[i]))
    }
  }

  n_pcd <- (length(plc) - 4L) %/% 12L
  cps <- vapply(seq_len(n_pcd + 1L), function(k) i32(plc, 4L * (k - 1L) + 1L), 0)
  base <- 4L * (n_pcd + 1L)
  paste0(vapply(seq_len(n_pcd), function(k) {
    pcd <- plc[(base + 8L * (k - 1L) + 1L):(base + 8L * k)]
    fc_comp <- i32(pcd, 3L) # fcCompressed at PCD offset 2
    comp <- bitwAnd(bitwShiftR(as.integer(min(fc_comp, 2147483647)), 30), 1L)
    fc <- fc_comp %% 2^30
    len <- cps[k + 1L] - cps[k]
    if (comp == 1L) {
      fc <- fc %/% 2
      iconv(rawToChar(wd[(fc + 1L):(fc + len)]), "CP1252", "UTF-8")
    } else {
      intToUtf8(as.integer(wd[(fc + 1L):(fc + 2L * len)][seq(1, 2L * len, 2L)]) +
        256L * as.integer(wd[(fc + 1L):(fc + 2L * len)][seq(2, 2L * len, 2L)]))
    }
  }, ""), collapse = "")
}
