# =============================================================================
#  data-raw/logo.R  --  generate the package hex logo (tidyverse paradigm)
# =============================================================================
#  A tidyverse hex sticker: a hexagon with a border colour and an image inside.
#  The image is the package idea in one picture: a DNA double helix whose two
#  strands swirl *around* the Tree of Life (data-raw/tree-of-life.svg), the ten
#  Sephiroth and the paths between them named in Hebrew.
#
#  The helix is drawn as a 3D object projected onto the sticker: a strand is
#  either in front of the diagram or behind it, and the drawing is split into a
#  back layer (painted before the tree, so its knockout hides the strands) and a
#  front layer (painted after the tree, so they cross over it).  Each backbone is
#  cut into the stretches that pass in front and the stretches that pass behind,
#  so the pair reads as a double helix winding around the diagram.
#
#  Written as plain SVG (no graphics device, no cairo needed) and optionally
#  rasterised to PNG with ImageMagick or librsvg.
#
#  Run:  Rscript data-raw/logo.R
#  Out:  man/figures/logo.svg  (+ logo.png if a rasteriser is available)
# =============================================================================

script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(sub("^--file=", "", script[[1L]])))) else getwd()

out_svg <- file.path(root, "man", "figures", "logo.svg")
out_png <- file.path(root, "man", "figures", "logo.png")
tree_svg <- file.path(root, "data-raw", "tree-of-life.svg")
dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)

## geometry ------------------------------------------------------------------
W <- 240; H <- 277                 # canvas (tidyverse hex aspect)
cx <- W / 2; cy <- H / 2; r <- 116 # hexagon centre + circumradius
ang <- (90 + 0:5 * 60) * pi / 180  # pointy-top hexagon
hexx <- cx + r * cos(ang)
hexy <- cy + r * sin(ang)
hex_pts <- paste(sprintf("%.1f,%.1f", hexx, hexy), collapse = " ")

## colours (base palette ~ KASP FAM/HEX + a tidy accent) ----------------------
ink     <- "#1b4965"               # border / text
accent  <- "#d1495b"               # "marker" accent
paper   <- "#fffdf7"               # hexagon fill; the tree's knockout colour

## the Tree of Life ----------------------------------------------------------
#  The source is black line art whose thick strokes are drawn twice: a wide
#  black path, then a narrower white path on top.  The white pass is a knockout,
#  so it has to become the hexagon fill, not a lighter ink -- otherwise the tree
#  paints paper-coloured boxes over everything behind it.
tree <- paste(readLines(tree_svg, warn = FALSE), collapse = "\n")
tree <- sub("(?s)^.*?<svg[^>]*>", "", tree, perl = TRUE)
tree <- sub("(?s)</svg>\\s*$", "", tree, perl = TRUE)
tree <- sub("(?s)<title>.*?</title>", "", tree, perl = TRUE)
tree <- gsub("#000000", ink, tree, fixed = TRUE)
tree <- gsub("#FFFFFF", paper, tree, fixed = TRUE)
tree <- trimws(tree)

#  Drawn extent of the source, in its own svg units (see .tmp/tree-bbox.py):
#  x 7.9..308.1, y 4.6..595.4, centred on (158, 300).  Keep the aspect.
tree_h  <- 138
tree_w  <- tree_h * (300.3 / 590.8)
#  The tree sits half a sephirah below the sticker centre: that is enough to take
#  the lower crossing pair off Malkuth's name, while the disc still clears the
#  wordmark cap by 14 units.
tree_cy <- 119
tree_tf <- sprintf(
  "translate(%.2f,%.2f)scale(%.6f)translate(-158,-300)",
  cx, tree_cy, tree_h / 590.8
)

## the helix ----------------------------------------------------------------
#  A helix is recognisable only if its projected aspect is about right: the ratio
#  of pitch to diameter sets the angle at which the strands cross the axis, and
#  that angle is what says "double helix" rather than "stack of flat loops".
#  B-DNA (3.4 nm pitch, 2 nm diameter) crosses at 28 degrees above horizontal,
#  i.e. |dx/dy| = 2*pi*amp*turns/span = 1.85 at the crossing.
#
#  At a true 10 nm scale the duplex would be 2.86 turns wide 2 nm, i.e. amp ~14:
#  a narrow ribbon that no longer winds around the Tree.  The composition needs
#  the strands out at amp 48, just clear of the sephirot (x 85..155) -- and at
#  that width 2.86 turns would cross at 8 degrees, flat and ribbon-like.  Cutting
#  to 1.65 turns buys 15.6 degrees at the same width.  So both dimensions are
#  exaggerated (width ~3.4x, pitch ~1.7x); the projected aspect goes from
#  pitch/amp 0.85 to 1.72, against 3.40 for B-DNA.
turns <- 1.65
y_top <- 46; y_bot <- 182          # helix vertical span; clears the wordmark
amp   <- 48                        # half-width of the projected cylinder
depth <- 0.10                      # perspective: the near side reads closer
#  Groove asymmetry.  The two backbones are not diametrically opposite around the
#  cylinder: at a base pair they subtend ~140 degrees across the minor groove and
#  ~220 across the major one -- that is where B-DNA's 12 A / 22 A groove widths
#  come from at a ~9.4 A phosphate radius.  Strand two therefore sits 40 degrees
#  off the half-turn (220 - 180).  Height is proportional to azimuth here, so
#  that phase is an axial stagger of 40/360 of the pitch -- about one base pair
#  (0.34 nm), which is how far the two glycosidic bonds of a pair really are
#  apart along the axis.  A full half-turn (0.5) would draw equal grooves.
groove <- 40 / 360                 # strand offset, as a fraction of the pitch
strand_w <- 6                      # backbone stroke width
n     <- 1200

pitch <- (y_bot - y_top) / turns
dy    <- groove * pitch            # axial offset between the two backbones

tt <- seq(0, 1, length.out = n)
yy <- y_top + tt * (y_bot - y_top)
th <- 2 * pi * turns * tt
#  The sine is negated so the near half-turn sits left of the axis: the front
#  strokes then descend to the left, i.e. climb to the right when read upwards,
#  the sense of a right-handed (B-DNA) screw.  Mirroring is also what decides
#  which side of the diagram each backbone passes in front of.
x1 <- cx - amp * sin(th) * (1 + depth * cos(th))
x2 <- cx + amp * sin(th) * (1 - depth * cos(th))
# The backbones leave a base pair at slightly different heights, strand one a
# little above strand two.  That stagger is what opens the alternating wide
# (major) and narrow (minor) grooves, so it scales with the pitch.
y1 <- yy - dy / 2
y2 <- yy + dy / 2
#  The backbones sit at opposite points of the same cylinder, so their depths
#  are opposite too.  Sharing one depth sign would put the two strands in front
#  over the same stretch, and on that stretch the mirrored one winds the wrong
#  way, i.e. it reads as a left-handed helix.
z1 <- cos(th)                      # > 0: strand one is in front
z2 <- -z1                          # > 0: strand two is in front

# One run per contiguous stretch of a strand that shares a depth sign.  Each
# run reaches one sample past its boundary so the two halves meet at the turn.
runs <- function(x, y, keep) {
  grp <- rle(keep)
  ends <- cumsum(grp$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  out <- character(0)
  for (i in which(grp$values)) {
    lo <- max(1L, starts[i] - 1L)
    hi <- min(length(x), ends[i] + 1L)
    out <- c(out, sprintf('<polyline points="%s"/>',
                          paste(sprintf("%.1f,%.1f", x[lo:hi], y[lo:hi]), collapse = " ")))
  }
  out
}
# The two backbones carry different colours, and every part is translucent -- the
# front stretches most opaque, the parts that pass behind the diagram faintest --
# so the helix reads as glass wound around the diagram rather than ink on top.
#
# Both ends of the helix land on the axis, i.e. inside the topmost sefirah,
# whose disc carries its Hebrew name.  The ends of a front run are therefore cut
# where they clear that disc (plus the round cap), so the strands read as
# coming from behind the diagram rather than painted over the name.
top_disc  <- 100 * 0.392 * tree_h / 590.8   # Keter's radius, in sticker units
cap_clear <- top_disc + strand_w / 2
strand <- function(x, y, front, col, zs) {
  keep <- if (front) zs >= 0 else zs < 0
  ends <- tt < 0.02 | tt > 0.98
  keep <- keep & !(ends & abs(x - cx) < cap_clear)
  sprintf('<g fill="none" stroke="%s" stroke-width="%d" stroke-linecap="round" opacity="%.2f">%s</g>',
          col, strand_w, if (front) 0.7 else 0.3,
          paste(runs(x, y, keep), collapse = ""))
}
strand_one <- function(front) strand(x1, y1, front, accent, z1)
strand_two <- function(front) strand(x2, y2, front, ink, z2)

svg <- c(
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.0f %.0f" width="%.0f" height="%.0f" role="img" aria-label="tidymarkers: a DNA helix wound around the Tree of Life">', W, H, W, H),
  '<style>text{font-family:Inter,Helvetica,Arial,sans-serif;}</style>',
  sprintf('<polygon points="%s" fill="%s" stroke="%s" stroke-width="7" stroke-linejoin="round"/>', hex_pts, paper, ink),
  strand_one(FALSE),
  strand_two(FALSE),
  sprintf('<g transform="%s">%s</g>', tree_tf, tree),
  strand_one(TRUE),
  strand_two(TRUE),
  sprintf('<text x="%.0f" y="%.0f" text-anchor="middle" font-size="22" font-weight="700"><tspan fill="%s">tidy</tspan><tspan fill="%s">markers</tspan></text>', cx, cy + 71.5, ink, accent),
  '</svg>')
writeLines(svg, out_svg)
message("wrote ", out_svg)

## optional raster ----------------------------------------------------------
#  Prefer ImageMagick (7 renamed the CLI to `magick`; on Windows a bare
#  `convert` is the System32 disk tool and must never be used), fall back to
#  librsvg's rsvg-convert, which ships with most desktops.
raster <- function(src, dest, w, h) {
  convert <- Sys.which("convert")
  if (nzchar(convert) && !grepl("system32[/\\\\]convert\\.exe$", convert, ignore.case = TRUE)) {
    return(system2(unname(convert), c("-background", "none", "-density", "300", src, dest)))
  }
  magick <- Sys.which("magick")
  if (nzchar(magick)) {
    return(system2(unname(magick), c("convert", "-background", "none", "-density", "300", src, dest)))
  }
  rsvg <- Sys.which("rsvg-convert")
  if (nzchar(rsvg)) {
    return(system2(unname(rsvg), c("-w", w, "-h", h, "-b", "none", src, "-o", dest)))
  }
  NA_integer_
}
st <- raster(out_svg, out_png, "1500", "1732")
if (isTRUE(st == 0L)) message("wrote ", out_png) else message("PNG raster skipped (no ImageMagick or rsvg-convert)")
