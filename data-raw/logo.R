# =============================================================================
#  data-raw/logo.R  --  generate the package hex logo (tidyverse paradigm)
# =============================================================================
#  A tidyverse hex sticker: a hexagon with a border colour and an image inside.
#  The image is a KASP idea: a DNA double helix whose "rungs" are coloured like
#  the four bases, with the family read-out (the allelic-discrimination clusters
#  FAM/HEX) suggested by the palette.  Written as plain SVG (no graphics device,
#  no cairo needed) and optionally rasterised to PNG with ImageMagick.
#
#  Run:  Rscript data-raw/logo.R
#  Out:  man/figures/logo.svg  (+ logo.png if ImageMagick is available)
# =============================================================================

out_svg <- here::here("man", "figures", "logo.svg")
out_png <- here::here("man", "figures", "logo.png")
fs::dir_create(fs::path_dir(out_svg))

## geometry ------------------------------------------------------------------
W <- 240; H <- 277                 # canvas (tidyverse hex aspect)
cx <- W / 2; cy <- H / 2; r <- 116 # hexagon centre + circumradius
ang <- (90 + 0:5 * 60) * pi / 180  # pointy-top hexagon
hexx <- cx + r * cos(ang)
hexy <- cy + r * sin(ang)
hex_pts <- paste(sprintf("%.1f,%.1f", hexx, hexy), collapse = " ")

## colours (base palette ~ KASP FAM/HEX + a tidy accent) ----------------------
ink     <- "#1b4965"               # border / text
accent  <- "#d1495b"               # "kasp" accent
bases   <- c(A = "#5fa8d3", C = "#8ac926", G = "#ffb703", T = "#d1495b")

## the helix ----------------------------------------------------------------
n  <- 400
tt <- seq(0, 1, length.out = n)
y0 <- cy - 84; y1 <- cy + 50          # helix vertical span
yy <- y0 + tt * (y1 - y0)
amp <- 34; turns <- 2.4
xL <- cx + amp * sin(tt * 2 * pi * turns)
xR <- cx - amp * sin(tt * 2 * pi * turns)
poly <- function(x, y) sprintf("%.1f,%.1f", x, y) |> paste(collapse = " ")

# rungs every ~12 points, coloured by a pseudo base
ri   <- seq(20, n - 20, by = 12)
base <- rep(c("A", "C", "G", "T"), length.out = length(ri))
rungs <- sprintf(
  '<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="5" stroke-linecap="round"/>',
  xL[ri], yy[ri], xR[ri], yy[ri], unname(bases[base]))

svg <- c(
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.0f %.0f" width="%.0f" height="%.0f" role="img" aria-label="tidykasp">', W, H, W, H),
  '<style>text{font-family:Inter,Helvetica,Arial,sans-serif;}</style>',
  sprintf('<polygon points="%s" fill="#fffdf7" stroke="%s" stroke-width="7" stroke-linejoin="round"/>', hex_pts, ink),
  paste0(rungs, collapse = ""),
  sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="6" stroke-linecap="round"/>', poly(xL, yy), ink),
  sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="6" stroke-linecap="round"/>', poly(xR, yy), accent),
  sprintf('<text x="%.0f" y="%.0f" text-anchor="middle" font-size="30" font-weight="700"><tspan fill="%s">tidy</tspan><tspan fill="%s">kasp</tspan></text>', cx, cy + 92, ink, accent),
  '</svg>')
writeLines(svg, out_svg)
message("wrote ", out_svg)

## optional raster (ImageMagick) ---------------------------------------------
if (nzchar(Sys.which("convert"))) {
  st <- system2("convert", c("-background", "none", "-density", "300", out_svg, out_png))
  if (st == 0L) message("wrote ", out_png) else message("PNG raster skipped (convert failed)")
}
