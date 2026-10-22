# Lecture: KASP marker design in R — an end-to-end walkthrough of the method

`kasp-slides.qmd` is a Quarto (reveal.js) deck. It is a **live run**: when it
renders, it actually executes the `tidymarkers` pipeline on the FRI locus
(*Arabidopsis thaliana*, `AT4G00650`) and shows the code, the table and the
number at every step.

## Presenting it

```bash
sh dev/lecture/render.sh   # = quarto render dev/lecture/kasp-slides.qmd + OligoArrayAux on PATH
# then open _site/dev/lecture/kasp-slides.html (self-contained, no network)
```

In the browser:

- `n` — **speaker notes** (almost every slide has `.notes`: what to say and what
  to point at);
- `s` — presenter mode (current + next slide + notes + timer);
- `f` — full screen; `o` — slide overview.

Also useful: `quarto preview dev/lecture/kasp-slides.qmd` (live update while you
edit the source) and `quarto render ... --to pdf` (a handout).

Requirements: the `quarto` CLI (>= 1.4), a restored `renv` plus the installed
package, `OligoArrayAux` on `PATH` (`render.sh` sets this up), and one network
fetch or warm cache (TAIR10 genome, Ensembl variants, FM index — all inside
`marker_cache_dir()`, never in the project).

## Structure (following the method's stages)

- **Stage 0.** What KASP is; what design has to solve.
- **Stage 1.** Target: variant quality, consequence, the 3' window
  (`pick_target_snp`).
- **Stage 2.** Design window and the strand hint (`strand_gc_hint`).
- **Stage 3.** The AS primer: the search (both strands × length × mismatch
  position), Tm/GC windows; why the chosen primer is not the "loudest" one
  (footprint matters more than raw discrimination).
- **Stage 4.** Discrimination and 3' stability (`DECIPHER`, `terminal_dg`).
- **Stage 5.** The common primer (Tm >= AS, amplicon <= 100 bp).
- **Stage 6.** Structural QC: hairpins, dimers, inverted repeats.
- **Stage 7.** Specificity: 3' seed (off-targets), in-silico PCR, **FM index**.
- **Stage 8.** Search animation (`design_search.gif`) — how the selection
  converges.
- **Stage 9.** Acceptance (14 criteria) and ordering (IDT/LGC).
- **Stage 10.** **The full method checklist** — every recommendation → the
  package function → our result (3 slides).

## Notes

- Heavy chunks are cached (`cache: true`), so a re-render is fast.
- Everything generated goes to `_site/` (HTML, chunk cache, figures); only the
  sources live in `dev/lecture/`. `_site/` and the cache are in `.gitignore`;
  `dev/` is in `.Rbuildignore`. The one exception is the search animation, which
  `render_design_search()` writes to `results/` — also ignored, and shared with
  `inst/scripts/run_kasp.R`.
- The data are real: FRI (`AT4G00650`), 1001 Genomes through the Ensembl REST
  API.
