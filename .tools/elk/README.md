# `.tools/elk` - the vendored decision-graph layout engine

Placement and edge routing for `plot_decision_tree()` / `decision_render()`
(`R/marker_decision_solver.R`). The engine is **vendored** in this repository:
a render must work offline and reproduce the committed figures exactly, so the
packages and their licences live in this tree.

| package | version | licence | role |
|---|---|---|---|
| `elkjs` | ^0.11.0 | EPL-2.0 | layered placement (two phases, see below) |
| `libavoid-js` | via `@mr_mint` | LGPL-2.1 / MIT | obstacle-avoiding orthogonal edge routing |
| `@mr_mint/elkjs-libavoid` | ^0.5.0 | MIT | the glue that feeds ELK boxes to libavoid |
| `opentype.js` | ^2.0.0 | MIT | measures label text from a real TTF |

Files: `render-layout.js` (the seam: one JSON document in - on stdin, or as a
path argument - and one out on stdout), `package.json`, `package-lock.json`,
`node_modules/` (licences ship inside each package).

## Why these

- **ELK, not Graphviz.** Graphviz's `dot` is a single ranking pass: no
  obstacle-avoiding routing, no control over the two-phase decomposition a stage
  band needs. ELK's layered algorithm exposes node placement
  (`NETWORK_SIMPLEX`), crossing minimisation (`LAYER_SWEEP`) and cycle breaking
  per phase.
- **libavoid, not ELK's router.** ELK's orthogonal router assumes the edges live
  inside the laid-out hierarchy; every attempt to keep a stage hierarchy either
  collapsed it (99 x 2 in) or lost the cross-stage edges (zero sections).
  libavoid routes a flat graph *around* rectangles, which is the real problem
  once the boxes are placed: it keeps a buffer from each box, nudges parallel
  segments apart, and emits orthogonal polylines.
- **Two phases.** Phase 1 lays each stage group out on its own
  (`direction=DOWN`: a gate chain is tall and narrow). Phase 2 treats those
  solved group boxes as leaves and lays them out `direction=RIGHT`: the stages
  read as columns and the page stays near square. A single hierarchy with
  `INCLUDE_CHILDREN` puts every stage on one layer; `SEPARATE_CHILDREN` loses the
  cross-stage edges.
- **Text from a real face.** `opentype.js` parses a TTF and measures each label
  at the figure's own point size, so a box is exactly as wide as the text it
  holds. Wrapping happens in the same process, against the same measurement.

## Rebuilding (only when bumping a version)

```bash
cd .tools/elk
npm install --no-audit --no-fund @mr_mint/elkjs-libavoid elkjs opentype.js
```

Then re-render and refresh the artifact hashes:

```bash
Rscript dev/render-decision.R
Rscript dev/artifacts.R --write && Rscript dev/check-artifacts.R
Rscript dev/test.R decision
```

`node` must be on `PATH` for a render. Without it, or with `TIDYKASP_NO_ELK=1`,
`decision_render()` falls back to a simple column placement and says so on the
console; the figure is still written.
