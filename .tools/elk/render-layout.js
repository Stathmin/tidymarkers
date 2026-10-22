// render-layout.js -- structure in, placement out.  ELK places, libavoid routes.
//
// in  : {nodes:[{id,label,group,kind}], edges:[{from,to,kind}],
//        style:{font_file,font_size,wrap,max_width,spacing,...}}       (points)
// out : {width,height,nodes:[{id,x,y,w,h}],edges:[{from,to,points}],
//        groups:[{id,x,y,w,h}]}                                        (points, y DOWN)
//
// This process owns *every* coordinate: text is measured with opentype.js from
// a real font file, the boxes are placed by elk "layered" (two phases, see
// below), and the edges are routed by libavoid around those boxes.  Nothing
// upstream of it knows a position, and the R side only formats what comes back.
//
// Text is wrapped here, not in R, so the line breaks and the measured width come
// from one measurement: a label can never be wrapped to a width the box does not
// reserve.
"use strict";
const fs = require("fs");
const ELK = require("elkjs/lib/elk.bundled.js");
const libavoid = require("@mr_mint/elkjs-libavoid");
const opentype = require("opentype.js");

function readStdin() {
  return new Promise((resolve, reject) => {
    let buf = "";
    process.stdin.setEncoding("utf8");
    process.stdin.on("data", (c) => (buf += c));
    process.stdin.on("end", () => resolve(buf));
    process.stdin.on("error", reject);
  });
}

// ---- text measurement -------------------------------------------------------
function loadFont(file) {
  const buf = fs.readFileSync(file);
  return opentype.parse(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
}

// width of one line in points, at `size` points (opentype's units are per em)
function lineWidth(font, text, size) {
  return font.getAdvanceWidth(String(text).replace(/\s+/g, " "), size);
}

// a bracketed tail such as "(n)" or "(gate)" is notation for the word before
// it: splitting between them leaves a stray "(", so they travel together
function wrapUnits(text) {
  const words = String(text).split(/[ ]+/).filter(Boolean);
  const units = [];
  for (const w of words) {
    const prev = units[units.length - 1];
    if (prev !== undefined && /^[(\[][^)\]]*[)\]]?$/.test(w)) {
      units[units.length - 1] = prev + " " + w;
    } else {
      units.push(w);
    }
  }
  return units;
}

function wrap(font, text, size, maxWidth) {
  const units = wrapUnits(text);
  if (!units.length) return [""];
  const n = units.length;
  const width = units.map((u) => lineWidth(font, u, size));
  const space = lineWidth(font, " ", size);
  const best = new Array(n + 1).fill(Infinity);
  const from = new Array(n + 1).fill(0);
  best[0] = 0;
  for (let j = 1; j <= n; j++) {
    let w = 0;
    for (let i = j - 1; i >= 0; i--) {
      w += width[i] + (i < j - 1 ? space : 0);
      if (w > maxWidth && j - i > 1) break;
      if (w > maxWidth) continue; // one unit wider than the box: unsplittable
      const slack = (maxWidth - w) / maxWidth;
      const cost = best[i] + slack * slack;
      if (cost < best[j]) {
        best[j] = cost;
        from[j] = i;
      }
    }
  }
  if (!isFinite(best[n])) return [units.join(" ")];
  const lines = [];
  for (let j = n; j > 0; j = from[j]) lines.unshift(units.slice(from[j], j).join(" "));
  return lines;
}

// ---- placement --------------------------------------------------------------
// Two phases, because a single hierarchy collapses (elk's layered algorithm
// puts stage containers on one layer when INCLUDE_CHILDREN is on, and drops the
// cross-stage edges when it is off):
//   phase 1  each group on its own, direction DOWN -- a gate chain is tall and
//            narrow, and laying it in isolation keeps every intra-group edge
//   phase 2  the solved bands, folded into rows; a row is a "DOWN" run with no
//            edges, so every band in it lands in one layer and runs across the
//            page.  Where the chain folds is chosen for A5's aspect.
const GROUP_OPTIONS = {
  "elk.algorithm": "layered",
  // top-to-bottom *within* a group: a gate chain reads as a chain
  "elk.direction": "DOWN",
  "elk.edgeRouting": "UNDEFINED",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
  "elk.layered.cycleBreaking.strategy": "GREEDY",
  "elk.layered.spacing.nodeNodeBetweenLayers": "10",
  "elk.aspectRatio": "0.95",
  "elk.layered.wrapping.strategy": "SINGLE_EDGE",
  "elk.layered.wrapping.cutting.strategy": "MSD",
  "elk.layered.wrapping.additionalEdgeSpacing": "12",
  "elk.spacing.nodeNode": "10",
  "elk.separateConnectedComponents": "false"
};

// a row is a "DOWN" run with no edges: every band lands in the one layer, and
// a layer runs across the page, which is exactly a row
const ROW_OPTIONS = {
  "elk.algorithm": "layered",
  "elk.direction": "DOWN",
  "elk.edgeRouting": "UNDEFINED",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
  "elk.layered.cycleBreaking.strategy": "DEPTH_FIRST",
  "elk.layered.spacing.nodeNodeBetweenLayers": "10",
  "elk.spacing.nodeNode": "10",
  "elk.separateConnectedComponents": "false"
};

function rect(x, y, w, h) {
  return { x: x || 0, y: y || 0, w: w || 0, h: h || 0 };
}

function grow(r, dx, dy) {
  return rect(r.x - dx, r.y - dy, r.w + 2 * dx, r.h + 2 * dy);
}

// The group bands are a chain, so the only placement decision left is where
// that chain folds into rows -- and that is what sets the page shape.  Ask for
// A5's aspect (420/595): a square page wastes a third of the sheet in print.
const A5 = 420 / 595;
const ROW_GAP = 10;

function foldRows(bands) {
  const n = bands.length;
  if (n < 2) return [bands];
  const options = [];
  for (let c1 = 1; c1 < n; c1++) {
    options.push([c1]);
    for (let c2 = c1 + 1; c2 < n; c2++) options.push([c1, c2]);
  }
  let best = null;
  for (const cuts of options) {
    const slices = [];
    let start = 0;
    for (const stop of [...cuts, n]) {
      slices.push(bands.slice(start, stop));
      start = stop;
    }
    const width = Math.max(...slices.map((r) => r.reduce((s, b) => s + b.w, 0) + ROW_GAP * (r.length - 1)));
    const height = slices.reduce((s, r) => s + Math.max(...r.map((b) => b.h)), 0) + ROW_GAP * (slices.length - 1);
    const score = Math.abs(width / height - A5);
    if (!best || score < best.score) best = { score, slices };
  }
  return best.slices;
}

async function place(elk, groups, edges) {
  // phase 1: a group on its own.  A gate chain is tall and narrow, so laying
  // it in isolation keeps every intra-group edge -- the bulk of the graph --
  // in one column with its wrap.
  const bands = [];
  for (const g of groups) {
    const memberIds = new Set(g.members.map((n) => n.id));
    const intra = edges
      .filter((e) => memberIds.has(e.from) && memberIds.has(e.to))
      .map((e, i) => ({ id: "i" + i + "|" + e.from + "|" + e.to, sources: [e.from], targets: [e.to] }));
    const laid = await elk.layout({
      id: "G::" + g.id,
      layoutOptions: GROUP_OPTIONS,
      children: g.members.map((n) => ({ id: n.id, width: n.w, height: n.h })),
      edges: intra
    });
    bands.push({ id: g.id, w: laid.width || 0, h: laid.height || 0, elk: laid, x: 0, y: 0 });
  }

  // phase 2: each row is its own elk run, so the order inside a row is the
  // solver's; only stacking one row under the other is arithmetic
  let y = 0;
  for (const row of foldRows(bands)) {
    const laid = await elk.layout({
      id: "row" + y,
      layoutOptions: ROW_OPTIONS,
      children: row.map((b) => ({ id: b.id, width: b.w, height: b.h })),
      edges: []
    });
    const byId = new Map((laid.children || []).map((c) => [c.id, c]));
    const rowH = laid.height || Math.max(...row.map((b) => b.h));
    for (const b of row) {
      const c = byId.get(b.id);
      b.x = Number.isFinite(c && c.x) ? c.x : 0;
      b.y = y + (Number.isFinite(c && c.y) ? c.y : 0);
    }
    y += rowH + ROW_GAP;
  }
  return bands;
}

(async () => {
  // the request arrives as a path argument (stdin when none): a JSON document
  // can exceed a pipe buffer, and a file path keeps the write side dumb
  const argFile = process.argv[2];
  const req = JSON.parse(argFile ? fs.readFileSync(argFile, "utf8") : await readStdin());
  const st = req.style || {};
  const size = Number(st.font_size || 9);
  const lineH = size * Number(st.line_height || 1.30);
  const padX = Number(st.pad_x || 7);
  const padY = Number(st.pad_y || 5);
  const fontFile = st.font_file || "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf";
  const boldFile = st.bold_file || fontFile;
  const maxText = Number(st.max_width || 96);
  const font = loadFont(fontFile);
  const bold = loadFont(boldFile);

  // ---- structure -> labelled boxes -----------------------------------------
  const nodes = new Map();
  for (const n of req.nodes) {
    const f = n.bold ? bold : font;
    const lines = wrap(f, n.label, size, maxText);
    const widths = lines.map((l) => lineWidth(f, l, size));
    const w = Math.max(...widths) + 2 * padX;
    const h = lines.length * lineH + 2 * padY;
    nodes.set(n.id, {
      id: n.id, lines, w, h, group: String(n.group), kind: n.kind,
      status: n.status || "neutral", bold: n.bold === true
    });
  }
  const edges = req.edges
    .filter((e) => nodes.has(e.from) && nodes.has(e.to) && e.from !== e.to)
    .map((e) => ({ from: e.from, to: e.to, kind: e.kind }));

  // groups in the order they were declared (the pipeline order); a group with
  // no member is dropped (an empty column would be a bare gap)
  const order = [];
  const byGroup = new Map();
  for (const n of nodes.values()) {
    if (!byGroup.has(n.group)) {
      byGroup.set(n.group, []);
      order.push(n.group);
    }
    byGroup.get(n.group).push(n);
  }
  const groups = order.map((id) => ({ id, members: byGroup.get(id) }));

  const elk = new ELK();
  const bands = await place(elk, groups, edges);

  // ---- phase 1 output -> absolute boxes ------------------------------------
  const placedGroups = new Map();
  for (const b of bands) {
    placedGroups.set(b.id, { x: b.x, y: b.y, w: b.w, h: b.h, elk: b.elk });
  }
  const placedNodes = [];
  for (const g of groups) {
    const pg = placedGroups.get(g.id);
    if (!pg) continue;
    const byId = new Map((pg.elk.children || []).map((c) => [c.id, c]));
    for (const m of g.members) {
      const c = byId.get(m.id);
      if (!c) continue;
      placedNodes.push({
        id: m.id,
        x: pg.x + c.x,
        y: pg.y + c.y,
        w: m.w,
        h: m.h,
        lines: m.lines,
        kind: m.kind,
        status: m.status,
        bold: m.bold
      });
    }
  }
  // elk sized the group from the boxes it placed; keep its extent (it already
  // includes the between-group spacing the root layout used)
  const groupBoxes = groups
    .filter((g) => placedGroups.has(g.id))
    .map((g) => {
      const pg = placedGroups.get(g.id);
      // a group placed at the origin comes back with null x/y from elk
      return { id: g.id, x: pg.x || 0, y: pg.y || 0, w: pg.w, h: pg.h };
    });

  // ---- ports ---------------------------------------------------------------
  // one pin at the middle of each side of every box, and every edge picks its
  // pair from the *placement*, never from libavoid geometry.  That is what
  // makes like edges look alike: all edges that leave the same side of a box
  // leave from the same point, so their first segment is a shared rail.
  //   same group, same column   -> S out, N in   (a group reads downwards)
  //   same group, other column  -> E out, W in   (a wrap continues rightward)
  //   different groups          -> the dominant axis of the two band centres
  const groupOfNode = new Map();
  for (const g of groups) for (const m of g.members) groupOfNode.set(m.id, g.id);
  const nodeById = new Map(placedNodes.map((n) => [n.id, n]));
  const boxById = new Map(groupBoxes.map((g) => [g.id, g]));

  // same column = the two centres sit inside half a box of each other.  A raw
  // centre index would split a column the moment two labels differ in width.
  // the x a vertical run uses inside a box.  Boxes in one column have
  // different widths, so their geometric centres differ by a point or two and
  // each chain edge would get a micro-jog; pinning N/S to the column axis
  // makes every edge in a column one straight line.
  const axisFor = new Map();
  {
    const perGroup = new Map();
    for (const n of placedNodes) {
      const gid = groupOfNode.get(n.id);
      if (!perGroup.has(gid)) perGroup.set(gid, []);
      perGroup.get(gid).push(n);
    }
    for (const list of perGroup.values()) {
      const sorted = [...list].sort((a, b) => (a.x + a.w / 2) - (b.x + b.w / 2));
      const clusters = [];
      for (const n of sorted) {
        const c = n.x + n.w / 2;
        const last = clusters[clusters.length - 1];
        if (last && c - last.last <= 24) {
          last.sum += c;
          last.count++;
          last.last = c;
          last.members.push(n);
        } else {
          clusters.push({ sum: c, count: 1, last: c, members: [n] });
        }
      }
      for (const cl of clusters) {
        const axis = cl.sum / cl.count;
        for (const n of cl.members) {
          axisFor.set(n.id, Math.max(2, Math.min(n.w - 2, axis - n.x)));
        }
      }
    }
  }

  const SIDE = {
    N: (n) => [axisFor.get(n.id), 0, "NORTH"],
    S: (n) => [axisFor.get(n.id), n.h, "SOUTH"],
    W: (n) => [0, n.h / 2, "WEST"],
    E: (n) => [n.w, n.h / 2, "EAST"]
  };
  const portsFor = (n) => Object.keys(SIDE).map((k) => {
    const [x, y, side] = SIDE[k](n);
    return { id: n.id + "#" + k, x, y, width: 0, height: 0, properties: { "port.side": side } };
  });

  function sides(a, b) {
    const ga = groupOfNode.get(a.id);
    const gb = groupOfNode.get(b.id);
    if (ga === gb) {
      const dx = (b.x + b.w / 2) - (a.x + a.w / 2);
      if (Math.abs(dx) < 0.5 * Math.min(a.w, b.w)) return ["S", "N"];
      return dx > 0 ? ["E", "W"] : ["W", "E"];
    }
    const A = boxById.get(ga);
    const B = boxById.get(gb);
    const dx = (B.x + B.w / 2) - (A.x + A.w / 2);
    const dy = (B.y + B.h / 2) - (A.y + A.h / 2);
    if (Math.abs(dx) >= Math.abs(dy)) return dx >= 0 ? ["E", "W"] : ["W", "E"];
    return dy >= 0 ? ["S", "N"] : ["N", "S"];
  }

  // ---- routing: libavoid, around the boxes ---------------------------------
  // `routeEdges` reads absolute coordinates off the graph; a flat graph of
  // boxes, their pins and their edges is the cleanest thing to hand it.
  const routeGraph = {
    id: "route",
    children: placedNodes.map((n) => ({
      id: n.id, x: n.x, y: n.y, width: n.w, height: n.h, ports: portsFor(n)
    })),
    edges: edges.map((e, i) => {
      const [so, ti] = sides(nodeById.get(e.from), nodeById.get(e.to));
      return {
        id: "e" + i,
        sources: [e.from + "#" + so],
        targets: [e.to + "#" + ti]
      };
    })
  };
  const routes = await libavoid.routeEdges(routeGraph, {
    routingType: "orthogonal",
    segmentPenalty: Number(st.segment_penalty || 8),
    crossingPenalty: Number(st.crossing_penalty || 20),
    shapeBufferDistance: Number(st.shape_buffer || 3),
    idealNudgingDistance: Number(st.nudge || 10),
    portDirectionPenalty: 100,
    nudgeOrthogonalSegmentsConnectedToShapes: true,
    nudgeOrthogonalTouchingColinearSegments: true,
    performUnifyingNudgingPreprocessingStep: true,
    selfLoopHandling: "skip"
  });

  // the router emits near-duplicate points (0.0-0.4 pt stubs) where a pin's
  // exit meets its bend, and those read as wobble.  Collapse sub-half-point
  // moves and drop corners that do not actually turn.
  function tidy(points) {
    const out = [];
    for (const p of points) {
      const q = out[out.length - 1];
      if (q && Math.abs(p.x - q.x) < 0.5 && Math.abs(p.y - q.y) < 0.5) continue;
      out.push(p);
    }
    const keep = [];
    for (let i = 0; i < out.length; i++) {
      const a = keep[keep.length - 1];
      const c = out[i + 1];
      if (a && c && (
        (Math.abs(a.x - out[i].x) < 0.5 && Math.abs(c.x - out[i].x) < 0.5) ||
        (Math.abs(a.y - out[i].y) < 0.5 && Math.abs(c.y - out[i].y) < 0.5)
      )) continue;
      keep.push(out[i]);
    }
    return keep.length >= 2 ? keep : out;
  }

  const outEdges = [];
  edges.forEach((e, i) => {
    const r = routes.get("e" + i);
    if (!r) return;
    const pts = tidy([r.sourcePoint, ...(r.bendPoints || []), r.targetPoint]
      .filter((p) => p && Number.isFinite(p.x) && Number.isFinite(p.y)))
      .map((p) => [p.x, p.y]);
    outEdges.push({ from: e.from, to: e.to, kind: e.kind, points: pts });
  });

  // ---- extent --------------------------------------------------------------
  // the page holds the boxes, the routes and a little air; libavoid may run a
  // route beside its boxes, so the routes count too
  const margin = Number(st.margin || 10);
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const n of placedNodes) {
    if (![n.x, n.y, n.w, n.h].every(Number.isFinite)) continue;
    x0 = Math.min(x0, n.x); y0 = Math.min(y0, n.y);
    x1 = Math.max(x1, n.x + n.w); y1 = Math.max(y1, n.y + n.h);
  }
  // an unrouted edge has no points; a null in the list would poison every min
  for (const e of outEdges) for (const p of e.points) {
    if (!p || !Number.isFinite(p[0]) || !Number.isFinite(p[1])) continue;
    x0 = Math.min(x0, p[0]); y0 = Math.min(y0, p[1]);
    x1 = Math.max(x1, p[0]); y1 = Math.max(y1, p[1]);
  }
  // group bands are drawn too, so they bound the page like a node does
  for (const g of groupBoxes) {
    if (![g.x, g.y, g.w, g.h].every(Number.isFinite)) continue;
    x0 = Math.min(x0, g.x); y0 = Math.min(y0, g.y);
    x1 = Math.max(x1, g.x + g.w); y1 = Math.max(y1, g.y + g.h);
  }
  // group origins come from elk as null when a group sits at 0; only shift real ones
  const shift = (v, o) => v - o;
  const nodesOut = placedNodes.map((n) => ({
    id: n.id, x: shift(n.x, x0 - margin), y: shift(n.y, y0 - margin),
    w: n.w, h: n.h, lines: n.lines, kind: n.kind, status: n.status, bold: n.bold
  }));
  const edgesOut = outEdges.map((e) => ({
    from: e.from, to: e.to, kind: e.kind,
    points: e.points.map((p) => [shift(p[0], x0 - margin), shift(p[1], y0 - margin)])
  }));
  const groupsOut = groupBoxes.map((g) => ({
    id: g.id, x: shift(g.x, x0 - margin), y: shift(g.y, y0 - margin),
    w: g.w, h: g.h
  }));

  process.stdout.write(JSON.stringify({
    width: x1 - x0 + 2 * margin,
    height: y1 - y0 + 2 * margin,
    font_size: size,
    line_height: lineH,
    nodes: nodesOut,
    edges: edgesOut,
    groups: groupsOut
  }));
})().catch((e) => {
  process.stderr.write(String((e && e.stack) || e) + "\n");
  process.exit(1);
});
