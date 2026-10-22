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
  return font.getAdvanceWidth(text, size);
}

function wrap(font, text, size, maxWidth) {
  const words = String(text).split(/[ ]+/).filter(Boolean);
  if (!words.length) return [""];
  const lines = [];
  let line = words[0];
  for (const w of words.slice(1)) {
    const cand = line + " " + w;
    if (lineWidth(font, cand, size) <= maxWidth) line = cand;
    else {
      lines.push(line);
      line = w;
    }
  }
  lines.push(line);
  return lines;
}

// ---- placement --------------------------------------------------------------
// Two phases, because a single hierarchy collapses (elk's layered algorithm
// puts stage containers on one layer when INCLUDE_CHILDREN is on, and drops the
// cross-stage edges when it is off):
//   phase 1  each group on its own, direction DOWN -- a gate chain is tall and
//            narrow, and laying it in isolation keeps every intra-group edge
//   phase 2  the solved group boxes, direction LEFT -- the groups read as
//            columns, which is the shape a page-width figure needs
const GROUP_OPTIONS = {
  "elk.algorithm": "layered",
  // top-to-bottom *within* a group: a gate chain reads as a chain
  "elk.direction": "DOWN",
  "elk.edgeRouting": "UNDEFINED",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
  "elk.layered.cycleBreaking.strategy": "GREEDY",
  "elk.layered.spacing.nodeNodeBetweenLayers": "18",
  "elk.spacing.nodeNode": "10",
  "elk.separateConnectedComponents": "false"
};

const ROOT_OPTIONS = {
  "elk.algorithm": "layered",
  "elk.direction": "RIGHT",
  "elk.edgeRouting": "UNDEFINED",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
  "elk.layered.cycleBreaking.strategy": "DEPTH_FIRST",
  "elk.layered.spacing.nodeNodeBetweenLayers": "70",
  "elk.spacing.nodeNode": "28",
  "elk.separateConnectedComponents": "false"
};

function rect(x, y, w, h) {
  return { x: x || 0, y: y || 0, w: w || 0, h: h || 0 };
}

function grow(r, dx, dy) {
  return rect(r.x - dx, r.y - dy, r.w + 2 * dx, r.h + 2 * dy);
}

async function place(elk, groups, edges) {
  const elkGroups = groups.map((g) => {
    const memberIds = new Set(g.members.map((n) => n.id));
    const intra = edges
      .filter((e) => memberIds.has(e.from) && memberIds.has(e.to))
      .map((e, i) => ({ id: "i" + i + "|" + e.from + "|" + e.to, sources: [e.from], targets: [e.to] }));
    return {
      id: "G::" + g.id,
      layoutOptions: GROUP_OPTIONS,
      children: g.members.map((n) => ({ id: n.id, width: n.w, height: n.h })),
      edges: intra
    };
  });
  // the group graph: one edge per pair of groups.  The layered algorithm needs
  // *some* linkage to order the groups left to right, and per-node cross edges
  // would multiply by the number of aligned nodes.
  const groupOf = new Map();
  for (const g of groups) for (const n of g.members) groupOf.set(n.id, g.id);
  const pair = new Map();
  for (const e of edges) {
    const a = groupOf.get(e.from);
    const b = groupOf.get(e.to);
    if (!a || !b || a === b) continue;
    const k = a + "\u0000" + b;
    if (!pair.has(k)) pair.set(k, { id: "p" + pair.size, from: a, to: b });
  }
  const groupEdges = [...pair.values()];
  // the root order follows the crossing edges: elk has no order to preserve
  // once the groups are disconnected
  const root = {
    id: "root",
    layoutOptions: ROOT_OPTIONS,
    children: elkGroups,
    edges: groupEdges.map((e) => ({ id: e.id, sources: ["G::" + e.from], targets: ["G::" + e.to] }))
  };
  return elk.layout(root);
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
    const w = Math.max(...lines.map((l) => lineWidth(f, l, size))) + 2 * padX;
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
  const laid = await place(elk, groups, edges);

  // ---- phase 1 output -> absolute boxes ------------------------------------
  const placedGroups = new Map();
  for (const c of laid.children || []) {
    const id = String(c.id).slice(3);
    placedGroups.set(id, { x: c.x, y: c.y, w: c.width, h: c.height, elk: c });
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
      return { id: g.id, x: pg.x, y: pg.y, w: pg.w, h: pg.h };
    });

  // ---- routing: libavoid, around the boxes ---------------------------------
  // `routeEdges` reads absolute coordinates off the graph; a flat graph of
  // boxes and edges is the cleanest thing to hand it (nested graphs are routed
  // relative to their owner's content area, which is not what we need here).
  const routeGraph = {
    id: "route",
    children: placedNodes.map((n) => ({ id: n.id, x: n.x, y: n.y, width: n.w, height: n.h })),
    edges: edges.map((e, i) => ({ id: "e" + i, source: e.from, target: e.to }))
  };
  const routes = await libavoid.routeEdges(routeGraph, {
    routingType: "orthogonal",
    segmentPenalty: Number(st.segment_penalty || 24),
    crossingPenalty: Number(st.crossing_penalty || 40),
    shapeBufferDistance: Number(st.shape_buffer || 7),
    idealNudgingDistance: Number(st.nudge || 7),
    portDirectionPenalty: 100,
    selfLoopHandling: "skip"
  });

  const outEdges = [];
  edges.forEach((e, i) => {
    const r = routes.get("e" + i);
    if (!r) return;
    const pts = [r.sourcePoint, ...(r.bendPoints || []), r.targetPoint].map((p) => [p.x, p.y]);
    outEdges.push({ from: e.from, to: e.to, kind: e.kind, points: pts });
  });

  // ---- extent --------------------------------------------------------------
  // the page holds the boxes, the routes and a little air; libavoid may run a
  // route beside its boxes, so the routes count too
  const margin = Number(st.margin || 14);
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const n of placedNodes) {
    x0 = Math.min(x0, n.x); y0 = Math.min(y0, n.y);
    x1 = Math.max(x1, n.x + n.w); y1 = Math.max(y1, n.y + n.h);
  }
  for (const e of outEdges) for (const p of e.points) {
    x0 = Math.min(x0, p[0]); y0 = Math.min(y0, p[1]);
    x1 = Math.max(x1, p[0]); y1 = Math.max(y1, p[1]);
  }
  for (const g of groupBoxes) {
    x0 = Math.min(x0, g.x); y0 = Math.min(y0, g.y);
    x1 = Math.max(x1, g.x + g.w); y1 = Math.max(y1, g.y + g.h);
  }
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
