// elk-layout.js -- ELK layered layout over a decision DAG, in two phases.
//
// in : {nodes:[{id,width,height,stage}], edges:[{id,from,to}]}   (points)
// out: {width,height,nodes,edges,stages}                          (inches, y DOWN)
//
// Phase 1 lays each stage out on its own with elk "layered" and direction DOWN:
// a stage is a gate chain (~6-23 nodes), so it wants a tall narrow band, and
// laying it in isolation keeps every intra-stage edge routable (the same edges
// inside a cross-stage hierarchy come back with zero sections).
//
// Phase 2 treats each solved stage band as one leaf box and lays those out with
// elk "layered" and direction LEFT: the stages read right to left at the call
// site, wide and short, which is the shape a README-width figure needs.
//
// Coordinates are assembled by hand from the two phases; elk never sees the
// hierarchy, so nothing can silently collapse two stages onto one layer.
"use strict";
const ELK = require("elkjs/lib/elk.bundled.js");
const SCALE = 72;

function readStdin() {
  return new Promise((resolve, reject) => {
    let buf = "";
    process.stdin.setEncoding("utf8");
    process.stdin.on("data", (c) => (buf += c));
    process.stdin.on("end", () => resolve(buf));
    process.stdin.on("error", reject);
  });
}

const STAGE_OPTIONS = {
  "elk.algorithm": "layered",
  "elk.direction": "DOWN",
  "elk.edgeRouting": "ORTHOGONAL",
  "elk.spacing.nodeNode": "8",
  "elk.layered.spacing.nodeNodeBetweenLayers": "12",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
  "elk.layered.cycleBreaking.strategy": "GREEDY",
  "elk.padding": "[top=4,left=4,bottom=4,right=4]",
  // sibling nodes (the alternate gate arms) may sit side by side instead of
  // being forced into one column, which is what keeps a band narrow
  "elk.portConstraints": "FREE",
  "elk.separateConnectedComponents": "false"
};

const ROOT_OPTIONS = {
  "elk.algorithm": "layered",
  "elk.direction": "LEFT",
  "elk.edgeRouting": "ORTHOGONAL",
  "elk.spacing.nodeNode": "6",
  "elk.layered.spacing.nodeNodeBetweenLayers": "18",
  "elk.layered.nodePlacement.strategy": "NETWORK_SIMPLEX",
  "elk.layered.crossingMinimization.strategy": "LAYER_SWEEP",
  "elk.layered.cycleBreaking.strategy": "DEPTH_FIRST",
  "elk.padding": "[top=4,left=4,bottom=4,right=4]",
  "elk.portConstraints": "FREE"
};

// Every section of every edge of a laid-out graph, in absolute points.
function collectEdges(node, ox, oy, out) {
  const x = ox + (node.x || 0);
  const y = oy + (node.y || 0);
  for (const c of node.children || []) collectEdges(c, x, y, out);
  for (const e of node.edges || []) {
    const pts = [];
    for (const s of e.sections || []) {
      if (s.startPoint) pts.push([s.startPoint.x + x, s.startPoint.y + y]);
      for (const b of s.bendPoints || []) pts.push([b.x + x, b.y + y]);
      if (s.endPoint) pts.push([s.endPoint.x + x, s.endPoint.y + y]);
    }
    if (pts.length) out.push({ id: e.id, from: e.sources[0], to: e.targets[0], points: pts });
  }
  return out;
}

function boxOf(n, ox, oy) {
  return { id: n.id, x: ox + (n.x || 0), y: oy + (n.y || 0), w: n.width, h: n.height };
}

(async () => {
  const g = JSON.parse(await readStdin());
  const elk = new ELK();

  // ---- phase 1: each stage, laid out on its own, top to bottom --------------
  const stageIds = [];
  const byStage = new Map();
  for (const n of g.nodes) {
    // ids may contain : (e.g. "as:gate:x"); never re-parse them
    const s = String(n.stage);
    if (!byStage.has(s)) { byStage.set(s, []); stageIds.push(s); }
    byStage.get(s).push(n);
  }
  const stageOf = new Map(g.nodes.map((n) => [n.id, String(n.stage)]));

  const bands = [];
  let cursorX = 0;
  for (const s of stageIds) {
    const members = byStage.get(s);
    const memberIds = new Set(members.map((n) => n.id));
    const intra = g.edges.filter((e) => memberIds.has(e.from) && memberIds.has(e.to));
    const sub = {
      id: "band",
      layoutOptions: STAGE_OPTIONS,
      children: members.map((n) => ({ id: n.id, width: n.width, height: n.height })),
      // `sources`/`targets` must be arrays: a bare string is dropped silently
      // by the JSON importer and the edge then fails "at least one source"
      edges: intra.map((e, i) => ({
        id: e.id || "i" + i + "::" + e.from + "->" + e.to,
        sources: [e.from],
        targets: [e.to]
      }))
    };
    const out = await elk.layout(sub);
    // Band padding.  The band must hold room for the stage label drawn above
    // it and for the routes that leave one band and enter the next, but a
    // large pad is what made the page 21.8 in wide (the bands are leaf boxes
    // in phase 2, so every point of pad multiplies by five).  24 pt keeps the
    // label clear of the first box and the page near square.
    const pad = Number(g.pad != null ? g.pad : 24);
    const w = Math.max(out.width || 0, ...members.map((n) => n.width)) + pad;
    const h = Math.max(out.height || 0, ...members.map((n) => n.height)) + pad;
    const bandX = cursorX;
    bands.push({
      id: s,
      x: bandX,
      y: 0,
      w,
      h,
      nodes: (out.children || []).map((c) => boxOf(c, bandX, 0)),
      edges: collectEdges(out, bandX, 0, [])
    });
    cursorX += w + 12; // provisional; phase 2 only uses the sizes
  }

  // ---- phase 2: the stage bands as leaf boxes, right to left ----------------
  // One root edge per *pair* of stages: a solver routes one line per edge, and
  // the per-node edges of a gate chain would otherwise all be drawn twice.
  const stageEdges = [];
  const seen = new Set();
  for (const e of g.edges) {
    const a = stageOf.get(e.from);
    const b = stageOf.get(e.to);
    if (!a || !b || a === b) continue;
    const k = a + ">" + b;
    if (seen.has(k)) continue;
    seen.add(k);
    stageEdges.push({ id: "s" + stageEdges.length, from: a, to: b });
  }
  const root = {
    id: "root",
    layoutOptions: ROOT_OPTIONS,
    children: bands.map((b) => ({ id: "B::" + b.id, width: b.w, height: b.h })),
    edges: stageEdges.map((e) => ({ id: e.id, sources: ["B::" + e.from], targets: ["B::" + e.to] }))
  };
  const laid = await elk.layout(root);
  const rootEdges = collectEdges(laid, 0, 0, []);
  const placed = new Map((laid.children || []).map((c) => [String(c.id).slice(3), c]));

  const nodes = [];
  const edges = [];
  for (const b of bands) {
    const p = placed.get(b.id) || { x: 0, y: 0 };
    const dx = (p.x || 0) - b.x;
    const dy = (p.y || 0) - b.y;
    for (const n of b.nodes) {
      nodes.push({ id: n.id, x: (n.x + dx) / SCALE, y: (n.y + dy) / SCALE });
    }
    for (const e of b.edges) {
      edges.push({
        id: e.id,
        from: e.from,
        to: e.to,
        points: e.points.map((q) => [(q[0] + dx) / SCALE, (q[1] + dy) / SCALE]),
        stage: b.id
      });
    }
  }
  // a root route joins the two band borders; carry it as a stage-level edge
  for (const e of rootEdges) {
    edges.push({
      id: e.id,
      from: e.from,
      to: e.to,
      points: e.points.map((q) => [q[0] / SCALE, q[1] / SCALE])
    });
  }
  const stages = bands.map((b) => {
    const p = placed.get(b.id) || { x: 0, y: 0 };
    return { id: b.id, x: (p.x || 0) / SCALE, y: (p.y || 0) / SCALE, w: b.w / SCALE, h: b.h / SCALE };
  });

  process.stdout.write(JSON.stringify({
    width: (laid.width || 0) / SCALE,
    height: (laid.height || 0) / SCALE,
    nodes,
    edges,
    stages
  }));
})().catch((e) => {
  process.stderr.write(String((e && e.stack) || e) + "\n");
  process.exit(1);
});
