import {api, LigandDepiction} from '../api';
import {MARKOV_COLORS, MODE_OTHER, MODE_PALETTE} from '../colors';
import {clear, h, hideTooltip, icon, pct, s, showTooltip} from '../dom';
import {BindingMode, MarkovEdge, MarkovNode, ThirdColor} from '../derive';
import {svgExportButtons} from '../export';
import {drawCategories, frameAxis, positionAt} from './barcode';
import {FeatureView, ResultsContext} from './context';
import {depiction, ligandThumbnail} from './ligand2d';
import {thresholdControl} from './modes';

const W = 900;
const H = 620;
const MAX_NODES = 24;
const LAYOUT_KEY = 'openmmdl-analysis-gui:markov-layout';
const STRUCTURES_KEY = 'openmmdl-analysis-gui:markov-structures';

function stored(key: string): string | null {
  try {
    return localStorage.getItem(key);
  } catch {
    return null;
  }
}

function store(key: string, value: string): void {
  try {
    localStorage.setItem(key, value);
  } catch {
    // View preferences are a convenience only.
  }
}

type Layout = 'force' | 'circle';

interface NodeView {
  node: MarkovNode;
  mode: BindingMode | undefined;
  x: number;
  y: number;
  r: number;
  /** Half width/height of a structure card; null for a plain circle of radius r. */
  card: {hw: number; hh: number} | null;
  g: SVGGElement;
}

interface EdgeView {
  edge: MarkovEdge;
  g: SVGGElement;
  paths: SVGPathElement[];
  label: SVGTextElement;
}

export function modeColor(mode: BindingMode | undefined): string {
  return mode && mode.rank < MODE_PALETTE.length ? MODE_PALETTE[mode.rank] : MODE_OTHER;
}

function mix(hex: string, target: string, t: number): string {
  const a = parseInt(hex.slice(1), 16);
  const b = parseInt(target.slice(1), 16);
  const ch = (n: number, shift: number) => (n >> shift) & 255;
  const c = [16, 8, 0].map((sh) => Math.round(ch(a, sh) + (ch(b, sh) - ch(a, sh)) * t));
  return `#${c.map((v) => v.toString(16).padStart(2, '0')).join('')}`;
}

/** Mean number of consecutive frames a mode lasts once entered. */
function dwell(labels: Int32Array, mode: number): number {
  let runs = 0;
  let frames = 0;
  for (let i = 0; i < labels.length; i++) {
    if (labels[i] !== mode) continue;
    frames++;
    if (i === 0 || labels[i - 1] !== mode) runs++;
  }
  return runs ? frames / runs : 0;
}

/** Deterministic force-directed placement (Fruchterman-Reingold with collision and gravity). */
function forceLayout(views: NodeView[], edges: MarkovEdge[]): void {
  const n = views.length;
  if (n === 1) {
    views[0].x = W / 2;
    views[0].y = H / 2;
    return;
  }
  const index = new Map(views.map((v, i) => [v.node.mode, i]));
  views.forEach((v, i) => {
    const angle = -Math.PI / 2 + (2 * Math.PI * i) / n;
    const radius = Math.min(W, H) * (0.18 + 0.22 * (i / n));
    v.x = W / 2 + radius * Math.cos(angle);
    v.y = H / 2 + radius * Math.sin(angle);
  });
  const k = Math.sqrt((W * H) / n) * 0.38;
  const maxPct = Math.max(1e-9, ...edges.map((e) => e.framesPct));
  for (let it = 0; it < 600; it++) {
    const temp = (W / 12) * (1 - it / 600) + 0.3;
    const dx = new Float64Array(n);
    const dy = new Float64Array(n);
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) {
        const a = views[i];
        const b = views[j];
        let ex = a.x - b.x;
        let ey = a.y - b.y;
        let d = Math.hypot(ex, ey);
        if (d < 0.01) {
          ex = 0.01 * (i - j);
          ey = 0.01;
          d = Math.hypot(ex, ey);
        }
        const gap = a.r + b.r + 36;
        const f = (k * k) / d + (d < gap ? (gap - d) * 3 : 0);
        dx[i] += (ex / d) * f;
        dy[i] += (ey / d) * f;
        dx[j] -= (ex / d) * f;
        dy[j] -= (ey / d) * f;
      }
    for (const e of edges) {
      const i = index.get(e.from);
      const j = index.get(e.to);
      if (i === undefined || j === undefined || i === j) continue;
      const ex = views[i].x - views[j].x;
      const ey = views[i].y - views[j].y;
      const d = Math.hypot(ex, ey) || 0.01;
      const f = ((d * d) / k) * (0.35 + e.framesPct / maxPct);
      dx[i] -= (ex / d) * f;
      dy[i] -= (ey / d) * f;
      dx[j] += (ex / d) * f;
      dy[j] += (ey / d) * f;
    }
    views.forEach((v, i) => {
      dx[i] += (W / 2 - v.x) * 0.12;
      dy[i] += (H / 2 - v.y) * 0.12;
      const len = Math.hypot(dx[i], dy[i]) || 1;
      const step = Math.min(len, temp);
      v.x += (dx[i] / len) * step;
      v.y += (dy[i] / len) * step;
    });
  }
  // Fit the result into the frame instead of pinning nodes against its edges.
  const minX = Math.min(...views.map((v) => v.x - v.r));
  const maxX = Math.max(...views.map((v) => v.x + v.r));
  const minY = Math.min(...views.map((v) => v.y - v.r));
  const maxY = Math.max(...views.map((v) => v.y + v.r));
  const scale = Math.min((W - 160) / Math.max(maxX - minX, 1), (H - 110) / Math.max(maxY - minY, 1), 1.5);
  const cx = (minX + maxX) / 2;
  const cy = (minY + maxY) / 2;
  for (const v of views) {
    v.x = W / 2 + (v.x - cx) * scale;
    v.y = H / 2 + (v.y - cy) * scale;
  }
  separate(views);
}

/** Push overlapping nodes apart (room for the percentage label below each), staying inside the frame. */
function separate(views: NodeView[]): void {
  for (let pass = 0; pass < 300; pass++) {
    let moved = false;
    for (let i = 0; i < views.length; i++)
      for (let j = i + 1; j < views.length; j++) {
        const a = views[i];
        const b = views[j];
        const ex = b.x - a.x;
        const ey = b.y - a.y;
        const d = Math.hypot(ex, ey) || 0.01;
        const gap = a.r + b.r + 30;
        if (d >= gap) continue;
        const push = (gap - d) / 2;
        a.x -= (ex / d) * push;
        a.y -= (ey / d) * push;
        b.x += (ex / d) * push;
        b.y += (ey / d) * push;
        moved = true;
      }
    for (const v of views) {
      v.x = Math.min(W - v.r - 50, Math.max(v.r + 50, v.x));
      v.y = Math.min(H - v.r - 40, Math.max(v.r + 30, v.y));
    }
    if (!moved) break;
  }
}

function circleLayout(views: NodeView[]): void {
  const R = views.length === 1 ? 0 : Math.min(W, H) / 2 - 90;
  views.forEach((v, i) => {
    const angle = -Math.PI / 2 + (2 * Math.PI * i) / views.length;
    v.x = W / 2 + R * Math.cos(angle);
    v.y = H / 2 + R * Math.sin(angle);
  });
}

/** Transition graph between binding modes, recomputed from the per-frame mode sequence. */
export class MarkovPanel {
  readonly el: HTMLElement;
  private graph = h('div', {class: 'oa-mk-wrap'});
  private side = h('div', {class: 'oa-legend'});
  private ethogram = h('canvas', {class: 'oa-ethogram'});
  private disposers: (() => void)[] = [];
  private resize = new ResizeObserver(() => this.drawEthogram());
  private layout: Layout;
  private svg: SVGSVGElement | null = null;
  private scene: SVGGElement | null = null;
  private nodes = new Map<number, NodeView>();
  private edges: EdgeView[] = [];
  private focus: number | null = null;
  private view = {x: 0, y: 0, k: 1};
  private positions = new Map<number, {x: number; y: number}>();
  private structures = stored(STRUCTURES_KEY) === 'true';
  private ligand: LigandDepiction | null = null;
  private tips = new Map<number, HTMLElement>();

  constructor(private ctx: ResultsContext) {
    this.layout = stored(LAYOUT_KEY) === 'circle' ? 'circle' : 'force';
    const out = h('output', null, `${ctx.minTransition}%`);
    const range = h('input', {type: 'range', min: 0, max: 20, step: 0.5, value: ctx.minTransition});
    let timer = 0;
    range.addEventListener('input', () => {
      out.textContent = `${range.value}%`;
      clearTimeout(timer);
      timer = window.setTimeout(() => ctx.setMinTransition(Number(range.value)), 100);
    });
    const layoutSeg = h('div', {class: 'oa-seg'});
    for (const [id, label] of [['force', 'Force'], ['circle', 'Circle']] as [Layout, string][])
      layoutSeg.append(h('button', {class: id === this.layout ? 'active' : '', onclick: () => {
        this.layout = id;
        store(LAYOUT_KEY, id);
        layoutSeg.querySelectorAll('button').forEach((b) => b.classList.toggle('active', b.textContent === label));
        this.positions.clear();
        this.renderGraph();
      }}, label));
    const canDraw = Boolean(ctx.summary.ligand) && !ctx.summary.peptide;
    const structuresBtn = h('button', {class: `oa-btn oa-btn-sm ${this.structures ? 'oa-btn-primary' : ''}`,
      title: 'Draw every binding mode as the ligand with its interactions',
      onclick: () => {
        this.structures = !this.structures;
        structuresBtn.classList.toggle('oa-btn-primary', this.structures);
        store(STRUCTURES_KEY, String(this.structures));
        this.positions.clear();
        this.withLigand(() => this.renderGraph());
      }}, 'Structures');
    const figures = ctx.data.figures.markov;
    this.el = h('div', null,
      h('div', {class: 'oa-panel-head'},
        h('div', null, h('h3', null, 'Markov chain of binding modes'),
          h('div', {class: 'oa-small oa-muted'},
            'A transition is drawn when it makes up at least the minimal share of frames (openmmdl analysis -m); its reverse direction is drawn too. Drag nodes, scroll to zoom, drag the background to pan.')),
        h('div', {class: 'oa-panel-tools'},
          h('label', {class: 'oa-range', style: {width: '240px'}},
            h('span', {class: 'oa-small oa-muted', style: {whiteSpace: 'nowrap'}}, 'Min transition'), range, out),
          layoutSeg,
          canDraw ? structuresBtn : null,
          h('button', {class: 'oa-btn oa-btn-sm', title: 'Recompute the layout and reset zoom',
            onclick: () => { this.positions.clear(); this.view = {x: 0, y: 0, k: 1}; this.renderGraph(); }}, icon('reset'), 'Reset'),
          svgExportButtons(() => this.svg, () => `${ctx.summary.name}_${ctx.schemaId}_markov_chain`))),
      h('div', {style: {marginBottom: '10px'}}, thresholdControl(ctx)),
      h('div', {class: 'oa-markov'},
        h('div', null, this.graph,
          h('div', {class: 'oa-small oa-muted', style: {marginTop: '12px'}}, 'Binding mode of every analysed frame (10 most frequent coloured)'),
          this.ethogram, frameAxis(ctx.frames)),
        h('div', null, this.side,
          figures.length ? h('div', {class: 'oa-small', style: {marginTop: '14px'}},
            h('div', {class: 'oa-muted'}, 'OpenMMDL figures from the run (large images):'),
            ...figures.map((f) => h('div', null, h('a', {href: api.fileUrl(ctx.root, f.path), target: '_blank'},
              `markov_chain_plot_${f.threshold}`)))) : null)));

    this.ethogram.addEventListener('mousemove', (e) => {
      const p = positionAt(this.ethogram, e, ctx.nFrames);
      const mode = ctx.modeAt(ctx.frames[p]);
      showTooltip(e, h('div', null, h('b', null, `Frame ${ctx.frames[p]}`), mode ? ` · ${mode.label} (${pct(mode.percent)})` : ''));
    });
    this.ethogram.addEventListener('mouseleave', hideTooltip);
    this.ethogram.addEventListener('click', (e) => ctx.setFrame(ctx.frames[positionAt(this.ethogram, e, ctx.nFrames)]));
    this.resize.observe(this.ethogram);
    this.disposers.push(
      ctx.on('markov', () => this.render()),
      ctx.on('mode', () => { this.applyFocus(); this.renderSide(); }),
      ctx.on('frame', () => this.drawEthogram()),
    );
    this.render();
    if (canDraw) this.withLigand(() => this.structures && this.renderGraph());
  }

  /** Load the ligand depiction once (for structure nodes and tooltips), then run then(). */
  private withLigand(then: () => void): void {
    const name = this.ctx.summary.ligand;
    if (this.ligand || !name || this.ctx.summary.peptide) {
      then();
      return;
    }
    depiction(this.ctx.root, name).then((data) => {
      this.ligand = data;
      then();
    }).catch(() => then());
  }

  private modeFeatures(mode: number): FeatureView[] {
    const m = this.modeByIndex(mode);
    return (m?.features ?? []).map((id) => this.ctx.byId.get(id)).filter((f): f is FeatureView => Boolean(f));
  }

  private modeByIndex(index: number): BindingMode | undefined {
    return this.ctx.modes.modes.find((m) => m.index === index);
  }

  private render(): void {
    this.renderGraph();
    this.renderSide();
    this.drawEthogram();
  }

  private drawEthogram(): void {
    const colors = new Map<number, string>();
    for (const m of this.ctx.modes.modes) colors.set(m.index, modeColor(m));
    const cursor = this.ctx.position(this.ctx.frame);
    drawCategories(this.ethogram, this.ctx.modes.labels, (label) => colors.get(label) ?? MODE_OTHER, cursor >= 0 ? cursor : null);
  }

  private renderGraph(): void {
    clear(this.graph);
    this.tips.clear();
    this.nodes.clear();
    this.edges = [];
    const {nodes: allNodes, edges: allEdges, hiddenModes} = this.ctx.markov;
    if (!allNodes.length) {
      this.graph.append(h('div', {class: 'oa-empty'}, 'No transition reaches the minimal share of frames. Lower the minimum.'));
      return;
    }
    const nodes = allNodes.slice(0, MAX_NODES);
    const kept = new Set(nodes.map((n) => n.mode));
    const edges = allEdges.filter((e) => kept.has(e.from) && kept.has(e.to));
    const maxOcc = Math.max(...nodes.map((n) => n.occurrences));
    const drawn = this.structures && this.ligand !== null;
    const views: NodeView[] = nodes.map((node) => {
      const share = Math.sqrt(node.occurrences / maxOcc);
      // Cards keep the depiction's aspect ratio; r (their half width) drives layout spacing.
      const hh = 25 + 15 * share;
      const card = drawn ? {hw: hh * 1.36, hh} : null;
      return {node, mode: this.modeByIndex(node.mode), x: 0, y: 0, card,
        r: card ? card.hw : 15 + 23 * share, g: s('g')};
    });
    if (this.layout === 'force') forceLayout(views, edges);
    else circleLayout(views);
    for (const v of views) {
      const saved = this.positions.get(v.node.mode);
      if (saved) Object.assign(v, saved);
    }

    const svg = s('svg', {class: 'oa-markov-svg oa-mk', viewBox: `0 0 ${W} ${H}`});
    const defs = this.defs();
    svg.append(defs);
    const scene = s('g');
    const background = s('rect', {x: -5 * W, y: -5 * H, width: 11 * W, height: 11 * H, fill: 'transparent'});
    scene.append(background);
    const edgeLayer = s('g');
    const nodeLayer = s('g');
    scene.append(edgeLayer, nodeLayer);
    svg.append(scene);
    this.svg = svg;
    this.scene = scene;
    for (const v of views) this.nodes.set(v.node.mode, v);

    const maxPct = Math.max(1e-9, ...edges.map((e) => e.framesPct));
    const pairs = new Set(edges.map((e) => `${e.from}>${e.to}`));
    for (const edge of edges) {
      const weight = edge.framesPct / maxPct;
      const from = this.nodes.get(edge.from)!;
      const color = mix(MARKOV_COLORS[from.node.color].color, '#0f172a', 0.35);
      const width = 1.2 + 4.5 * weight;
      const hit = s('path', {fill: 'none', stroke: 'transparent', 'stroke-width': Math.max(14, width + 8)});
      const vis = s('path', {class: 'vis', fill: 'none', stroke: color, 'stroke-width': width, 'stroke-linecap': 'round',
        'marker-end': `url(#oa-mk-arrow-${from.node.color})`, opacity: 0.35 + 0.55 * weight});
      const label = s('text', {class: 'edge-label', 'text-anchor': 'middle'}, pct(edge.probability, 0));
      const g = s('g', {class: 'edge'}, hit, vis, label);
      g.dataset.from = String(edge.from);
      g.dataset.to = String(edge.to);
      g.dataset.bidir = String(pairs.has(`${edge.to}>${edge.from}`));
      g.addEventListener('mouseenter', () => this.highlightEdge(g, true));
      g.addEventListener('mousemove', (e) => showTooltip(e, this.edgeTip(edge)));
      g.addEventListener('mouseleave', () => { this.highlightEdge(g, false); hideTooltip(); });
      g.addEventListener('click', (e) => {
        e.stopPropagation();
        const target = this.modeByIndex(edge.to);
        if (target) this.ctx.selectMode(target);
      });
      edgeLayer.append(g);
      this.edges.push({edge, g, paths: [hit, vis], label});
    }

    views.forEach((v) => {
      const {node} = v;
      const g = v.g;
      g.setAttribute('class', 'node');
      const color = MARKOV_COLORS[node.color].color;
      const ring = 2 * Math.PI * (v.r + 6);
      const body = s('g', {class: 'body'},
        node.selfLoop ? s('path', {class: 'self', fill: 'none', stroke: mix(color, '#0f172a', 0.35), 'stroke-width': 2,
          'marker-end': `url(#oa-mk-arrow-${node.color})`}) : s('g'),
        ...(v.card ? this.structureCard(v, v.card, defs) : [
          s('circle', {class: 'halo', r: v.r + 11, fill: 'none'}),
          s('circle', {r: v.r + 6, fill: 'none', stroke: '#e2e8f0', 'stroke-width': 4}),
          s('circle', {r: v.r + 6, fill: 'none', stroke: mix(color, '#0f172a', 0.25), 'stroke-width': 4, 'stroke-linecap': 'round',
            'stroke-dasharray': `${(ring * node.selfProbability) / 100} ${ring}`, transform: 'rotate(-90)'}),
          s('circle', {class: 'disc', r: v.r, fill: `url(#oa-mk-grad-${node.color})`, stroke: mix(color, '#0f172a', 0.3),
            'stroke-width': 1.2, filter: 'url(#oa-mk-shadow)'}),
          s('text', {class: 'num', y: 5, 'text-anchor': 'middle'}, String(node.mode + 1))]),
        s('text', {class: 'cap', y: (v.card ? v.card.hh : v.r) + (v.card ? 18 : 24), 'text-anchor': 'middle'}, pct(node.framesPct)));
      g.append(body);
      g.addEventListener('mouseenter', () => this.hoverNode(node.mode));
      g.addEventListener('mousemove', (e) => showTooltip(e, this.nodeTip(node)));
      g.addEventListener('mouseleave', () => { this.hoverNode(null); hideTooltip(); });
      this.dragNode(g, v);
      nodeLayer.append(g);
    });
    this.positionAll();
    this.wirePanZoom(svg, background);
    this.applyView();
    this.graph.append(svg);
    this.applyFocus();
    const hidden = hiddenModes + Math.max(0, allNodes.length - nodes.length);
    if (hidden)
      this.graph.append(h('div', {class: 'oa-small oa-muted', style: {marginTop: '6px'}},
        `${hidden} rarer mode${hidden > 1 ? 's' : ''} not drawn (no transition above the minimum${allNodes.length > MAX_NODES ? `, or beyond the ${MAX_NODES} most frequent` : ''}).`));
  }

  /** A card showing the ligand with this mode's interactions, edged in the node colour, with a number badge
   * and a bar along the bottom for how often the mode persists from one frame to the next. */
  private structureCard(v: NodeView, card: {hw: number; hh: number}, defs: SVGDefsElement): SVGElement[] {
    const {node} = v;
    const {hw, hh} = card;
    const color = MARKOV_COLORS[node.color].color;
    const clipId = `oa-mk-clip-${node.mode}`;
    defs.append(s('clipPath', {id: clipId}, s('rect', {x: -hw, y: -hh, width: 2 * hw, height: 2 * hh, rx: 12})));
    const thumb = ligandThumbnail(this.ligand!, this.modeFeatures(node.mode));
    for (const [k, value] of Object.entries({x: -hw + 3, y: -hh + 3, width: 2 * hw - 6, height: 2 * hh - 10,
      preserveAspectRatio: 'xMidYMid meet'})) thumb.setAttribute(k, String(value));
    const barW = (2 * hw - 20) * (node.selfProbability / 100);
    return [
      s('rect', {class: 'halo', x: -hw - 6, y: -hh - 6, width: 2 * hw + 12, height: 2 * hh + 12, rx: 16, fill: 'none'}),
      s('rect', {class: 'disc', x: -hw, y: -hh, width: 2 * hw, height: 2 * hh, rx: 12, fill: '#ffffff', stroke: color,
        'stroke-width': 3, filter: 'url(#oa-mk-shadow)'}),
      s('g', {'clip-path': `url(#${clipId})`, 'pointer-events': 'none'}, thumb),
      s('rect', {x: -hw + 10, y: hh - 7, width: 2 * hw - 20, height: 3, rx: 1.5, fill: '#e2e8f0'}),
      s('rect', {x: -hw + 10, y: hh - 7, width: barW, height: 3, rx: 1.5, fill: mix(color, '#0f172a', 0.25)}),
      s('circle', {cx: hw - 4, cy: -hh + 4, r: 11, fill: `url(#oa-mk-grad-${node.color})`, stroke: mix(color, '#0f172a', 0.3)}),
      s('text', {class: 'num', x: hw - 4, y: -hh + 8, 'text-anchor': 'middle', style: 'font-size: 11px'}, String(node.mode + 1)),
    ];
  }

  /** Distance from a node's centre to its outline (plus margin) in direction (ux, uy). */
  private rim(v: NodeView, ux: number, uy: number): number {
    if (!v.card) return v.r + 6;
    return Math.min(v.card.hw / Math.max(Math.abs(ux), 1e-6), v.card.hh / Math.max(Math.abs(uy), 1e-6)) + 2;
  }

  private defs(): SVGDefsElement {
    const defs = s('defs');
    defs.append(s('filter', {id: 'oa-mk-shadow', x: '-30%', y: '-30%', width: '160%', height: '160%'},
      s('feDropShadow', {dx: 0, dy: 2, stdDeviation: 2.5, 'flood-color': '#0f172a', 'flood-opacity': 0.22})));
    for (const [key, info] of Object.entries(MARKOV_COLORS) as [ThirdColor, {color: string}][]) {
      defs.append(s('radialGradient', {id: `oa-mk-grad-${key}`, cx: '35%', cy: '30%', r: '75%'},
        s('stop', {offset: '0%', 'stop-color': mix(info.color, '#ffffff', 0.55)}),
        s('stop', {offset: '100%', 'stop-color': info.color})));
      defs.append(s('marker', {id: `oa-mk-arrow-${key}`, viewBox: '0 0 10 10', refX: 8, refY: 5, markerUnits: 'userSpaceOnUse',
        markerWidth: 11, markerHeight: 11, orient: 'auto-start-reverse'},
      s('path', {d: 'M0,1 L9,5 L0,9 L2.5,5 z', fill: mix(info.color, '#0f172a', 0.35)})));
    }
    return defs;
  }

  private positionAll(): void {
    for (const v of this.nodes.values()) {
      v.g.setAttribute('transform', `translate(${v.x.toFixed(1)},${v.y.toFixed(1)})`);
      const self = v.g.querySelector<SVGPathElement>('path.self');
      if (self) self.setAttribute('d', this.loopPath(v));
    }
    for (const view of this.edges) this.placeEdge(view);
  }

  /** Self-loop drawn on the outside of the node, away from the graph centre. */
  private loopPath(v: NodeView): string {
    let angle = Math.atan2(v.y - H / 2, v.x - W / 2);
    if (this.nodes.size === 1) angle = -Math.PI / 2;
    const r = this.rim(v, Math.cos(angle), Math.sin(angle));
    const a0 = angle - 0.45;
    const a1 = angle + 0.45;
    const p0 = [Math.cos(a0) * r, Math.sin(a0) * r];
    const p1 = [Math.cos(a1) * (r + 3), Math.sin(a1) * (r + 3)];
    const reach = r + 34 + v.r * 0.3;
    const c0 = [Math.cos(angle - 0.55) * reach, Math.sin(angle - 0.55) * reach];
    const c1 = [Math.cos(angle + 0.55) * reach, Math.sin(angle + 0.55) * reach];
    return `M${p0[0]},${p0[1]} C${c0[0]},${c0[1]} ${c1[0]},${c1[1]} ${p1[0]},${p1[1]}`;
  }

  private placeEdge(view: EdgeView): void {
    const a = this.nodes.get(view.edge.from)!;
    const b = this.nodes.get(view.edge.to)!;
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    const len = Math.hypot(dx, dy) || 1;
    const bend = view.g.dataset.bidir === 'true' ? 0.2 : 0.08;
    const mx = (a.x + b.x) / 2 + (-dy / len) * len * bend;
    const my = (a.y + b.y) / 2 + (dx / len) * len * bend;
    const unit = (p: NodeView) => {
      const l = Math.hypot(mx - p.x, my - p.y) || 1;
      return [(mx - p.x) / l, (my - p.y) / l] as const;
    };
    const start = towards(a, mx, my, this.rim(a, ...unit(a)) + 1);
    const end = towards(b, mx, my, this.rim(b, ...unit(b)) + 3);
    const d = `M${start.x.toFixed(1)},${start.y.toFixed(1)} Q${mx.toFixed(1)},${my.toFixed(1)} ${end.x.toFixed(1)},${end.y.toFixed(1)}`;
    view.paths.forEach((p) => p.setAttribute('d', d));
    view.label.setAttribute('x', String(0.25 * start.x + 0.5 * mx + 0.25 * end.x));
    view.label.setAttribute('y', String(0.25 * start.y + 0.5 * my + 0.25 * end.y - 4));
  }

  private toSvg(event: MouseEvent): {x: number; y: number} {
    const rect = this.svg!.getBoundingClientRect();
    const sx = ((event.clientX - rect.left) / rect.width) * W;
    const sy = ((event.clientY - rect.top) / rect.height) * H;
    return {x: (sx - this.view.x) / this.view.k, y: (sy - this.view.y) / this.view.k};
  }

  private dragNode(g: SVGGElement, v: NodeView): void {
    g.addEventListener('pointerdown', (down) => {
      down.stopPropagation();
      g.setPointerCapture(down.pointerId);
      const start = this.toSvg(down);
      const origin = {x: v.x, y: v.y};
      let moved = false;
      const move = (e: PointerEvent) => {
        const p = this.toSvg(e);
        if (!moved && Math.hypot(p.x - start.x, p.y - start.y) < 3) return;
        moved = true;
        v.x = origin.x + p.x - start.x;
        v.y = origin.y + p.y - start.y;
        this.positions.set(v.node.mode, {x: v.x, y: v.y});
        this.positionAll();
        hideTooltip();
      };
      const up = () => {
        g.removeEventListener('pointermove', move);
        g.removeEventListener('pointerup', up);
        if (!moved && v.mode) this.ctx.selectMode(v.mode);
      };
      g.addEventListener('pointermove', move);
      g.addEventListener('pointerup', up);
    });
  }

  private wirePanZoom(svg: SVGSVGElement, background: SVGRectElement): void {
    svg.addEventListener('wheel', (e) => {
      e.preventDefault();
      const rect = svg.getBoundingClientRect();
      const sx = ((e.clientX - rect.left) / rect.width) * W;
      const sy = ((e.clientY - rect.top) / rect.height) * H;
      const k = Math.min(4, Math.max(0.4, this.view.k * Math.exp(-e.deltaY * 0.0015)));
      this.view.x = sx - ((sx - this.view.x) * k) / this.view.k;
      this.view.y = sy - ((sy - this.view.y) * k) / this.view.k;
      this.view.k = k;
      this.applyView();
    }, {passive: false});
    background.addEventListener('pointerdown', (down) => {
      background.setPointerCapture(down.pointerId);
      const rect = svg.getBoundingClientRect();
      const scale = W / rect.width;
      const origin = {...this.view};
      let moved = false;
      const move = (e: PointerEvent) => {
        const ddx = (e.clientX - down.clientX) * scale;
        const ddy = (e.clientY - down.clientY) * scale;
        if (Math.hypot(ddx, ddy) > 3) moved = true;
        this.view.x = origin.x + ddx;
        this.view.y = origin.y + ddy;
        this.applyView();
      };
      const up = () => {
        background.removeEventListener('pointermove', move);
        background.removeEventListener('pointerup', up);
        if (!moved) {
          this.focus = null;
          this.applyFocus();
        }
      };
      background.addEventListener('pointermove', move);
      background.addEventListener('pointerup', up);
      svg.classList.add('panning');
      background.addEventListener('pointerup', () => svg.classList.remove('panning'), {once: true});
    });
  }

  private applyView(): void {
    this.scene?.setAttribute('transform', `translate(${this.view.x},${this.view.y}) scale(${this.view.k})`);
  }

  private hoverNode(mode: number | null): void {
    this.focus = mode;
    this.applyFocus();
  }

  /** Light up a node and its transitions; everything else fades. Falls back to the selected mode. */
  private applyFocus(): void {
    if (!this.svg) return;
    const center = this.focus ?? this.ctx.selectedMode?.index ?? null;
    const shown = center !== null && this.nodes.has(center);
    this.svg.classList.toggle('focus', shown && this.focus !== null);
    this.svg.classList.toggle('select', shown && this.focus === null);
    const near = new Set<number>(shown ? [center!] : []);
    for (const view of this.edges) {
      const on = shown && (view.edge.from === center || view.edge.to === center);
      view.g.classList.toggle('hl', on);
      view.g.classList.toggle('out', on && view.edge.from === center);
      if (on) {
        near.add(view.edge.from);
        near.add(view.edge.to);
      }
    }
    for (const [mode, v] of this.nodes) {
      v.g.classList.toggle('hl', near.has(mode));
      v.g.classList.toggle('selected', mode === this.ctx.selectedMode?.index);
    }
  }

  private highlightEdge(g: SVGGElement, on: boolean): void {
    g.classList.toggle('hover', on);
  }

  /** Built once per node and render: the depiction is too heavy to rebuild on every mouse move. */
  private nodeTip(node: MarkovNode): HTMLElement {
    const cached = this.tips.get(node.mode);
    if (cached) return cached;
    const mode = this.modeByIndex(node.mode);
    const tip = h('div', null, h('b', null, mode?.label ?? ''),
      this.ligand ? h('div', {class: 'oa-mk-tip-structure'}, ligandThumbnail(this.ligand, this.modeFeatures(node.mode))) : null,
      h('div', null, `Occurs in ${pct(node.framesPct)} of frames (${node.occurrences})`),
      h('div', null, `Stays in this mode: ${pct(node.selfProbability)} of its frames`),
      h('div', null, `Mean dwell: ${dwell(this.ctx.modes.labels, node.mode).toFixed(1)} frames`),
      h('div', {style: {opacity: '0.75'}}, MARKOV_COLORS[node.color].label),
      h('div', {style: {opacity: '0.75'}}, 'Click to select · drag to move'));
    this.tips.set(node.mode, tip);
    return tip;
  }

  private edgeTip(edge: MarkovEdge): HTMLElement {
    const from = this.modeByIndex(edge.from)?.label;
    const to = this.modeByIndex(edge.to)?.label;
    const back = this.ctx.markov.edges.find((e) => e.from === edge.to && e.to === edge.from);
    return h('div', null, h('b', null, `${from} → ${to}`),
      h('div', null, `${pct(edge.probability)} of the frames leaving ${from}`),
      h('div', null, `${edge.count} transition${edge.count === 1 ? '' : 's'} · ${pct(edge.framesPct, 2)} of all frames`),
      back ? h('div', {style: {opacity: '0.75'}}, `Back: ${pct(back.probability)} (${back.count})`) : null,
      h('div', {style: {opacity: '0.75'}}, `Click to go to ${to}`));
  }

  private renderSide(): void {
    clear(this.side);
    this.side.append(h('div', {class: 'oa-h4'}, 'Node colour'));
    for (const info of Object.values(MARKOV_COLORS))
      this.side.append(h('div', {class: 'oa-legend-row'}, h('span', {class: 'sw', style: {background: info.color}}), info.label));
    this.side.append(h('div', {class: 'oa-small oa-muted'},
      'Node size = occurrence. Ring = how often the mode persists from one frame to the next. Edge width = share of frames.'));
    const mode = this.ctx.selectedMode;
    if (!mode) {
      this.side.append(h('div', {class: 'oa-small oa-muted', style: {marginTop: '12px'}}, 'Click a node to inspect a mode.'));
      return;
    }
    const node = this.ctx.markov.nodes.find((n) => n.mode === mode.index);
    const list = (title: string, items: {other: number; share: number}[]) => {
      if (!items.length) return null;
      const box = h('div', {class: 'oa-mk-list'}, h('div', {class: 'oa-small oa-muted'}, title));
      for (const item of items.slice(0, 8)) {
        const other = this.modeByIndex(item.other);
        box.append(h('button', {class: 'oa-mk-item', type: 'button', onclick: () => other && this.ctx.selectMode(other)},
          h('span', {class: 'sw', style: {background: modeColor(other)}}),
          h('span', {class: 'name'}, other?.label.replace('Binding_Mode_', 'Mode ') ?? ''),
          h('span', {class: 'bar'}, h('span', {style: {width: `${Math.min(100, item.share)}%`}})),
          h('span', {class: 'val'}, pct(item.share, 0))));
      }
      return box;
    };
    const edges = this.ctx.markov.edges;
    const incoming = edges.filter((e) => e.to === mode.index);
    const arrivals = incoming.reduce((sum, e) => sum + e.count, 0) || 1;
    this.side.append(
      h('div', {class: 'oa-mk-card'},
        h('div', {class: 'oa-h4', style: {margin: '0 0 4px'}}, mode.label),
        h('div', {class: 'oa-small'}, `${pct(mode.percent)} of frames · stays ${pct(node?.selfProbability ?? 0)} · mean dwell ${dwell(this.ctx.modes.labels, mode.index).toFixed(1)} frames`),
        list('Goes to (probability of leaving to)', edges.filter((e) => e.from === mode.index)
          .sort((a, b) => b.probability - a.probability).map((e) => ({other: e.to, share: e.probability}))),
        list('Comes from (share of arrivals)', incoming.sort((a, b) => b.count - a.count)
          .map((e) => ({other: e.from, share: (100 * e.count) / arrivals})))));
  }

  dispose(): void {
    this.resize.disconnect();
    this.disposers.forEach((d) => d());
  }
}

function towards(p: {x: number; y: number}, tx: number, ty: number, dist: number): {x: number; y: number} {
  const dx = tx - p.x;
  const dy = ty - p.y;
  const len = Math.hypot(dx, dy) || 1;
  return {x: p.x + (dx / len) * dist, y: p.y + (dy / len) * dist};
}
