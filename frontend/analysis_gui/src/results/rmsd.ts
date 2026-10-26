import {api} from '../api';
import {clear, h, hideTooltip, s, showTooltip} from '../dom';
import {downloadCsv, exportButtons, svgExportButtons} from '../export';
import {ResultsContext} from './context';

const COLORS = ['#4f46e5', '#0ea5e9', '#f97316', '#10b981', '#e11d48'];
const W = 900;
const H = 300;
const M = {left: 48, right: 16, top: 12, bottom: 32};

export function seriesLabel(name: string): string {
  if (name === 'backbone') return 'Backbone';
  if (name === 'protein') return 'Protein';
  if (name === 'nucleicbackbone') return 'Nucleic backbone';
  if (name === 'nucleic') return 'Nucleic acid';
  if (name.startsWith('resname ')) return `Ligand (${name.slice(8)})`;
  if (name.startsWith('chainID ')) return `Peptide (chain ${name.slice(8)})`;
  return name;
}

function niceStep(range: number, ticks: number): number {
  const raw = range / ticks;
  const mag = 10 ** Math.floor(Math.log10(raw || 1));
  const norm = raw / mag;
  return (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
}

export function stats(values: (number | null)[]): {mean: number; median: number; max: number; last: number} {
  const v = values.filter((x): x is number => x !== null);
  const sorted = [...v].sort((a, b) => a - b);
  return {
    mean: v.reduce((a, b) => a + b, 0) / Math.max(v.length, 1),
    median: sorted[Math.floor(sorted.length / 2)] ?? 0,
    max: sorted[sorted.length - 1] ?? 0,
    last: v[v.length - 1] ?? 0,
  };
}

/** RMSD over time (RMSD/RMSD_over_time.csv) as an interactive line chart. */
export class RmsdPanel {
  readonly el: HTMLElement;
  private chart = h('div');
  private hidden = new Set<string>();
  private cursor: SVGLineElement | null = null;
  private xOf: ((frame: number) => number) | null = null;
  private disposers: (() => void)[] = [];

  constructor(private ctx: ResultsContext) {
    const rmsd = ctx.summary.rmsd;
    const matrix = ctx.summary.rmsdMatrix;
    const matrixFile = ctx.summary.files.find((f) => f.path === matrix);
    this.el = h('div', null,
      h('div', {class: 'oa-panel-head'},
        h('div', null, h('h3', null, 'RMSD over time'),
          h('div', {class: 'oa-small oa-muted'}, 'Fitted on the backbone of frame 0. Click the chart to show that frame in 3D.')),
        rmsd ? h('div', {class: 'oa-panel-tools'},
          svgExportButtons(() => this.chart.querySelector('svg'), () => `${ctx.summary.name}_rmsd`),
          exportButtons({CSV: () => downloadCsv([['frame', ...rmsd.series.map((sr) => seriesLabel(sr.name))],
            ...rmsd.frames.map((f, i) => [f, ...rmsd.series.map((sr) => sr.values[i])])], `${ctx.summary.name}_rmsd`)})) : null),
      rmsd ? this.chart : h('div', {class: 'oa-empty'}, 'No RMSD/RMSD_over_time.csv in this folder.'),
      rmsd ? this.statsTable() : null,
      matrix ? h('div', {style: {marginTop: '18px'}},
        h('div', {class: 'oa-h4'}, 'RMSD between frames'),
        matrixFile?.inline === false ? h('a', {href: api.fileUrl(ctx.root, matrix), target: '_blank'}, 'Open the RMSD matrix figure') :
          h('img', {class: 'oa-fig', src: api.fileUrl(ctx.root, matrix), alt: 'RMSD between frames', loading: 'lazy'})) : null);
    if (rmsd) this.draw();
    this.disposers.push(ctx.on('frame', () => this.moveCursor()));
  }

  private statsTable(): HTMLElement {
    const rows = this.ctx.summary.rmsd!.series.map((series) => {
      const st = stats(series.values);
      return h('tr', null, h('td', null, seriesLabel(series.name)),
        ...[st.mean, st.median, st.max, st.last].map((v) => h('td', null, `${v.toFixed(2)} Å`)));
    });
    return h('table', {class: 'oa-kv', style: {marginTop: '10px', maxWidth: '640px'}},
      h('tr', null, ...['', 'Mean', 'Median', 'Max', 'Last frame'].map((t) => h('td', {style: {fontWeight: '700'}}, t))),
      ...rows);
  }

  private draw(): void {
    clear(this.chart);
    const rmsd = this.ctx.summary.rmsd!;
    const series = rmsd.series.filter((sr) => !this.hidden.has(sr.name));
    const frames = rmsd.frames;
    const x0 = frames[0];
    const x1 = frames[frames.length - 1] || 1;
    const yMax = Math.max(0.5, ...series.flatMap((sr) => sr.values.filter((v): v is number => v !== null))) * 1.08;
    const xOf = (f: number) => M.left + ((f - x0) / Math.max(x1 - x0, 1)) * (W - M.left - M.right);
    const yOf = (v: number) => H - M.bottom - (v / yMax) * (H - M.top - M.bottom);
    this.xOf = xOf;

    const svg = s('svg', {class: 'oa-chart', viewBox: `0 0 ${W} ${H}`});
    const yStep = niceStep(yMax, 5);
    for (let v = 0; v <= yMax; v += yStep) {
      svg.append(s('line', {class: 'grid', x1: M.left, x2: W - M.right, y1: yOf(v), y2: yOf(v)}));
      svg.append(s('text', {x: M.left - 6, y: yOf(v) + 4, 'text-anchor': 'end'}, v.toFixed(yStep < 1 ? 1 : 0)));
    }
    const xStep = niceStep(x1 - x0, 8);
    for (let f = Math.ceil(x0 / xStep) * xStep; f <= x1; f += xStep)
      svg.append(s('text', {x: xOf(f), y: H - M.bottom + 16, 'text-anchor': 'middle'}, String(f)));
    svg.append(s('line', {class: 'axis', x1: M.left, x2: W - M.right, y1: H - M.bottom, y2: H - M.bottom}));
    svg.append(s('text', {x: (W + M.left) / 2, y: H - 2, 'text-anchor': 'middle'}, 'Frame'));
    svg.append(s('text', {x: 12, y: (H - M.bottom) / 2, transform: `rotate(-90 12 ${(H - M.bottom) / 2})`, 'text-anchor': 'middle'}, 'RMSD (Å)'));

    rmsd.series.forEach((sr, i) => {
      if (this.hidden.has(sr.name)) return;
      let d = '';
      sr.values.forEach((v, j) => {
        if (v === null) return;
        d += `${d ? 'L' : 'M'}${xOf(frames[j]).toFixed(1)},${yOf(v).toFixed(1)}`;
      });
      svg.append(s('path', {d, fill: 'none', stroke: COLORS[i % COLORS.length], 'stroke-width': 1.8}));
    });

    this.cursor = s('line', {y1: M.top, y2: H - M.bottom, stroke: '#0f172a', 'stroke-width': 1, 'stroke-dasharray': '3 3'});
    const hover = s('line', {y1: M.top, y2: H - M.bottom, stroke: '#94a3b8', 'stroke-width': 1, visibility: 'hidden'});
    svg.append(this.cursor, hover);
    const overlay = s('rect', {x: M.left, y: M.top, width: W - M.left - M.right, height: H - M.top - M.bottom,
      fill: 'transparent', style: 'cursor: crosshair'});
    const frameAt = (e: MouseEvent) => {
      const rect = svg.getBoundingClientRect();
      const x = ((e.clientX - rect.left) / rect.width) * W;
      const t = (x - M.left) / (W - M.left - M.right);
      const target = x0 + t * (x1 - x0);
      let best = 0;
      for (let j = 1; j < frames.length; j++) if (Math.abs(frames[j] - target) < Math.abs(frames[best] - target)) best = j;
      return best;
    };
    overlay.addEventListener('mousemove', (e) => {
      const j = frameAt(e);
      hover.setAttribute('x1', String(xOf(frames[j])));
      hover.setAttribute('x2', String(xOf(frames[j])));
      hover.setAttribute('visibility', 'visible');
      showTooltip(e, h('div', null, h('b', null, `Frame ${frames[j]}`),
        ...rmsd.series.filter((sr) => !this.hidden.has(sr.name)).map((sr) =>
          h('div', null, `${seriesLabel(sr.name)}: ${sr.values[j]?.toFixed(2) ?? '–'} Å`))));
    });
    overlay.addEventListener('mouseleave', () => { hover.setAttribute('visibility', 'hidden'); hideTooltip(); });
    overlay.addEventListener('click', (e) => this.ctx.setFrame(frames[frameAt(e)]));
    svg.append(overlay);

    const legend = h('div', {class: 'oa-chips', style: {marginBottom: '8px'}});
    rmsd.series.forEach((sr, i) => {
      const chip = h('button', {class: `oa-chip ${this.hidden.has(sr.name) ? 'off' : ''}`},
        h('span', {class: 'sw', style: {background: COLORS[i % COLORS.length]}}), seriesLabel(sr.name));
      chip.addEventListener('click', () => {
        if (this.hidden.has(sr.name)) this.hidden.delete(sr.name);
        else this.hidden.add(sr.name);
        this.draw();
      });
      legend.append(chip);
    });
    this.chart.append(legend, svg);
    this.moveCursor();
  }

  private moveCursor(): void {
    if (!this.cursor || !this.xOf) return;
    const x = this.xOf(this.ctx.frame);
    this.cursor.setAttribute('x1', String(x));
    this.cursor.setAttribute('x2', String(x));
  }

  dispose(): void {
    this.disposers.forEach((d) => d());
  }
}
