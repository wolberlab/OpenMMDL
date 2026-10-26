import {TYPES, typeInfo} from '../colors';
import {clear, h, hideTooltip, icon, pct, showTooltip} from '../dom';
import {canvasToPng, download, downloadCsv, exportButtons, safeName} from '../export';
import {drawBarcode, frameAxis, positionAt} from './barcode';
import {FeatureView, ResultsContext} from './context';

type SortKey = 'occurrence' | 'residue' | 'type' | 'first';

interface Row {
  el: HTMLElement;
  canvas: HTMLCanvasElement;
  feature: FeatureView;
}

function firstFrame(f: FeatureView): number {
  return f.runs.length ? f.runs[0] : Number.MAX_SAFE_INTEGER;
}

function residueKey(f: FeatureView): string {
  const r = f.residues[0];
  return r ? `${r.chain}:${String(r.resnr).padStart(6, '0')}` : `~${f.partner}`;
}

/** Barcode timeline of every interaction in df_all.csv. */
export class InteractionsPanel {
  readonly el: HTMLElement;
  private list = h('div', {class: 'oa-int-list'});
  private countLabel = h('span', {class: 'oa-muted oa-small'});
  private rows: Row[] = [];
  private enabled = new Set(TYPES.map((t) => t.key));
  private search = '';
  private sort: SortKey = 'occurrence';
  private minPct = 0;
  private onlyPassing = false;
  private disposers: (() => void)[] = [];
  private resize = new ResizeObserver(() => this.redraw());

  constructor(private ctx: ResultsContext, private filterType: ((type: string) => boolean) | null = null) {
    const present = new Map<string, number>();
    for (const f of ctx.features) present.set(f.type, (present.get(f.type) ?? 0) + 1);
    const chips = h('div', {class: 'oa-chips'});
    for (const type of TYPES) {
      const n = present.get(type.key);
      if (!n || (filterType && !filterType(type.key))) continue;
      const chip = h('button', {class: 'oa-chip', title: 'Click to toggle, double-click to show only this type'},
        h('span', {class: 'sw', style: {background: type.color}}), type.label, h('span', {class: 'n'}, String(n)));
      chip.addEventListener('click', () => {
        if (this.enabled.has(type.key)) this.enabled.delete(type.key);
        else this.enabled.add(type.key);
        chip.classList.toggle('off', !this.enabled.has(type.key));
        this.render();
      });
      chip.addEventListener('dblclick', () => {
        this.enabled = new Set([type.key]);
        chips.querySelectorAll('.oa-chip').forEach((c) => c.classList.toggle('off', c !== chip));
        this.render();
      });
      chips.append(chip);
    }

    const searchInput = h('input', {class: 'oa-input', placeholder: 'Filter by residue, atom or type',
      style: {width: '220px'}, oninput: () => { this.search = searchInput.value.trim().toLowerCase(); this.render(); }});
    const sortSelect = h('select', {class: 'oa-select', style: {width: 'auto'},
      onchange: () => { this.sort = sortSelect.value as SortKey; this.render(); }},
    h('option', {value: 'occurrence'}, 'Sort: occurrence'), h('option', {value: 'residue'}, 'Sort: residue'),
    h('option', {value: 'type'}, 'Sort: type'), h('option', {value: 'first'}, 'Sort: first frame'));
    const minOut = h('output', null, '0%');
    const minRange = h('input', {type: 'range', min: 0, max: 100, step: 1, value: 0,
      oninput: () => { this.minPct = Number(minRange.value); minOut.textContent = `${this.minPct}%`; this.render(); }});
    const passingBox = h('input', {type: 'checkbox',
      onchange: () => { this.onlyPassing = passingBox.checked; this.render(); }});

    this.el = h('div', null,
      h('div', {class: 'oa-panel-head'},
        h('div', null, h('h3', null, filterType ? 'Water bridges' : 'Interactions over the trajectory'), this.countLabel),
        h('div', {class: 'oa-panel-tools'},
          h('div', {class: 'oa-search', style: {display: 'flex', gap: '6px', alignItems: 'center'}}, icon('search'), searchInput),
          sortSelect,
          h('label', {class: 'oa-range', title: 'Hide interactions present in fewer frames', style: {width: '200px'}},
            h('span', {class: 'oa-small oa-muted'}, 'Min'), minRange, minOut),
          h('label', {class: 'oa-check oa-small', title: 'Only interactions that pass the binding-mode threshold'},
            passingBox, 'Above mode threshold'),
          exportButtons({PNG: () => this.exportPng(), CSV: () => this.exportCsv()}))),
      chips.childElementCount > 1 ? h('div', {style: {marginBottom: '10px'}}, chips) : null,
      h('div', {class: 'oa-int-axis'}, h('span', null, 'Interaction'), frameAxis(ctx.frames), h('span', {style: {textAlign: 'right'}}, 'Frames')),
      this.list);

    this.resize.observe(this.list);
    this.disposers.push(
      ctx.on('frame', () => this.redraw()),
      ctx.on('features', () => this.markSelected()),
      ctx.on('modes', () => { if (this.onlyPassing) this.render(); }),
    );
    this.render();
  }

  private visible(): FeatureView[] {
    const passing = this.ctx.passing();
    const minCount = (this.minPct / 100) * this.ctx.nFrames;
    const items = this.ctx.features.filter((f) =>
      this.enabled.has(f.type) && (!this.filterType || this.filterType(f.type)) &&
      f.count >= minCount && (!this.onlyPassing || passing.has(f.id)) &&
      (!this.search || `${f.label} ${f.sublabel} ${f.id}`.toLowerCase().includes(this.search)));
    const typeOrder = new Map(TYPES.map((t, i) => [t.key, i]));
    const compare: Record<SortKey, (a: FeatureView, b: FeatureView) => number> = {
      occurrence: (a, b) => b.count - a.count,
      residue: (a, b) => residueKey(a).localeCompare(residueKey(b)) || b.count - a.count,
      type: (a, b) => (typeOrder.get(a.type) ?? 99) - (typeOrder.get(b.type) ?? 99) || b.count - a.count,
      first: (a, b) => firstFrame(a) - firstFrame(b) || b.count - a.count,
    };
    return items.sort(compare[this.sort]);
  }

  render(): void {
    clear(this.list);
    this.rows = [];
    const items = this.visible();
    this.countLabel.textContent = `${items.length} of ${this.ctx.features.filter((f) => !this.filterType || this.filterType(f.type)).length} shown · bars mark the frames where an interaction is present`;
    if (!items.length) {
      this.list.append(h('div', {class: 'oa-empty'}, 'No interactions match the filters.'));
      return;
    }
    for (const feature of items) {
      const type = typeInfo(feature.type);
      const canvas = h('canvas', {class: 'oa-barcode'});
      const el = h('div', {class: 'oa-int-row', title: feature.id},
        h('div', {class: 'oa-int-label'}, h('b', null, feature.label), h('span', null, feature.sublabel)),
        canvas,
        h('div', {class: 'oa-int-pct', style: {color: type.color}}, pct((100 * feature.count) / this.ctx.nFrames, 0)));
      el.addEventListener('click', (e) => {
        if (e.target === canvas) return;
        this.ctx.selectFeatures([feature.id], e.ctrlKey || e.metaKey || e.shiftKey);
      });
      el.addEventListener('mouseenter', () => this.ctx.hoverFeatures([feature.id]));
      el.addEventListener('mouseleave', () => this.ctx.hoverFeatures([]));
      canvas.addEventListener('mousemove', (e) => {
        const p = positionAt(canvas, e, this.ctx.nFrames);
        const frame = this.ctx.frames[p];
        const mode = this.ctx.modeAt(frame);
        showTooltip(e, h('div', null, h('b', null, `Frame ${frame}`), ` · ${feature.flags[p] ? 'present' : 'absent'}`,
          mode ? h('div', null, mode.label) : null, h('div', {class: 'oa-muted'}, 'Click to show this frame in 3D')));
      });
      canvas.addEventListener('mouseleave', hideTooltip);
      canvas.addEventListener('click', (e) => {
        this.ctx.selectFeatures([feature.id], e.ctrlKey || e.metaKey || e.shiftKey);
        this.ctx.setFrame(this.ctx.frames[positionAt(canvas, e, this.ctx.nFrames)]);
      });
      this.rows.push({el, canvas, feature});
      this.list.append(el);
    }
    this.markSelected();
    this.redraw();
  }

  private markSelected(): void {
    for (const row of this.rows) row.el.classList.toggle('selected', this.ctx.selected.has(row.feature.id));
  }

  redraw(): void {
    const cursor = this.ctx.position(this.ctx.frame);
    for (const row of this.rows)
      drawBarcode(row.canvas, row.feature.flags, typeInfo(row.feature.type).color, cursor >= 0 ? cursor : null);
  }

  private fileName(): string {
    return `${this.ctx.summary.name}_${this.filterType ? 'water_bridges' : 'interactions'}_${this.ctx.schemaId}`;
  }

  private exportCsv(): void {
    const frames = this.ctx.frames;
    downloadCsv([
      ['interaction', 'label', 'type', 'ligand_atoms', 'frames_present', 'occurrence_percent', 'present_in_frames'],
      ...this.rows.map(({feature: f}) => [f.id, f.label, typeInfo(f.type).label, f.ligandAtoms.join(' '), f.count,
        ((100 * f.count) / this.ctx.nFrames).toFixed(2), frames.filter((_, i) => f.flags[i]).join(' ')]),
    ], this.fileName());
  }

  /** The barcode chart of the rows currently shown, as one image. */
  private async exportPng(): Promise<void> {
    const rowH = 30;
    const labelW = 420;
    const barW = 1400;
    const top = 40;
    const width = labelW + barW + 120;
    const height = top + this.rows.length * rowH + 20;
    const n = this.ctx.nFrames;
    const blob = await canvasToPng(width, height, (g) => {
      g.font = '600 15px Inter, Segoe UI, Arial, sans-serif';
      g.fillStyle = '#64748b';
      g.textAlign = 'center';
      for (let t = 0; t <= 6; t++) {
        const p = Math.min(n - 1, Math.round((t / 6) * (n - 1)));
        g.fillText(String(this.ctx.frames[p]), labelW + ((p + 0.5) / n) * barW, 24);
      }
      this.rows.forEach(({feature: f}, r) => {
        const y = top + r * rowH;
        const color = typeInfo(f.type).color;
        g.textAlign = 'left';
        g.fillStyle = '#0f172a';
        g.font = '600 14px Inter, Segoe UI, Arial, sans-serif';
        g.fillText(f.label, 8, y + 13);
        g.fillStyle = '#64748b';
        g.font = '12px Inter, Segoe UI, Arial, sans-serif';
        g.fillText(f.sublabel.slice(0, 60), 8, y + 26);
        g.fillStyle = '#f1f5f9';
        g.fillRect(labelW, y + 4, barW, rowH - 8);
        g.fillStyle = color;
        for (let i = 0; i < n; i++) if (f.flags[i]) g.fillRect(labelW + (i / n) * barW, y + 4, Math.max(1, barW / n), rowH - 8);
        g.textAlign = 'right';
        g.font = '700 14px Inter, Segoe UI, Arial, sans-serif';
        g.fillText(`${Math.round((100 * f.count) / n)}%`, width - 10, y + 20);
      });
    });
    download(blob, `${safeName(this.fileName())}.png`);
  }

  dispose(): void {
    this.resize.disconnect();
    this.disposers.forEach((d) => d());
  }
}
