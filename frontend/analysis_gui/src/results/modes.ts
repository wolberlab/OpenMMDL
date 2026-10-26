import {TYPES, typeInfo} from '../colors';
import {clear, h, hideTooltip, pct, showTooltip} from '../dom';
import {canvasToPng, download, downloadCsv, exportButtons, safeName, standaloneSvg, svgToPng} from '../export';
import {BindingMode} from '../derive';
import {drawBarcode, frameAxis, positionAt} from './barcode';
import {FeatureView, ResultsContext} from './context';
import {Ligand2D} from './ligand2d';
import {LigandNetwork} from './network2d';

const TOP = 10;
const VIEW_KEY = 'openmmdl-analysis-gui:mode-2d';
type View2D = 'depiction' | 'network';

function modeFlags(mode: BindingMode, n: number): Uint8Array {
  const flags = new Uint8Array(n);
  for (const p of mode.positions) flags[p] = 1;
  return flags;
}

function storedView(): View2D {
  try {
    return localStorage.getItem(VIEW_KEY) === 'network' ? 'network' : 'depiction';
  } catch {
    return 'depiction';
  }
}

/** Threshold slider shared by the Binding modes and Markov panels. */
export function thresholdControl(ctx: ResultsContext): HTMLElement {
  const out = h('output', null, `${ctx.threshold}%`);
  const info = h('span', {class: 'oa-small oa-muted'});
  const range = h('input', {type: 'range', min: 0, max: 100, step: 1, value: ctx.threshold});
  const update = () => {
    info.textContent = `${ctx.modes.passing.length} interactions pass · ${ctx.modes.modes.length} modes`;
  };
  let timer = 0;
  range.addEventListener('input', () => {
    out.textContent = `${range.value}%`;
    clearTimeout(timer);
    timer = window.setTimeout(() => ctx.setThreshold(Number(range.value)), 120);
  });
  const reset = h('button', {class: 'oa-btn oa-btn-sm', title: `Back to the run's threshold (${ctx.summary.threshold}%)`,
    onclick: () => { range.value = String(ctx.summary.threshold); out.textContent = `${range.value}%`; ctx.setThreshold(ctx.summary.threshold); }},
  'Reset');
  ctx.on('modes', () => {
    range.value = String(ctx.threshold);
    out.textContent = `${ctx.threshold}%`;
    update();
  });
  update();
  return h('div', {class: 'oa-panel-tools'},
    h('label', {class: 'oa-range', style: {width: '260px'},
      title: 'An interaction defines binding modes when it is present in at least this share of frames (openmmdl analysis -b)'},
    h('span', {class: 'oa-small oa-muted', style: {whiteSpace: 'nowrap'}}, 'Mode threshold'), range, out),
    reset, info);
}

/** Binding modes recomputed in the browser, with interactions that can be switched on and off. */
export class ModesPanel {
  readonly el: HTMLElement;
  private strip = h('div', {class: 'oa-mode-strip-row'});
  private stripInfo = h('div', {class: 'oa-small oa-muted'});
  private detail = h('div', {class: 'oa-mode-detail'});
  private slot3d = h('div', {class: 'oa-mode-3d'});
  private views = h('div', {class: 'oa-mode-views'});
  private view2dHost = h('div');
  private chipBox = h('div', {class: 'oa-mode-chips'});
  private typeBox = h('div', {class: 'oa-chips'});
  private chipCount = h('span', {class: 'oa-small oa-muted'});
  private timeline = h('canvas', {class: 'oa-barcode', style: {height: '18px'}});
  private showAll = false;
  private view2d: View2D = storedView();
  private cards = new Map<number, {el: HTMLElement; canvas: HTMLCanvasElement; mode: BindingMode}>();
  private chips = new Map<string, HTMLElement>();
  private ligand: Ligand2D;
  private network: LigandNetwork;
  private emphasis: string[] = [];
  private disposers: (() => void)[] = [];
  private resize = new ResizeObserver(() => this.drawStrips());

  constructor(private ctx: ResultsContext) {
    const toggle = (ids: string[]) => ctx.toggleHidden(ids);
    this.ligand = new Ligand2D(ctx.root, ctx.summary.peptide ? null : ctx.summary.ligand, toggle, (ids) => this.hover(ids));
    this.network = new LigandNetwork(ctx, toggle);
    this.el = h('div', null,
      h('div', {class: 'oa-panel-head'},
        h('div', null, h('h3', null, 'Binding modes'),
          h('div', {class: 'oa-small oa-muted'},
            'Frames with the same set of interactions above the threshold form one mode, as in openmmdl analysis. Recomputed live from df_all.csv.')),
        thresholdControl(ctx)),
      this.stripInfo, this.strip, this.detail);
    this.resize.observe(this.strip);
    this.timeline.addEventListener('mousemove', (e) => {
      const p = positionAt(this.timeline, e, ctx.nFrames);
      const mode = ctx.modeAt(ctx.frames[p]);
      showTooltip(e, h('div', null, h('b', null, `Frame ${ctx.frames[p]}`), mode ? ` · ${mode.label}` : '',
        h('div', {class: 'oa-muted'}, 'Click to show this frame')));
    });
    this.timeline.addEventListener('mouseleave', hideTooltip);
    this.timeline.addEventListener('click', (e) => ctx.setFrame(ctx.frames[positionAt(this.timeline, e, ctx.nFrames)]));
    this.disposers.push(
      ctx.on('modes', () => this.render()),
      ctx.on('mode', () => { this.markSelected(); this.renderDetail(); }),
      ctx.on('hidden', () => this.refreshToggles()),
      ctx.on('frame', () => this.drawStrips()),
    );
    this.render();
  }

  /** Host the shared 3D view next to the 2D one (full-width layout); null gives the slot back. */
  setEmbedded3D(el: HTMLElement | null): void {
    clear(this.slot3d);
    if (el) this.slot3d.append(el);
    this.views.classList.toggle('with-3d', Boolean(el));
  }

  /** True when the run's own top-10 table (same threshold) groups frames differently. */
  private differsFromRun(): boolean {
    const top = this.ctx.data.topModes;
    if (!top?.length || this.ctx.threshold !== this.ctx.summary.threshold) return false;
    return top.some((t) => {
      const mine = this.ctx.modes.modes.find((m) => m.label === t.mode);
      return !mine || Math.abs(mine.percent - t.percent) > 1e-6;
    });
  }

  private current(): BindingMode | null {
    return this.ctx.selectedMode ?? this.ctx.modes.modes[0] ?? null;
  }

  render(): void {
    clear(this.strip);
    this.cards.clear();
    const modes = this.ctx.modes.modes;
    const shown = this.showAll ? modes : modes.slice(0, TOP);
    const covered = shown.reduce((sum, m) => sum + m.count, 0);
    this.stripInfo.replaceChildren(
      `${shown.length === modes.length ? 'All' : `Top ${shown.length} of`} ${modes.length} modes · ${pct((100 * covered) / this.ctx.nFrames)} of frames`,
      this.differsFromRun() ? h('div', {class: 'oa-note warn', style: {margin: '6px 0'}},
        'The run’s own top-10 table and figures group frames differently: openmmdl analysis leaves out one interaction column when it builds modes. The modes here use every interaction above the threshold.') : '');
    for (const mode of shown) {
      const canvas = h('canvas', {class: 'oa-mode-strip'});
      const el = h('button', {class: 'oa-mode-card', type: 'button', onclick: () => this.ctx.selectMode(mode)},
        h('div', {class: 'oa-mode-top'}, h('span', {class: 'oa-mode-name'}, `#${mode.rank + 1} · ${mode.label.replace('Binding_Mode_', 'Mode ')}`),
          h('span', {class: 'oa-mode-pct'}, pct(mode.percent))),
        h('div', {class: 'oa-mode-sub'}, `${mode.count} frames · ${mode.features.length} interactions`),
        canvas);
      this.cards.set(mode.index, {el, canvas, mode});
      this.strip.append(el);
    }
    if (modes.length > TOP)
      this.strip.append(h('button', {class: 'oa-btn oa-btn-sm', style: {alignSelf: 'center'},
        onclick: () => { this.showAll = !this.showAll; this.render(); }},
      this.showAll ? `Top ${TOP} only` : `All ${modes.length} modes`));
    this.markSelected();
    this.drawStrips();
    this.renderDetail();
  }

  private markSelected(): void {
    const current = this.current();
    for (const [index, card] of this.cards) {
      card.el.classList.toggle('selected', index === current?.index);
      if (index === current?.index) card.el.scrollIntoView({block: 'nearest', inline: 'nearest'});
    }
  }

  private drawStrips(): void {
    const cursor = this.ctx.position(this.ctx.frame);
    const c = cursor >= 0 ? cursor : null;
    for (const card of this.cards.values()) drawBarcode(card.canvas, modeFlags(card.mode, this.ctx.nFrames), '#4f46e5', c);
    const mode = this.current();
    if (mode) drawBarcode(this.timeline, modeFlags(mode, this.ctx.nFrames), '#4f46e5', c);
  }

  private modeFeatures(mode: BindingMode): FeatureView[] {
    const typeOrder = new Map(TYPES.map((t, i) => [t.key, i]));
    return mode.features.map((id) => this.ctx.byId.get(id)!).filter(Boolean)
      .sort((a, b) => (typeOrder.get(a.type) ?? 99) - (typeOrder.get(b.type) ?? 99) || b.count - a.count);
  }

  private renderDetail(): void {
    clear(this.detail);
    const mode = this.current();
    if (!mode) {
      this.detail.append(h('div', {class: 'oa-empty'}, 'No binding modes: no frames were analysed.'));
      return;
    }
    const features = this.modeFeatures(mode);
    const rep = this.ctx.frames[mode.representative];
    const all = features.map((f) => f.id);

    clear(this.chipBox);
    this.chips.clear();
    for (const f of features) {
      const chip = this.chip(f);
      this.chips.set(f.id, chip);
      this.chipBox.append(chip);
    }
    if (!features.length) this.chipBox.append(h('span', {class: 'oa-muted'}, 'No interaction above the threshold in these frames.'));

    clear(this.typeBox);
    for (const type of TYPES) {
      const ids = features.filter((f) => f.type === type.key).map((f) => f.id);
      if (!ids.length) continue;
      const chip = h('button', {class: 'oa-chip', dataset: {type: type.key}, title: `Switch all ${type.label} interactions of this mode`,
        onclick: () => this.ctx.toggleHidden(ids)},
      h('span', {class: 'sw', style: {background: type.color}}), type.label, h('span', {class: 'n'}, String(ids.length)));
      this.typeBox.append(chip);
    }

    const seg = h('div', {class: 'oa-seg'});
    for (const [id, label] of [['depiction', 'Depiction'], ['network', 'ProLIF network']] as [View2D, string][])
      seg.append(h('button', {class: id === this.view2d ? 'active' : '', onclick: () => this.setView(id)}, label));

    clear(this.view2dHost);
    this.view2dHost.append(this.view2d === 'network' ? this.network.el : this.ligand.el);
    const toggles = h('div', {class: 'oa-mode-toggles'},
      h('div', {class: 'oa-mode-toggles-head'},
        h('span', {class: 'oa-h4', style: {margin: '0'}}, 'Interactions'), this.chipCount,
        h('span', {style: {flex: '1'}}),
        h('button', {class: 'oa-btn oa-btn-sm', onclick: () => this.ctx.toggleHidden(all, false)}, 'All on'),
        h('button', {class: 'oa-btn oa-btn-sm', onclick: () => this.ctx.toggleHidden(all, true)}, 'All off')),
      this.typeBox, this.chipBox,
      h('div', {class: 'oa-small oa-muted'}, 'Click to switch an interaction on or off, double-click to show it alone. The 2D, network and 3D views follow.'));
    this.views.replaceChildren(h('div', {class: 'oa-mode-2d'}, this.view2dHost), this.slot3d, toggles);

    this.detail.append(
      h('div', {class: 'oa-mode-detail-head'},
        h('div', null, h('div', {class: 'oa-title', style: {fontSize: '18px'}}, mode.label),
          h('div', {class: 'oa-muted oa-small'},
            `Rank ${mode.rank + 1} · ${pct(mode.percent)} of frames (${mode.count}) · first at frame ${this.ctx.frames[mode.first]} · representative frame ${rep}`)),
        h('div', {style: {display: 'flex', gap: '6px', flexWrap: 'wrap'}},
          seg,
          h('button', {class: 'oa-btn oa-btn-sm', onclick: () => this.ctx.selectMode(mode)}, 'Representative frame'),
          h('button', {class: 'oa-btn oa-btn-sm', title: 'Next frame of this mode after the current one',
            onclick: () => this.nextOccurrence(mode)}, 'Next occurrence'),
          exportButtons({
            ...(this.view2d === 'depiction' ? {'2D SVG': () => this.export2D(mode, 'svg')} : {}),
            '2D PNG': () => this.export2D(mode, 'png'),
            'Modes CSV': () => this.exportModes(),
          }))),
      this.views,
      h('div', null, h('div', {class: 'oa-small oa-muted', style: {marginBottom: '4px'}}, 'Where this mode occurs (click to show a frame)'),
        this.timeline, frameAxis(this.ctx.frames)));
    if (this.ctx.selectedMode?.index !== mode.index)
      this.ctx.selectMode(mode, this.ctx.modeAt(this.ctx.frame)?.index !== mode.index);
    this.refreshToggles();
    this.drawStrips();
  }

  private chip(f: FeatureView): HTMLElement {
    const type = typeInfo(f.type);
    const chip = h('button', {class: 'oa-chip', title: `${f.id}\nPresent in ${pct((100 * f.count) / this.ctx.nFrames)} of all frames`},
      h('span', {class: 'sw', style: {background: type.color}}), f.label, h('span', {class: 'n'}, f.sublabel));
    let timer = 0;
    chip.addEventListener('click', () => {
      clearTimeout(timer);
      timer = window.setTimeout(() => this.ctx.toggleHidden([f.id]), 220);
    });
    chip.addEventListener('dblclick', () => {
      clearTimeout(timer);
      const mode = this.current();
      if (!mode) return;
      this.ctx.toggleHidden(mode.features.filter((id) => id !== f.id), true);
      this.ctx.toggleHidden([f.id], false);
    });
    chip.addEventListener('mouseenter', () => this.hover([f.id]));
    chip.addEventListener('mouseleave', () => this.hover([]));
    return chip;
  }

  private hover(ids: string[]): void {
    this.emphasis = ids;
    this.ctx.hoverFeatures(ids);
    this.paint2D();
  }

  private paint2D(): void {
    const mode = this.current();
    if (!mode) return;
    const features = this.modeFeatures(mode);
    if (this.view2d === 'network') this.network.show(features, this.ctx.hidden, this.emphasis);
    else this.ligand.update(features, this.ctx.hidden, this.emphasis);
  }

  private refreshToggles(): void {
    const mode = this.current();
    if (!mode) return;
    for (const [id, chip] of this.chips) chip.classList.toggle('off', this.ctx.hidden.has(id));
    for (const chip of this.typeBox.querySelectorAll<HTMLElement>('.oa-chip')) {
      const ids = mode.features.filter((id) => this.ctx.byId.get(id)?.type === chip.dataset.type);
      chip.classList.toggle('off', ids.every((id) => this.ctx.hidden.has(id)));
    }
    const on = mode.features.filter((id) => !this.ctx.hidden.has(id)).length;
    this.chipCount.textContent = `${on} of ${mode.features.length} on`;
    this.paint2D();
  }

  private setView(view: View2D): void {
    this.view2d = view;
    try {
      localStorage.setItem(VIEW_KEY, view);
    } catch {
      // Remembering the 2D style is a convenience only.
    }
    this.renderDetail();
  }

  private async export2D(mode: BindingMode, format: 'svg' | 'png'): Promise<void> {
    const name = safeName(`${this.ctx.summary.name}_${this.ctx.schemaId}_${mode.label}_${this.view2d}`);
    if (this.view2d === 'network') {
      const img = new Image();
      img.src = await this.network.image();
      await img.decode();
      download(await canvasToPng(img.width, img.height, (g) => g.drawImage(img, 0, 0)), `${name}.png`);
      return;
    }
    const svg = this.ligand.el.querySelector('svg');
    if (!svg) throw new Error('the ligand is not drawn yet');
    if (format === 'svg') download(new Blob([standaloneSvg(svg).text], {type: 'image/svg+xml'}), `${name}.svg`);
    else download(await svgToPng(svg), `${name}.png`);
  }

  private exportModes(): void {
    downloadCsv([
      ['rank', 'binding_mode', 'frames', 'percent', 'first_frame', 'representative_frame', 'interactions', 'frame_list'],
      ...this.ctx.modes.modes.map((m) => [m.rank + 1, m.label, m.count, m.percent.toFixed(2), this.ctx.frames[m.first],
        this.ctx.frames[m.representative],
        m.features.map((id) => {
          const f = this.ctx.byId.get(id);
          return f ? `${f.label} (${f.sublabel})` : id;
        }).join('; '),
        m.positions.map((p) => this.ctx.frames[p]).join(' ')]),
    ], `${this.ctx.summary.name}_${this.ctx.schemaId}_binding_modes_threshold${this.ctx.threshold}`);
  }

  private nextOccurrence(mode: BindingMode): void {
    const now = this.ctx.position(this.ctx.frame);
    const next = mode.positions.find((p) => p > now) ?? mode.positions[0];
    if (this.ctx.selectedMode?.index !== mode.index) this.ctx.selectMode(mode, false);
    this.ctx.setFrame(this.ctx.frames[next]);
  }

  dispose(): void {
    this.resize.disconnect();
    this.network.dispose();
    this.disposers.forEach((d) => d());
  }
}
