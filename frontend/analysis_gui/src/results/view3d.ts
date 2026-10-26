import * as NGL from 'ngl';
import {api, FrameInteraction} from '../api';
import {typeInfo} from '../colors';
import {clear, h, icon, toast} from '../dom';
import {download, safeName} from '../export';
import {ResultsContext} from './context';

type Rgb = [number, number, number];

interface CloudSet {
  coordinates: number[][];
  color: Rgb;
  radius: number;
}

const CLOUD_LABELS: Record<string, string> = {
  hydrophobic: 'Hydrophobic', acceptor: 'H-bond acceptor', donor: 'H-bond donor',
  waterbridge: 'Water bridge', negative_ionizable: 'Negative ionizable',
  positive_ionizable: 'Positive ionizable', pistacking: 'π-stacking', pication: 'π-cation',
  halogen: 'Halogen', metal: 'Metal',
};

function hexToRgb(hex: string): Rgb {
  const n = parseInt(hex.slice(1), 16);
  return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
}

function residueSele(residues: {resnr: number; chain: string}[]): string {
  const parts = residues.map((r) => (r.chain && r.chain !== 'nan' ? `${r.resnr}:${r.chain}` : `${r.resnr}`));
  return parts.length ? `(${[...new Set(parts)].join(' or ')}) and not hydrogen` : 'none';
}

const FRAME_CACHE = 64;

/** Docked NGL view of the analysed complex, following the results context. */
export class View3D {
  readonly el: HTMLElement;
  private stage: NGL.Stage | null = null;
  private comp: NGL.StructureComponent | null = null;
  private cartoonRep: NGL.RepresentationElement | null = null;
  private siteRep: NGL.RepresentationElement | null = null;
  private focusRep: NGL.RepresentationElement | null = null;
  private labelRep: NGL.RepresentationElement | null = null;
  private waterRep: NGL.RepresentationElement | null = null;
  private lineComp: NGL.Component | null = null;
  private cloudComps = new Map<string, NGL.Component>();
  private clouds: Record<string, CloudSet> | null = null;
  private extraComp: NGL.Component | null = null;
  private frameCache = new Map<number, Float32Array>();
  private interactionCache = new Map<number, FrameInteraction[]>();
  private wanted = -1;
  private shown = -1;
  private playing = false;
  private fps = 5;
  private disposers: (() => void)[] = [];
  private stageHost = h('div');
  private caption = h('div', {class: 'oa-3d-caption'});
  private message = h('div', {class: 'oa-3d-msg'});
  private slider = h('input', {type: 'range', min: 0, max: 0, step: 1, value: 0});
  private frameLabel = h('span', {class: 'frame'});
  private playBtn = h('button', {class: 'oa-btn oa-btn-sm oa-btn-icon', title: 'Play / pause'}, icon('play'));
  private showWaters = true;
  private showLines = true;
  private resizeObserver: ResizeObserver;
  private ligandSele: string;

  constructor(private ctx: ResultsContext) {
    this.ligandSele = ctx.summary.peptide ? `:${ctx.summary.peptide}` : `[${ctx.summary.ligand ?? 'UNK'}]`;
    const nFrames = ctx.summary.structure?.frames ?? 0;
    this.slider.max = String(Math.max(0, nFrames - 1));
    this.slider.addEventListener('input', () => ctx.setFrame(Number(this.slider.value)));
    this.playBtn.addEventListener('click', () => this.togglePlay());
    const speed = h('select', {class: 'oa-select', style: {width: 'auto', padding: '3px 6px'}, title: 'Playback speed',
      onchange: () => (this.fps = Number(speed.value))},
    ...[2, 5, 10, 20].map((fps) => h('option', {value: fps, selected: fps === this.fps}, `${fps} fps`)));
    const step = (delta: number) => () => ctx.setFrame(Math.min(Math.max(ctx.frame + delta, 0), nFrames - 1));

    this.el = h('div', {class: 'oa-card oa-dock'},
      h('div', {class: 'oa-3d-toolbar'},
        h('span', {class: 'oa-3d-title'}, '3D view'), h('span', {class: 'grow'}),
        this.toggleButton('Waters', () => this.showWaters, (v) => { this.showWaters = v; this.waterRep?.setVisibility(v); this.drawInteractions(); }),
        this.toggleButton('Contacts', () => this.showLines, (v) => { this.showLines = v; this.drawInteractions(); }),
        this.cloudMenu(),
        h('button', {class: 'oa-btn oa-btn-sm', title: 'Center on ligand', onclick: () => this.center()}, icon('reset'), 'Center'),
        h('button', {class: 'oa-btn oa-btn-sm', title: 'Save the 3D view as a PNG image', onclick: () => this.snapshot()},
          icon('download'), 'PNG')),
      h('div', {class: 'oa-3d-stage'}, this.stageHost, this.caption, this.message),
      nFrames > 1 ? h('div', {class: 'oa-3d-player'},
        h('button', {class: 'oa-btn oa-btn-sm oa-btn-icon', title: 'Previous frame', onclick: step(-1)}, icon('prev')),
        this.playBtn,
        h('button', {class: 'oa-btn oa-btn-sm oa-btn-icon', title: 'Next frame', onclick: step(1)}, icon('next')),
        this.slider, this.frameLabel, speed) : null);

    this.resizeObserver = new ResizeObserver(() => this.stage?.handleResize());
    this.resizeObserver.observe(this.stageHost);
    this.disposers.push(
      ctx.on('frame', (f) => this.showFrame(f)),
      ctx.on('features', () => { this.updateResidues(true); this.drawInteractions(); }),
      ctx.on('modes', () => { this.updateResidues(false); this.drawInteractions(); }),
      ctx.on('hidden', () => { this.updateResidues(false); this.drawInteractions(); }),
      ctx.on('hoverFeatures', (ids) => this.drawInteractions(ids)),
    );
  }

  private toggleButton(label: string, get: () => boolean, set: (v: boolean) => void): HTMLElement {
    const btn = h('button', {class: 'oa-btn oa-btn-sm'}, label);
    const sync = () => btn.classList.toggle('oa-btn-primary', get());
    btn.addEventListener('click', () => { set(!get()); sync(); });
    sync();
    return btn;
  }

  private cloudMenu(): HTMLElement {
    const wrap = h('div', {class: 'oa-menu'});
    if (!this.ctx.summary.clouds) return wrap;
    const pop = h('div', {class: 'oa-menu-pop', style: {display: 'none'}});
    const btn = h('button', {class: 'oa-btn oa-btn-sm', title: 'Interaction point clouds over all frames',
      onclick: async (e: MouseEvent) => {
        e.stopPropagation();
        const open = pop.style.display === 'none';
        pop.style.display = open ? 'flex' : 'none';
        if (open && !this.clouds) await this.fillCloudMenu(pop);
      }}, icon('layers'), 'Clouds');
    document.addEventListener('click', (e) => { if (!wrap.contains(e.target as Node)) pop.style.display = 'none'; });
    wrap.append(btn, pop);
    return wrap;
  }

  private async fillCloudMenu(pop: HTMLElement): Promise<void> {
    pop.append(h('div', {class: 'head'}, 'Loading…'));
    try {
      const response = await fetch(api.fileUrl(this.ctx.root, this.ctx.summary.clouds!));
      this.clouds = await response.json();
    } catch {
      clear(pop);
      pop.append(h('div', {class: 'head'}, 'Could not read clouds.json'));
      return;
    }
    clear(pop);
    pop.append(h('div', {class: 'head'}, 'Point clouds (all frames)'));
    for (const [key, cloud] of Object.entries(this.clouds!)) {
      if (!cloud.coordinates.length) continue;
      const [r, g, b] = cloud.color.map((c) => Math.round(c * 255));
      const box = h('input', {type: 'checkbox', onchange: () => this.setCloud(key, box.checked)});
      pop.append(h('label', null, box, h('span', {class: 'sw', style: {background: `rgb(${r},${g},${b})`}}),
        `${CLOUD_LABELS[key] ?? key}`, h('span', {class: 'oa-muted'}, ` ${cloud.coordinates.length}`)));
    }
  }

  private setCloud(key: string, visible: boolean): void {
    if (!this.stage || !this.clouds) return;
    const existing = this.cloudComps.get(key);
    if (existing) {
      existing.setVisibility(visible);
      return;
    }
    if (!visible) return;
    const cloud = this.clouds[key];
    const shape = new NGL.Shape(`cloud-${key}`);
    const radius = Math.max(cloud.radius, 0.15);
    for (const p of cloud.coordinates) shape.addSphere(p as Rgb, cloud.color, radius, CLOUD_LABELS[key] ?? key);
    const comp = this.stage.addComponentFromObject(shape) as NGL.Component;
    comp.addRepresentation('buffer', {opacity: 0.65});
    this.cloudComps.set(key, comp);
  }

  /** Create the NGL stage once the element is in the document. */
  async mount(): Promise<void> {
    const structure = this.ctx.summary.structure;
    if (!structure) {
      this.message.textContent = 'No structure in this results folder (Visualization/interacting_waters.pdb or complex.pdb).';
      return;
    }
    this.message.append(h('div', null, h('div', {class: 'oa-spinner'}), 'Loading structure…'));
    this.stage = new NGL.Stage(this.stageHost, {backgroundColor: 'white', tooltip: true});
    try {
      const comp = await this.stage.loadFile(api.fileUrl(this.ctx.root, structure.pdb),
        {ext: 'pdb', name: structure.pdb, defaultRepresentation: false}) as NGL.StructureComponent;
      this.comp = comp;
      this.cartoonRep = comp.addRepresentation('cartoon', {sele: 'protein or nucleic', color: '#a5b4fc', opacity: 0.9});
      comp.addRepresentation('licorice', {sele: `${this.ligandSele} and not hydrogen`, colorScheme: 'element',
        colorValue: '#f97316', multipleBond: 'symmetric', radiusSize: 0.22});
      this.siteRep = comp.addRepresentation('licorice', {sele: 'none', colorScheme: 'element',
        colorValue: '#cbd5e1', radiusSize: 0.12, opacity: 0.85});
      this.focusRep = comp.addRepresentation('licorice', {sele: 'none', colorScheme: 'element',
        colorValue: '#4f46e5', radiusSize: 0.2});
      this.labelRep = comp.addRepresentation('label', {sele: 'none', labelType: 'res', color: '#1e1b4b',
        showBackground: true, backgroundColor: '#ffffff', backgroundOpacity: 0.75, zOffset: 2, radiusSize: 1.1});
      this.waterRep = comp.addRepresentation('ball+stick', {sele: 'water and not hydrogen', color: '#38bdf8',
        radiusSize: 0.18, visible: this.showWaters});
      this.message.remove();
      this.updateResidues(false);
      this.center(0);
      await this.showFrame(this.ctx.frame);
    } catch (e) {
      clear(this.message);
      this.message.textContent = `Could not load the structure: ${(e as Error).message}`;
    }
  }

  async snapshot(): Promise<void> {
    if (!this.stage) return;
    try {
      // The stage may have been created while its tab was hidden; size it to what is on screen now.
      this.stage.handleResize();
      const blob = await this.stage.makeImage({factor: 2, antialias: true, trim: false, transparent: false}) as Blob;
      download(blob, `${safeName(`${this.ctx.summary.name}_frame${this.ctx.frame}`)}.png`);
    } catch (e) {
      toast(`Could not save the image: ${(e as Error).message}`, 'error');
    }
  }

  center(duration = 600): void {
    if (!this.comp) return;
    const site = this.siteSelection();
    this.comp.autoView(site === 'none' ? this.ligandSele : `${this.ligandSele} or ${site}`, duration);
  }

  private siteSelection(): string {
    const residues = [...this.ctx.visible()].flatMap((id) => this.ctx.byId.get(id)?.residues ?? []);
    return residueSele(residues);
  }

  private updateResidues(focus: boolean): void {
    if (!this.comp) return;
    const shown = (ids: Iterable<string>) => [...ids].filter((id) => !this.ctx.hidden.has(id))
      .flatMap((id) => this.ctx.byId.get(id)?.residues ?? []);
    this.siteRep?.setSelection(residueSele(shown(this.ctx.passing())));
    const selected = shown(this.ctx.selected);
    const focusSele = residueSele(selected);
    this.focusRep?.setSelection(focusSele);
    this.labelRep?.setSelection(selected.length ? `(${focusSele}) and .CA` : 'none');
    this.cartoonRep?.setParameters({opacity: selected.length ? 0.45 : 0.9});
    if (focus && selected.length) this.comp.autoView(`${this.ligandSele} or ${focusSele}`, 600);
  }

  async showFrame(frame: number): Promise<void> {
    const structure = this.ctx.summary.structure;
    this.slider.value = String(frame);
    this.updateCaption(frame);
    if (!this.comp || !structure?.traj) return;
    this.wanted = frame;
    let coords = this.frameCache.get(frame);
    if (!coords) {
      try {
        coords = await api.coords(this.ctx.root, structure.pdb, structure.traj, frame);
      } catch (e) {
        toast((e as Error).message, 'error');
        return;
      }
      this.frameCache.set(frame, coords);
      if (this.frameCache.size > FRAME_CACHE) this.frameCache.delete(this.frameCache.keys().next().value!);
    }
    if (this.wanted !== frame) return;
    if (coords.length !== this.comp.structure.atomCount * 3) {
      this.caption.textContent = 'Trajectory and structure atom counts differ; showing the static structure.';
      return;
    }
    this.comp.structure.updatePosition(coords);
    this.comp.updateRepresentations({position: true});
    this.shown = frame;
    await this.drawInteractions();
  }

  private updateCaption(frame: number): void {
    const n = this.ctx.summary.structure?.frames ?? 0;
    this.frameLabel.textContent = `Frame ${frame} / ${Math.max(n - 1, 0)}`;
    const mode = this.ctx.modeAt(frame);
    const parts = [`Frame ${frame}`];
    if (mode) parts.push(`${mode.label} (${mode.percent.toFixed(1)}%)`);
    else if (this.ctx.position(frame) < 0) parts.push('not analysed');
    const total = this.ctx.selected.size;
    if (total) {
      const shown = [...this.ctx.selected].filter((id) => !this.ctx.hidden.has(id)).length;
      parts.push(shown === total ? `${total} interaction${total > 1 ? 's' : ''}` : `${shown} of ${total} interactions on`);
    }
    this.caption.textContent = parts.join(' · ');
  }

  private async frameRows(frame: number): Promise<FrameInteraction[]> {
    let rows = this.interactionCache.get(frame);
    if (!rows) {
      rows = await api.frameInteractions(this.ctx.root, this.ctx.schemaDir, frame);
      this.interactionCache.set(frame, rows);
      if (this.interactionCache.size > FRAME_CACHE)
        this.interactionCache.delete(this.interactionCache.keys().next().value!);
    }
    return rows;
  }

  /** Dashed contact lines of the shown frame, limited to selected (or hovered, or passing) features. */
  async drawInteractions(hovered: string[] = []): Promise<void> {
    if (!this.stage) return;
    const frame = this.shown >= 0 ? this.shown : this.ctx.frame;
    this.updateCaption(this.ctx.frame);
    let rows: FrameInteraction[] = [];
    if (this.showLines && this.ctx.position(frame) >= 0) {
      try {
        rows = await this.frameRows(frame);
      } catch {
        rows = [];
      }
    }
    const visible = hovered.length ? new Set(hovered) : this.ctx.visible();
    const shape = new NGL.Shape('contacts');
    let count = 0;
    for (const row of rows) {
      const id = row.features.find((f) => visible.has(f));
      if (!id) continue;
      const feature = this.ctx.byId.get(id);
      const type = typeInfo(feature?.type ?? '');
      const color = hexToRgb(type.color);
      const name = feature ? `${feature.label} · ${feature.sublabel}` : id;
      if (row.water && this.showWaters) {
        dashed(shape, row.ligand, row.water, color, name);
        dashed(shape, row.water, row.partner, color, name);
        shape.addSphere(row.water, [0.22, 0.74, 0.97], 0.35, 'bridging water');
      } else {
        dashed(shape, row.ligand, row.partner, color, name);
      }
      count++;
    }
    if (this.lineComp) this.stage.removeComponent(this.lineComp);
    this.lineComp = null;
    if (count) {
      this.lineComp = this.stage.addComponentFromObject(shape) as NGL.Component;
      this.lineComp.addRepresentation('buffer', {});
    }
  }

  togglePlay(): void {
    this.playing = !this.playing;
    clear(this.playBtn);
    this.playBtn.append(icon(this.playing ? 'pause' : 'play'));
    if (this.playing) this.tick();
  }

  private tick(): void {
    if (!this.playing) return;
    const n = this.ctx.summary.structure?.frames ?? 0;
    const started = performance.now();
    this.ctx.setFrame((this.ctx.frame + 1) % Math.max(n, 1));
    const wait = Math.max(0, 1000 / this.fps - (performance.now() - started));
    setTimeout(() => this.tick(), wait);
  }

  /** Show an extra PDB (e.g. stable-water clusters) as spheres; null removes it. */
  async showExtra(path: string | null, color = '#ef4444'): Promise<void> {
    if (!this.stage) return;
    if (this.extraComp) this.stage.removeComponent(this.extraComp);
    this.extraComp = null;
    if (!path) return;
    const comp = await this.stage.loadFile(api.fileUrl(this.ctx.root, path), {ext: 'pdb', defaultRepresentation: false});
    if (!comp) return;
    this.extraComp = comp as NGL.Component;
    (comp as NGL.StructureComponent).addRepresentation('spacefill', {color, radiusScale: 0.45, opacity: 0.8});
    (comp as NGL.StructureComponent).autoView(800);
  }

  dispose(): void {
    this.playing = false;
    this.disposers.forEach((d) => d());
    this.resizeObserver.disconnect();
    this.stage?.dispose();
    this.stage = null;
  }
}

function dashed(shape: NGL.Shape, a: Rgb, b: Rgb, color: Rgb, name: string): void {
  const length = Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
  const dashes = Math.max(2, Math.round(length / 0.45));
  for (let i = 0; i < dashes; i += 2) {
    const t0 = i / dashes;
    const t1 = Math.min((i + 1) / dashes, 1);
    const p0: Rgb = [a[0] + (b[0] - a[0]) * t0, a[1] + (b[1] - a[1]) * t0, a[2] + (b[2] - a[2]) * t0];
    const p1: Rgb = [a[0] + (b[0] - a[0]) * t1, a[1] + (b[1] - a[1]) * t1, a[2] + (b[2] - a[2]) * t1];
    shape.addCylinder(p0, p1, color, 0.07, name);
  }
}
