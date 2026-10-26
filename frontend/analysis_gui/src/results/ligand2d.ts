import {api, LigandDepiction} from '../api';
import {typeInfo} from '../colors';
import {clear, h, hideTooltip, s, showTooltip} from '../dom';
import {FeatureView} from './context';

const cache = new Map<string, Promise<LigandDepiction>>();

export function depiction(root: string, ligand: string): Promise<LigandDepiction> {
  const key = `${root}|${ligand}`;
  if (!cache.has(key)) {
    const promise = api.ligand2d(root, ligand);
    promise.catch(() => cache.delete(key));
    cache.set(key, promise);
  }
  return cache.get(key)!;
}

type Atom = LigandDepiction['atoms'][number];

function featuresByAtom(features: FeatureView[]): Map<number, FeatureView[]> {
  const bySerial = new Map<number, FeatureView[]>();
  for (const f of features)
    for (const serial of f.ligandAtoms) {
      if (!bySerial.has(serial)) bySerial.set(serial, []);
      bySerial.get(serial)!.push(f);
    }
  return bySerial;
}

/** The depiction SVG with an empty halo layer behind the drawing. */
function parseDepiction(data: LigandDepiction): {svg: SVGSVGElement; halos: SVGGElement} {
  const holder = h('div');
  holder.innerHTML = data.svg;
  const svg = holder.querySelector('svg')!;
  svg.removeAttribute('width');
  svg.removeAttribute('height');
  if (!svg.getAttribute('viewBox')) svg.setAttribute('viewBox', `0 0 ${data.width} ${data.height}`);
  const halos = s('g', {class: 'halos'});
  const background = svg.querySelector('rect');
  if (background?.nextSibling) svg.insertBefore(halos, background.nextSibling);
  else svg.prepend(halos);
  return {svg, halos};
}

/** One coloured halo per interaction type on each interacting atom; switched-off ones as dashed rings. */
function drawHalos(layer: SVGGElement, data: LigandDepiction, features: FeatureView[],
  hidden: Set<string> = new Set(), emphasized: Set<string> = new Set()): void {
  clear(layer);
  const bySerial = featuresByAtom(features);
  for (const atom of data.atoms) {
    const atFeatures = bySerial.get(atom.serial);
    if (!atFeatures) continue;
    const types = [...new Set(atFeatures.filter((f) => !hidden.has(f.id)).map((f) => f.type))];
    if (!types.length)
      layer.append(s('circle', {cx: atom.x, cy: atom.y, r: 13, fill: 'none', stroke: '#94a3b8',
        'stroke-width': 1.5, 'stroke-dasharray': '3 3'}));
    types.forEach((type, k) => {
      layer.append(s('circle', {cx: atom.x, cy: atom.y, r: Math.max(15 - k * 4.5, 5),
        fill: typeInfo(type).halo, 'fill-opacity': 0.75}));
    });
    if (atFeatures.some((f) => emphasized.has(f.id)))
      layer.append(s('circle', {cx: atom.x, cy: atom.y, r: 17.5, fill: 'none', stroke: '#0f172a', 'stroke-width': 2}));
  }
}

/** A standalone depiction with the halos of these interactions, for small multiples (Markov nodes, tooltips). */
export function ligandThumbnail(data: LigandDepiction, features: FeatureView[]): SVGSVGElement {
  const {svg, halos} = parseDepiction(data);
  drawHalos(halos, data, features);
  // Crop to the molecule (with room for halos and atom labels) so small copies stay legible.
  const pad = 22;
  const xs = data.atoms.map((a) => a.x);
  const ys = data.atoms.map((a) => a.y);
  if (xs.length) {
    const x0 = Math.max(0, Math.min(...xs) - pad);
    const y0 = Math.max(0, Math.min(...ys) - pad);
    svg.setAttribute('viewBox', `${x0} ${y0} ${Math.min(data.width, Math.max(...xs) + pad) - x0} ${Math.min(data.height, Math.max(...ys) + pad) - y0}`);
  }
  return svg;
}

/** RDKit depiction of the ligand with interaction halos per atom; click an atom to switch its interactions. */
export class Ligand2D {
  readonly el = h('div', {class: 'oa-lig2d'});
  private halos: SVGGElement | null = null;
  private data: LigandDepiction | null = null;
  private features: FeatureView[] = [];
  private hidden = new Set<string>();
  private emphasized = new Set<string>();

  constructor(root: string, ligand: string | null, private onToggle: (ids: string[]) => void,
    private onHover: (ids: string[]) => void) {
    if (!ligand) {
      this.el.append(h('div', {class: 'oa-empty'}, 'No small-molecule ligand (peptide analysis).'));
      return;
    }
    this.el.append(h('div', {class: 'oa-empty', style: {border: '0'}}, h('div', {class: 'oa-spinner'}), 'Drawing ligand…'));
    depiction(root, ligand).then((data) => {
      this.data = data;
      clear(this.el);
      const {svg, halos} = parseDepiction(data);
      this.halos = halos;
      this.el.append(svg, h('div', {class: 'oa-lig2d-legend'}, 'Click an atom to switch its interactions on or off'));
      svg.addEventListener('mousemove', (e) => this.hover(e, svg));
      svg.addEventListener('mouseleave', () => { hideTooltip(); this.onHover([]); svg.style.cursor = ''; });
      svg.addEventListener('click', (e) => {
        const atom = this.atomAt(e, svg);
        const ids = atom ? (this.atomFeatures().get(atom.serial) ?? []).map((f) => f.id) : [];
        if (ids.length) this.onToggle(ids);
      });
      this.paint();
    }).catch((e: Error) => {
      clear(this.el);
      this.el.append(h('div', {class: 'oa-empty'}, `Could not draw the ligand: ${e.message}`));
    });
  }

  update(features: FeatureView[], hidden: Set<string>, emphasized: string[] = []): void {
    this.features = features;
    this.hidden = hidden;
    this.emphasized = new Set(emphasized);
    this.paint();
  }

  private atomFeatures(): Map<number, FeatureView[]> {
    return featuresByAtom(this.features);
  }

  private paint(): void {
    if (this.halos && this.data) drawHalos(this.halos, this.data, this.features, this.hidden, this.emphasized);
  }

  private atomAt(event: MouseEvent, svg: SVGSVGElement): Atom | null {
    if (!this.data) return null;
    const rect = svg.getBoundingClientRect();
    const scale = this.data.width / rect.width;
    const x = (event.clientX - rect.left) * scale;
    const y = (event.clientY - rect.top) * scale;
    let best: Atom | null = null;
    let bestDist = 18 * 18;
    for (const atom of this.data.atoms) {
      const d = (atom.x - x) ** 2 + (atom.y - y) ** 2;
      if (d < bestDist) {
        bestDist = d;
        best = atom;
      }
    }
    return best;
  }

  private hover(event: MouseEvent, svg: SVGSVGElement): void {
    const atom = this.atomAt(event, svg);
    const features = atom ? this.atomFeatures().get(atom.serial) ?? [] : [];
    svg.style.cursor = features.length ? 'pointer' : '';
    this.onHover(features.filter((f) => !this.hidden.has(f.id)).map((f) => f.id));
    if (!atom) {
      hideTooltip();
      return;
    }
    showTooltip(event, h('div', null, h('b', null, `${atom.symbol} · atom ${atom.serial}`),
      ...features.map((f) => h('div', {style: {opacity: this.hidden.has(f.id) ? '0.55' : '1'}},
        `${typeInfo(f.type).label} · ${f.label}${this.hidden.has(f.id) ? ' (off)' : ''}`)),
      features.length ? h('div', {class: 'oa-muted'}, 'Click to switch on/off') : null));
  }
}
