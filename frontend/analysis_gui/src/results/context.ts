import {Feature, SchemaData, Summary} from '../api';
import {typeInfo} from '../colors';
import {BindingMode, decodeRuns, deriveModes, markov, MarkovResult, ModeResult} from '../derive';

export interface FeatureView extends Feature {
  flags: Uint8Array;
  label: string;
  sublabel: string;
}

type Events = {
  frame: number;
  features: Set<string>;
  hidden: Set<string>;
  mode: BindingMode | null;
  modes: ModeResult;
  markov: MarkovResult;
  hoverFeatures: string[];
};

type Listener<T> = (value: T) => void;

function residueLabel(f: Feature): string {
  if (!f.residues.length) return f.partner;
  const names = f.residues.map((r) => `${r.resname}${r.resnr}${r.chain && r.chain !== 'nan' ? `.${r.chain}` : ''}`);
  return names.length > 3 ? `${names.slice(0, 3).join(', ')} +${names.length - 3}` : names.join(', ');
}

function atomLabel(f: Feature): string {
  if (!f.ligandAtoms.length) return f.partner === 'ligand' ? 'ligand' : '';
  const atoms = f.ligandAtoms.length > 3 ? `${f.ligandAtoms.slice(0, 3).join(',')}…` : f.ligandAtoms.join(',');
  return `atom${f.ligandAtoms.length > 1 ? 's' : ''} ${atoms}`;
}

/** Shared state of one opened results folder and schema. Panels read it and listen to changes. */
export class ResultsContext {
  readonly features: FeatureView[];
  readonly byId: Map<string, FeatureView>;
  readonly frames: number[];
  readonly nFrames: number;
  modes: ModeResult;
  markov: MarkovResult;
  threshold: number;
  minTransition: number;
  frame: number;
  selected = new Set<string>();
  /** Interactions switched off by the user; hidden in 2D, network and 3D. */
  hidden = new Set<string>();
  selectedMode: BindingMode | null = null;
  private listeners = new Map<keyof Events, Set<Listener<unknown>>>();

  constructor(readonly summary: Summary, readonly schemaId: string, readonly schemaDir: string,
    readonly data: SchemaData) {
    this.frames = data.frames;
    this.nFrames = data.frames.length;
    const byResidue = schemaId !== 'ligand';
    this.features = data.features.map((f) => {
      const type = typeInfo(f.type).label + (f.group ? ` · ${f.group}` : '');
      const label = byResidue ? residueLabel(f) : `Ligand ${atomLabel(f)}`;
      const sublabel = byResidue ? `${type} · ${atomLabel(f)}` : `${type} · ${residueLabel(f)}`;
      return {...f, flags: decodeRuns(f.runs, this.nFrames), label, sublabel};
    });
    this.byId = new Map(this.features.map((f) => [f.id, f]));
    this.threshold = summary.threshold;
    this.minTransition = summary.minTransition;
    this.modes = deriveModes(this.features, this.nFrames, this.threshold);
    this.markov = markov(this.modes.labels, this.minTransition);
    this.frame = this.frames[0] ?? 0;
  }

  get root(): string {
    return this.summary.root;
  }

  on<K extends keyof Events>(event: K, listener: Listener<Events[K]>): () => void {
    if (!this.listeners.has(event)) this.listeners.set(event, new Set());
    this.listeners.get(event)!.add(listener as Listener<unknown>);
    return () => this.listeners.get(event)?.delete(listener as Listener<unknown>);
  }

  private emit<K extends keyof Events>(event: K, value: Events[K]): void {
    this.listeners.get(event)?.forEach((listener) => listener(value));
  }

  /** Index of a trajectory frame in the analysed frame list, or -1. */
  position(frame: number): number {
    let lo = 0;
    let hi = this.frames.length - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      if (this.frames[mid] === frame) return mid;
      if (this.frames[mid] < frame) lo = mid + 1;
      else hi = mid - 1;
    }
    return -1;
  }

  /** Binding mode of a trajectory frame, if it was analysed. */
  modeAt(frame: number): BindingMode | null {
    const p = this.position(frame);
    if (p < 0) return null;
    const index = this.modes.labels[p];
    return this.modes.modes.find((m) => m.index === index) ?? null;
  }

  setFrame(frame: number): void {
    if (frame === this.frame) return;
    this.frame = frame;
    this.emit('frame', frame);
  }

  selectFeatures(ids: string[], additive = false): void {
    if (!additive) this.selected.clear();
    for (const id of ids) {
      if (additive && this.selected.has(id)) this.selected.delete(id);
      else this.selected.add(id);
    }
    this.emit('features', this.selected);
  }

  hoverFeatures(ids: string[]): void {
    this.emit('hoverFeatures', ids);
  }

  /** Switch interactions off (or back on); with no explicit state, flips them as a group. */
  toggleHidden(ids: string[], hide?: boolean): void {
    const target = hide ?? !ids.every((id) => this.hidden.has(id));
    for (const id of ids) {
      if (target) this.hidden.add(id);
      else this.hidden.delete(id);
    }
    this.emit('hidden', this.hidden);
  }

  selectMode(mode: BindingMode | null, showFrame = true): void {
    this.selectedMode = mode;
    this.emit('mode', mode);
    if (mode) {
      this.selected = new Set(mode.features);
      this.emit('features', this.selected);
      if (showFrame) this.setFrame(this.frames[mode.representative]);
    }
  }

  setThreshold(threshold: number): void {
    this.threshold = threshold;
    const previous = this.selectedMode;
    this.modes = deriveModes(this.features, this.nFrames, threshold);
    this.emit('modes', this.modes);
    this.markov = markov(this.modes.labels, this.minTransition);
    this.emit('markov', this.markov);
    if (previous) this.selectMode(this.modeAt(this.frames[previous.representative]), false);
  }

  setMinTransition(value: number): void {
    this.minTransition = value;
    this.markov = markov(this.modes.labels, value);
    this.emit('markov', this.markov);
  }

  /** Features above the binding-mode threshold. */
  passing(): Set<string> {
    return new Set(this.modes.passing);
  }

  /** What the 3D view shows: the active selection, or else every passing interaction left on. */
  visible(): Set<string> {
    const base = this.selected.size ? this.selected : this.passing();
    return new Set([...base].filter((id) => !this.hidden.has(id)));
  }
}
