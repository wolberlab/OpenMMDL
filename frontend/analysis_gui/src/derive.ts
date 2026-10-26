/** Binding modes and Markov transitions recomputed from df_all.csv presence data.
 * Rules follow openmmdl_analysis (BindingModeProcesser threshold, first-appearance labels,
 * MarkovChainAnalysis edge filter and node colouring). */

export function decodeRuns(runs: number[], n: number): Uint8Array {
  const flags = new Uint8Array(n);
  for (let i = 0; i < runs.length; i += 2) flags.fill(1, runs[i], runs[i] + runs[i + 1]);
  return flags;
}

export interface ModeInput {
  id: string;
  flags: Uint8Array;
  count: number;
}

export interface BindingMode {
  /** Order of first appearance, 0-based; the label is Binding_Mode_{index + 1}. */
  index: number;
  label: string;
  /** Rank by frame count, 0 = most frequent. */
  rank: number;
  features: string[];
  positions: number[];
  count: number;
  percent: number;
  first: number;
  /** Frame position whose full contact vector is closest to the mode's mean. */
  representative: number;
}

export interface ModeResult {
  passing: string[];
  /** Sorted by rank. */
  modes: BindingMode[];
  /** Mode index per frame position. */
  labels: Int32Array;
}

/** Same filter as the CLI: a feature counts when present in >= threshold% of the analysed frames. */
export function passingFeatures<T extends ModeInput>(features: T[], nFrames: number, thresholdPct: number): T[] {
  const limit = (thresholdPct / 100) * nFrames;
  return features.filter((f) => f.count >= limit);
}

export function deriveModes(features: ModeInput[], nFrames: number, thresholdPct: number): ModeResult {
  const passing = passingFeatures(features, nFrames, thresholdPct);
  const byKey = new Map<string, number>();
  const labels = new Int32Array(nFrames);
  const modes: BindingMode[] = [];
  for (let i = 0; i < nFrames; i++) {
    let key = '';
    for (const f of passing) key += f.flags[i] ? '1' : '0';
    let m = byKey.get(key);
    if (m === undefined) {
      m = modes.length;
      byKey.set(key, m);
      modes.push({
        index: m, label: `Binding_Mode_${m + 1}`, rank: 0,
        features: passing.filter((f) => f.flags[i]).map((f) => f.id),
        positions: [], count: 0, percent: 0, first: i, representative: i,
      });
    }
    labels[i] = m;
    modes[m].positions.push(i);
  }
  for (const mode of modes) {
    mode.count = mode.positions.length;
    mode.percent = (100 * mode.count) / nFrames;
    mode.representative = representative(mode.positions, features);
  }
  const ranked = [...modes].sort((a, b) => b.count - a.count || a.index - b.index);
  ranked.forEach((mode, rank) => (mode.rank = rank));
  return {passing: passing.map((f) => f.id), modes: ranked, labels};
}

/** Position minimising the L1 distance to the mean contact vector over all features. */
export function representative(positions: number[], features: ModeInput[]): number {
  if (positions.length <= 2) return positions[0];
  const mean = features.map((f) => {
    let sum = 0;
    for (const p of positions) sum += f.flags[p];
    return sum / positions.length;
  });
  let best = positions[0];
  let bestDist = Infinity;
  for (const p of positions) {
    let dist = 0;
    for (let j = 0; j < features.length; j++) dist += Math.abs(features[j].flags[p] - mean[j]);
    if (dist < bestDist - 1e-9) {
      bestDist = dist;
      best = p;
    }
  }
  return best;
}

export type ThirdColor = 'first' | 'second' | 'third' | 'even' | 'notTop10';

export interface MarkovNode {
  mode: number;
  occurrences: number;
  framesPct: number;
  selfCount: number;
  selfProbability: number;
  selfFramesPct: number;
  color: ThirdColor;
  top10: boolean;
  /** Has a kept self-loop. */
  selfLoop: boolean;
}

export interface MarkovEdge {
  from: number;
  to: number;
  count: number;
  /** Share of all frames, the quantity the min-transition threshold applies to. */
  framesPct: number;
  /** count / occurrences(from). */
  probability: number;
}

export interface MarkovResult {
  nodes: MarkovNode[];
  edges: MarkovEdge[];
  hiddenModes: number;
}

export function markov(labels: Int32Array, minTransitionPct: number): MarkovResult {
  const n = labels.length;
  const occurrences = new Map<number, number>();
  for (const label of labels) occurrences.set(label, (occurrences.get(label) ?? 0) + 1);

  const transitions = new Map<string, number>();
  const selfLoops = new Map<number, number>();
  for (let i = 0; i < n - 1; i++) {
    const a = labels[i];
    const b = labels[i + 1];
    if (a === b) selfLoops.set(a, (selfLoops.get(a) ?? 0) + 1);
    else transitions.set(`${a}>${b}`, (transitions.get(`${a}>${b}`) ?? 0) + 1);
  }

  const edgeMap = new Map<string, MarkovEdge>();
  const addEdge = (from: number, to: number) => {
    const key = `${from}>${to}`;
    if (edgeMap.has(key)) return;
    const count = transitions.get(key) ?? 0;
    edgeMap.set(key, {from, to, count, framesPct: (100 * count) / n,
      probability: (100 * count) / (occurrences.get(from) ?? 1)});
  };
  for (const [key, count] of transitions) {
    if ((100 * count) / n < minTransitionPct) continue;
    const [a, b] = key.split('>').map(Number);
    addEdge(a, b);
    addEdge(b, a);
  }
  const edges = [...edgeMap.values()].filter((e) => e.count > 0);

  const top10 = new Set(
    [...occurrences.entries()].sort((x, y) => y[1] - x[1] || x[0] - y[0]).slice(0, 10).map(([m]) => m));
  const part = Math.floor(n / 3);
  const end1 = part + (n % 3);
  const end2 = end1 + part;
  const thirds = new Map<number, [number, number, number]>();
  labels.forEach((label, i) => {
    const t = thirds.get(label) ?? [0, 0, 0];
    t[i < end1 ? 0 : i < end2 ? 1 : 2]++;
    thirds.set(label, t);
  });

  const inGraph = new Set<number>();
  for (const e of edges) {
    inGraph.add(e.from);
    inGraph.add(e.to);
  }
  const keptSelf = new Set<number>();
  for (const [mode, count] of selfLoops) {
    if ((100 * count) / n >= minTransitionPct) {
      keptSelf.add(mode);
      inGraph.add(mode);
    }
  }

  const nodes: MarkovNode[] = [...inGraph].map((mode) => {
    const occ = occurrences.get(mode) ?? 0;
    const selfCount = selfLoops.get(mode) ?? 0;
    const [p1, p2, p3] = (thirds.get(mode) ?? [0, 0, 0]).map((c) => c / Math.max(occ, 1));
    const color: ThirdColor = !top10.has(mode) ? 'notTop10' :
      p1 > 0.5 ? 'first' : p2 > 0.5 ? 'second' : p3 > 0.5 ? 'third' : 'even';
    return {
      mode, occurrences: occ, framesPct: (100 * occ) / n, selfCount,
      selfProbability: (100 * selfCount) / Math.max(occ, 1), selfFramesPct: (100 * selfCount) / n,
      color, top10: top10.has(mode), selfLoop: keptSelf.has(mode),
    };
  }).sort((a, b) => b.occurrences - a.occurrences);
  return {nodes, edges, hiddenModes: occurrences.size - nodes.length};
}
