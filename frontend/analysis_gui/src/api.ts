export interface OptionField {
  dest: string;
  flag: string;
  default: unknown;
  help: string;
  required: boolean;
  choices: string[] | null;
  label: string;
  group: string;
  kind: 'file' | 'text' | 'bool' | 'int' | 'float' | 'choice';
  ext?: string[];
  min?: number;
  max?: number;
}

export interface OptionSchema {
  fields: OptionField[];
  groups: {id: string; label: string}[];
}

export interface CommandPreview {
  argv: string[];
  command: string;
  errors: Record<string, string>;
  warnings: string[];
}

export interface DirEntry {name: string; path: string; dir: boolean; size: number | null}
export interface DirListing {path: string; parent: string | null; entries: DirEntry[]; isResults: boolean}
export interface Place {label: string; path: string}

export interface TopologyInfo {
  supported: boolean;
  atoms?: number;
  waters?: number;
  ions?: number;
  ligands: {name: string; copies: number; heavyAtoms: number; chains: string[]}[];
  suggestedLigand?: string | null;
  chains: {id: string; residues: number; nucleic: boolean}[];
  sdfNearby?: string[];
}

export interface SimulationCandidate {label: string; topology: string; trajectory: string; ligands: string[]}

export interface RunSnapshot {
  id: string;
  status: 'running' | 'finished' | 'failed' | 'cancelled';
  returncode: number | null;
  command: string;
  outputDir: string;
  elapsed: number;
  /** Form values of the run (older servers do not send them). */
  values?: Record<string, unknown>;
  lines: string[];
  /** Seconds since the start at which each line arrived. */
  times?: number[];
  next: number;
  progress: string;
}

export interface FileEntry {
  path: string;
  size: number;
  ext: string;
  image?: boolean;
  inline?: boolean;
  width?: number;
  height?: number;
}

export interface StableWaterLevel {
  name: string;
  minFrames: number;
  clusters: number;
  representative: string | null;
  residues: string | null;
}

export interface Summary {
  root: string;
  name: string;
  run: {command: string; values: Record<string, unknown>; status?: string; started?: number; ended?: number | null} | null;
  ligand: string | null;
  peptide: string | null;
  threshold: number;
  minTransition: number;
  schemas: {id: string; dir: string}[];
  rmsd: {frames: number[]; series: {name: string; values: (number | null)[]}[]} | null;
  rmsdMatrix: string | null;
  structure: {pdb: string; traj: string | null; frames?: number} | null;
  clouds: string | null;
  interactingWaters: number[];
  stableWaters: {dir: string; eps: string; levels: StableWaterLevel[]; csv: string | null}[];
  log: {path: string; messages: {level: string; text: string}[]} | null;
  files: FileEntry[];
}

export interface Residue {resnr: number; resname: string; chain: string}

export interface Feature {
  id: string;
  type: string;
  partner: string;
  ligandAtoms: number[];
  group: string | null;
  residues: Residue[];
  count: number;
  runs: number[];
}

export interface TopMode {
  mode: string;
  firstFrame: number | null;
  frames: number[];
  representativeFrame: number | null;
  percent: number;
}

export interface SchemaData {
  frames: number[];
  features: Feature[];
  topModes: TopMode[] | null;
  figures: {
    arranged: string | null;
    individual: Record<string, string>;
    markov: {threshold: number; path: string}[];
    ligandNumbering: string | null;
    modePdbs: Record<string, string>;
  };
}

export interface FrameInteraction {
  features: string[];
  ligand: [number, number, number];
  partner: [number, number, number];
  water: [number, number, number] | null;
}

export interface LigandDepiction {
  svg: string;
  width: number;
  height: number;
  atoms: {serial: number; symbol: string; x: number; y: number}[];
}

async function request<T>(url: string, init?: RequestInit): Promise<T> {
  const response = await fetch(url, init);
  const text = await response.text();
  let body: unknown = null;
  try {
    body = text ? JSON.parse(text) : null;
  } catch {
    throw new Error(`${response.status} ${response.statusText}: ${text.slice(0, 200)}`);
  }
  if (!response.ok) {
    const error = new Error((body as {error?: string})?.error ?? `${response.status} ${response.statusText}`);
    (error as Error & {body?: unknown}).body = body;
    throw error;
  }
  return body as T;
}

function query(params: Record<string, string | number>): string {
  return Object.entries(params).map(([k, v]) => `${k}=${encodeURIComponent(String(v))}`).join('&');
}

const post = (body: unknown): RequestInit => ({
  method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify(body),
});

export interface UploadResult {path: string; size: number; outputDir: string}

export const api = {
  /** Stream a dropped file to <outputDir>/inputs on the server (a default folder when empty). */
  upload(file: File, outputDir: string, onProgress: (fraction: number) => void): Promise<UploadResult> {
    return new Promise((resolve, reject) => {
      const xhr = new XMLHttpRequest();
      xhr.open('POST', `/api/upload?${query({name: file.name, outputDir})}`);
      xhr.upload.onprogress = (e) => e.lengthComputable && onProgress(e.loaded / e.total);
      xhr.onload = () => {
        let body: {error?: string} & Partial<UploadResult> = {};
        try {
          body = JSON.parse(xhr.responseText);
        } catch {
          // Non-JSON error page.
        }
        if (xhr.status === 200) resolve(body as UploadResult);
        else reject(new Error(body.error ?? `Upload of ${file.name} failed (${xhr.status})`));
      };
      xhr.onerror = () => reject(new Error(`Upload of ${file.name} failed`));
      xhr.send(file);
    });
  },
  options: () => request<OptionSchema>('/api/options'),
  command: (values: Record<string, unknown>, outputDir: string) =>
    request<CommandPreview>('/api/command', post({values, outputDir})),
  places: () => request<Place[]>('/api/fs/places'),
  list: (path: string) => request<DirListing>(`/api/fs/list?${query({path})}`),
  inspectTopology: (path: string) => request<TopologyInfo>(`/api/inspect/topology?${query({path})}`),
  inspectSimulation: (path: string) => request<SimulationCandidate[]>(`/api/inspect/simulation?${query({path})}`),
  startRun: (values: Record<string, unknown>, outputDir: string, copyTopology: boolean) =>
    request<RunSnapshot>('/api/runs', post({values, outputDir, copyTopology})),
  run: (id: string, since: number) => request<RunSnapshot>(`/api/runs/${id}?since=${since}`),
  cancelRun: (id: string) => request<RunSnapshot>(`/api/runs/${id}/cancel`, {method: 'POST'}),
  summary: (root: string) => request<Summary>(`/api/results/summary?${query({root})}`),
  schema: (root: string, dir: string) => request<SchemaData>(`/api/results/schema?${query({root, dir})}`),
  frameInteractions: (root: string, dir: string, frame: number) =>
    request<FrameInteraction[]>(`/api/results/frame-interactions?${query({root, dir, frame})}`),
  ligand2d: (root: string, ligand: string) =>
    request<LigandDepiction>(`/api/results/ligand2d?${query({root, ligand})}`),
  async lignetwork(root: string, dir: string, ligand: string, features: string[]): Promise<string> {
    const response = await fetch('/api/results/lignetwork', post({root, dir, ligand, features}));
    const text = await response.text();
    if (!response.ok) {
      let message = text;
      try {
        message = JSON.parse(text).error ?? text;
      } catch {
        // Plain-text error body.
      }
      throw new Error(message);
    }
    return text;
  },
  async coords(root: string, pdb: string, traj: string, frame: number): Promise<Float32Array> {
    const response = await fetch(`/api/results/coords?${query({root, pdb, traj, frame})}`);
    if (!response.ok) throw new Error(`Could not load frame ${frame}`);
    return new Float32Array(await response.arrayBuffer());
  },
  fileUrl: (root: string, path: string) => `/api/results/file?${query({root, path})}`,
  zipUrl: (root: string, path = '', figuresOnly = false) =>
    `/api/results/zip?${query({root, path, kind: figuresOnly ? 'figures' : 'all'})}`,
};
