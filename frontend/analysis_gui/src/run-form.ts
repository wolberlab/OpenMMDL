import {api, CommandPreview, OptionField, OptionSchema, RunSnapshot, TopologyInfo} from './api';
import {clear, formatBytes, h, icon, toast} from './dom';
import {pickPath} from './picker';

interface FormState {
  values: Record<string, unknown>;
  outputDir: string;
  copyTopology: boolean;
  advancedOpen: boolean;
}

/** Shown up front; everything else lives under "Advanced options". */
const ESSENTIAL = ['topology', 'trajectory', 'ligand_sdf', 'ligand_name'];
const ADVANCED_GROUPS: Record<string, string> = {
  ligand: 'Ligand and receptor',
  analysis: 'Analysis',
  outputs: 'Extra outputs',
  advanced: 'Other',
};
/** Which field a dropped file fills, by extension. */
const DROP_TARGETS: [RegExp, string][] = [
  [/\.(pdb|tpr)$/i, 'topology'],
  [/\.(dcd|xtc|trr)$/i, 'trajectory'],
  [/\.sdf$/i, 'ligand_sdf'],
  [/\.csv$/i, 'dataframe'],
];

let schemaPromise: Promise<OptionSchema> | null = null;
let saved: FormState | null = null;

function defaults(schema: OptionSchema): Record<string, unknown> {
  const values: Record<string, unknown> = {};
  for (const f of schema.fields) values[f.dest] = f.kind === 'bool' ? Boolean(f.default) : f.default ?? '';
  return values;
}

function sep(path: string): string {
  return path.includes('\\') ? '\\' : '/';
}

function dirname(path: string): string {
  const i = Math.max(path.lastIndexOf('/'), path.lastIndexOf('\\'));
  return i > 0 ? path.slice(0, i) : path;
}

function basename(path: string): string {
  return path.split(/[\\/]/).pop() ?? path;
}

function highlightCommand(command: string): HTMLElement {
  const el = h('div', {class: 'oa-command'});
  for (const token of command.split(/(\s+)/)) {
    if (/^-{1,2}[a-z_]+$/i.test(token)) el.append(h('span', {class: 'flag'}, token));
    else el.append(token);
  }
  return el;
}

function isDefault(field: OptionField, value: unknown): boolean {
  if (field.kind === 'bool') return Boolean(value) === Boolean(field.default);
  if (value === '' || value === null || value === undefined) return field.default === null || field.default === '';
  return String(value) === String(field.default) || Number(value) === Number(field.default);
}

/** Highlight an element while files are dragged over it; hand dropped files to onFiles. */
function dropTarget(el: HTMLElement, onFiles: (files: File[]) => void): void {
  let depth = 0;
  el.addEventListener('dragenter', (e) => {
    e.preventDefault();
    depth++;
    el.classList.add('over');
  });
  el.addEventListener('dragover', (e) => {
    e.preventDefault();
    if (e.dataTransfer) e.dataTransfer.dropEffect = 'copy';
  });
  el.addEventListener('dragleave', () => {
    depth = Math.max(0, depth - 1);
    if (!depth) el.classList.remove('over');
  });
  el.addEventListener('drop', (e) => {
    e.preventDefault();
    e.stopPropagation();
    depth = 0;
    el.classList.remove('over');
    const files = [...(e.dataTransfer?.files ?? [])];
    if (files.length) onFiles(files);
  });
}

// Dropping a file next to a drop zone must not make the browser open it instead.
window.addEventListener('dragover', (e) => e.preventDefault());
window.addEventListener('drop', (e) => e.preventDefault());

/** Configure and start `openmmdl analysis`. */
export function runFormPage(onStarted: (run: RunSnapshot) => void, prefill: Record<string, unknown> | null): {el: HTMLElement; dispose(): void} {
  const el = h('div', null, h('div', {class: 'oa-card'}, h('div', {class: 'oa-empty', style: {border: '0'}},
    h('div', {class: 'oa-spinner'}), 'Loading options…')));
  schemaPromise ??= api.options();
  let disposed = false;
  schemaPromise.then((schema) => {
    if (!disposed) build(el, schema, onStarted, prefill);
  }).catch((e: Error) => {
    clear(el);
    el.append(h('div', {class: 'oa-note error'}, e.message));
  });
  return {el, dispose: () => { disposed = true; }};
}

function build(root: HTMLElement, schema: OptionSchema, onStarted: (run: RunSnapshot) => void,
  prefill: Record<string, unknown> | null): void {
  if (!saved || prefill)
    saved = {values: {...defaults(schema), ...(prefill ?? {})}, outputDir: '', copyTopology: true, advancedOpen: false};
  const state = saved;
  const fields = new Map(schema.fields.map((f) => [f.dest, f]));
  const errorEls = new Map<string, HTMLElement>();
  const inputs = new Map<string, HTMLInputElement | HTMLSelectElement>();
  const suggestions = new Map<string, HTMLElement>();
  const touched = new Set<string>();
  let preview: CommandPreview | null = null;
  let timer = 0;

  const commandHost = h('div');
  const notes = h('div', {class: 'oa-notes'});
  const changedBadge = h('span', {class: 'oa-badge oa-badge-muted'});
  const uploads = h('div', {class: 'oa-uploads'});
  const detected = h('div', {class: 'oa-detect'});
  const runBtn = h('button', {class: 'oa-btn oa-btn-primary', style: {width: '100%', padding: '10px'}, onclick: start},
    icon('run'), 'Run analysis');

  function set(dest: string, value: unknown, fromInput = false): void {
    state.values[dest] = value;
    touched.add(dest);
    const input = inputs.get(dest);
    if (input && !fromInput) {
      if (input instanceof HTMLInputElement && input.type === 'checkbox') input.checked = Boolean(value);
      else input.value = value === null || value === undefined ? '' : String(value);
    }
    if (dest === 'topology') onTopology(String(value ?? ''));
    if (dest === 'stable_water_analysis') syncWaterEps();
    schedule();
  }

  function setOutputDir(path: string): void {
    state.outputDir = path;
    const out = inputs.get('output_dir');
    if (out) out.value = path;
  }

  function schedule(): void {
    clearTimeout(timer);
    timer = window.setTimeout(refresh, 200);
  }

  async function refresh(): Promise<void> {
    updateChanged();
    try {
      preview = await api.command(state.values, state.outputDir);
    } catch (e) {
      notes.replaceChildren(h('div', {class: 'oa-note error'}, (e as Error).message));
      return;
    }
    clear(commandHost);
    commandHost.append(highlightCommand(preview.command));
    for (const [dest, errorEl] of errorEls) {
      const message = touched.has(dest) ? preview.errors[dest] ?? '' : '';
      errorEl.textContent = message;
      inputs.get(dest)?.classList.toggle('invalid', Boolean(message));
    }
    const missing = Object.keys(preview.errors).map((dest) =>
      dest === 'output_dir' ? 'output folder' : fields.get(dest)?.label.toLowerCase() ?? dest);
    notes.replaceChildren(
      ...preview.warnings.map((w) => h('div', {class: 'oa-note warn'}, w)),
      missing.length ? h('div', {class: 'oa-note info'}, `Still needed: ${missing.join(', ')}.`) :
        h('div', {class: 'oa-note ok'}, 'Ready to run.'));
    runBtn.disabled = missing.length > 0;
  }

  function updateChanged(): void {
    const changed = schema.fields.filter((f) => !ESSENTIAL.includes(f.dest) && !isDefault(f, state.values[f.dest])).length +
      (state.copyTopology ? 0 : 1);
    changedBadge.textContent = changed ? `${changed} changed` : '';
    changedBadge.style.display = changed ? '' : 'none';
  }

  async function start(): Promise<void> {
    runBtn.disabled = true;
    try {
      const run = await api.startRun(state.values, state.outputDir, state.copyTopology);
      onStarted(run);
    } catch (e) {
      toast((e as Error).message, 'error');
      runBtn.disabled = false;
    }
  }

  let lastTopology = '';
  async function onTopology(path: string): Promise<void> {
    if (!path || path === lastTopology) return;
    lastTopology = path;
    if (!state.outputDir) setOutputDir(`${dirname(path)}${sep(path)}openmmdl_analysis`);
    let info: TopologyInfo;
    try {
      info = await api.inspectTopology(path);
    } catch {
      return;
    }
    if (!info.supported) return;
    const ligandBox = suggestions.get('ligand_name')!;
    clear(ligandBox);
    for (const lig of info.ligands.slice(0, 6)) {
      ligandBox.append(h('button', {type: 'button', class: state.values.ligand_name === lig.name ? 'active' : '',
        title: `${lig.heavyAtoms} heavy atoms · ${lig.copies} cop${lig.copies > 1 ? 'ies' : 'y'} in chain ${lig.chains.join(', ')}`,
        onclick: () => { set('ligand_name', lig.name); lastTopology = ''; onTopology(path); }}, lig.name));
    }
    if (!state.values.ligand_name && info.suggestedLigand) set('ligand_name', info.suggestedLigand);
    const chainBox = suggestions.get('peptide');
    if (chainBox) {
      clear(chainBox);
      for (const chain of info.chains.filter((c) => c.id.trim() && c.residues <= 60))
        chainBox.append(h('button', {type: 'button', onclick: () => set('peptide', chain.id)}, `chain ${chain.id} (${chain.residues})`));
    }
    const sdfBox = suggestions.get('ligand_sdf')!;
    clear(sdfBox);
    for (const sdf of info.sdfNearby ?? [])
      sdfBox.append(h('button', {type: 'button', title: sdf, onclick: () => set('ligand_sdf', sdf)}, basename(sdf)));
    if (info.chains.length && info.chains.every((c) => c.nucleic) && !state.values.receptor_nucleic) set('receptor_nucleic', true);
  }

  function syncWaterEps(): void {
    const eps = inputs.get('water_eps');
    if (eps) eps.disabled = !state.values.stable_water_analysis;
  }

  /** Upload dropped files one by one; the first sets the output folder when none is chosen yet. */
  async function uploadFiles(files: File[], only?: string): Promise<void> {
    const plan: [File, string][] = [];
    for (const file of files) {
      const dest = only ?? DROP_TARGETS.find(([re]) => re.test(file.name))?.[1];
      const field = dest ? fields.get(dest) : undefined;
      if (!dest || !field || (field.ext && !field.ext.some((e) => file.name.toLowerCase().endsWith(e)))) {
        toast(`${file.name}: not a topology, trajectory, SDF or CSV file`, 'error');
        continue;
      }
      plan.push([file, dest]);
    }
    for (const [file, dest] of plan) {
      const bar = h('div', {style: {width: '0%'}});
      const status = h('span', {class: 'oa-small oa-muted'}, 'uploading…');
      const row = h('div', {class: 'oa-upload'},
        h('span', {class: 'name', title: file.name}, file.name),
        h('span', {class: 'oa-small oa-muted'}, `${fields.get(dest)!.label} · ${formatBytes(file.size)}`),
        h('div', {class: 'oa-progress'}, bar), status);
      uploads.append(row);
      try {
        const result = await api.upload(file, state.outputDir, (f) => {
          bar.style.width = `${Math.round(f * 100)}%`;
          status.textContent = `${Math.round(f * 100)}%`;
        });
        if (!state.outputDir) setOutputDir(result.outputDir);
        bar.style.width = '100%';
        status.textContent = '✓';
        row.classList.add('done');
        set(dest, result.path);
      } catch (e) {
        status.textContent = (e as Error).message;
        row.classList.add('failed');
      }
    }
  }

  function info(field: OptionField): HTMLElement {
    return h('span', {class: 'oa-info', title: `${field.help} (${field.flag})`}, 'ⓘ');
  }

  function fileField(field: OptionField): HTMLElement {
    const input = h('input', {class: 'oa-input oa-input-mono', spellcheck: false,
      placeholder: `Drop or browse ${field.ext?.join(', ')}`,
      value: String(state.values[field.dest] ?? ''), onchange: () => set(field.dest, input.value.trim(), true)});
    inputs.set(field.dest, input);
    const browse = h('button', {class: 'oa-btn', type: 'button', title: 'Pick a file on this computer (it is not copied)',
      onclick: async () => {
        const path = await pickPath({title: `Select ${field.label.toLowerCase()}`, mode: 'file', ext: field.ext,
          start: String(state.values[field.dest] || state.values.topology || '') || undefined});
        if (path) set(field.dest, path);
      }}, icon('folder'), 'Browse');
    const row = h('div', {class: 'oa-path-row oa-drop-field'}, input, browse);
    dropTarget(row, (files) => uploadFiles(files.slice(0, 1), field.dest));
    return row;
  }

  function widget(field: OptionField): HTMLElement {
    const value = state.values[field.dest];
    if (field.kind === 'file') return fileField(field);
    if (field.kind === 'bool') {
      const box = h('input', {type: 'checkbox', checked: Boolean(value), onchange: () => set(field.dest, box.checked, true)});
      inputs.set(field.dest, box);
      return h('label', {class: 'oa-check', title: `${field.help} (${field.flag})`}, box, field.label);
    }
    if (field.kind === 'choice') {
      const select = h('select', {class: 'oa-select', onchange: () => set(field.dest, select.value, true)},
        ...(field.choices ?? []).map((c) => h('option', {value: c, selected: String(value) === c}, c)));
      inputs.set(field.dest, select);
      return select;
    }
    const input = h('input', {class: 'oa-input', value: value === null || value === undefined ? '' : String(value),
      type: field.kind === 'int' || field.kind === 'float' ? 'number' : 'text',
      step: field.kind === 'float' ? 'any' : undefined, min: field.min, max: field.max,
      placeholder: field.dest === 'ligand_name' ? 'e.g. UNK' : field.kind === 'text' ? 'optional' : undefined,
      oninput: () => set(field.dest, field.dest === 'ligand_name' ? input.value.trim().toUpperCase() : input.value.trim(), true)});
    inputs.set(field.dest, input);
    return input;
  }

  function fieldEl(field: OptionField, wide = false): HTMLElement {
    const error = h('div', {class: 'oa-field-error'});
    errorEls.set(field.dest, error);
    const suggest = h('div', {class: 'oa-suggest'});
    suggestions.set(field.dest, suggest);
    if (field.kind === 'bool') return h('div', {class: 'oa-field'}, widget(field), error);
    return h('div', {class: `oa-field ${wide ? 'wide' : ''}`},
      h('div', {class: 'oa-field-label'}, field.label, field.required || field.dest === 'ligand_name' ? ' *' : '', info(field)),
      widget(field), suggest, error);
  }

  function dropZone(): HTMLElement {
    const picker = h('input', {type: 'file', multiple: true, style: {display: 'none'},
      accept: '.pdb,.tpr,.dcd,.xtc,.trr,.sdf,.csv',
      onchange: () => { uploadFiles([...(picker.files ?? [])]); picker.value = ''; }});
    const zone = h('div', {class: 'oa-drop'},
      icon('download'),
      h('div', {class: 'oa-drop-title'}, 'Drop your files here'),
      h('div', {class: 'oa-small oa-muted'}, 'Topology (.pdb, .tpr), trajectory (.dcd, .xtc, .trr) and ligand (.sdf) go to the right fields.'),
      h('div', {class: 'oa-drop-actions'},
        h('button', {class: 'oa-btn oa-btn-sm', type: 'button', onclick: () => picker.click()}, 'Choose files…'),
        h('button', {class: 'oa-btn oa-btn-sm', type: 'button', onclick: async () => {
          const path = await pickPath({title: 'Select an OpenMMDL simulation folder', mode: 'folder'});
          if (path) detect(path);
        }}, icon('folder'), 'Use a simulation folder…')),
      picker);
    dropTarget(zone, (files) => uploadFiles(files));
    return zone;
  }

  async function detect(path: string): Promise<void> {
    clear(detected);
    const candidates = await api.inspectSimulation(path).catch(() => []);
    if (!candidates.length) {
      detected.append(h('div', {class: 'oa-note warn'}, `No analysis-ready files found in ${basename(path)} (Final_Output or MD_Postprocessing).`));
      return;
    }
    for (const c of candidates) {
      const files = [c.topology, c.trajectory, ...c.ligands.slice(0, 1)].map(basename).join(' · ');
      detected.append(h('div', {class: 'oa-detect-item', title: [c.topology, c.trajectory, ...c.ligands].join('\n')},
        h('div', null, h('b', null, c.label), h('small', null, files)),
        h('button', {class: 'oa-btn oa-btn-sm oa-btn-primary', type: 'button', onclick: () => {
          if (!state.outputDir) setOutputDir(`${path}${sep(path)}Analysis`);
          set('trajectory', c.trajectory);
          if (c.ligands[0]) set('ligand_sdf', c.ligands[0]);
          set('topology', c.topology);
          clear(detected);
        }}, 'Use')));
    }
  }

  function outputField(): HTMLElement {
    const input = h('input', {class: 'oa-input oa-input-mono', value: state.outputDir, spellcheck: false,
      placeholder: 'Results go here (created if missing)',
      onchange: () => { state.outputDir = input.value.trim(); touched.add('output_dir'); schedule(); }});
    inputs.set('output_dir', input);
    const error = h('div', {class: 'oa-field-error'});
    errorEls.set('output_dir', error);
    return h('div', {class: 'oa-field wide'},
      h('div', {class: 'oa-field-label'}, 'Output folder *',
        h('span', {class: 'oa-info', title: 'All files and folders of the analysis are written here; type a new subfolder name to create it.'}, 'ⓘ')),
      h('div', {class: 'oa-path-row'}, input, h('button', {class: 'oa-btn', type: 'button', onclick: async () => {
        const path = await pickPath({title: 'Select the output folder', mode: 'folder', start: state.outputDir || undefined});
        if (path) {
          setOutputDir(path);
          touched.add('output_dir');
          schedule();
        }
      }}, icon('folder'), 'Browse')),
      error);
  }

  function advanced(): HTMLElement {
    const sections: HTMLElement[] = [];
    for (const [group, label] of Object.entries(ADVANCED_GROUPS)) {
      const groupFields = schema.fields.filter((f) => f.group === group && !ESSENTIAL.includes(f.dest)) as OptionField[];
      const bools = groupFields.filter((f) => f.kind === 'bool');
      const others = groupFields.filter((f) => f.kind !== 'bool');
      const extra = group === 'advanced' ? [h('label', {class: 'oa-check',
        title: 'The analysis rewrites its topology in place (water names, renumbering); the copy goes to <output>/inputs.'},
      h('input', {type: 'checkbox', checked: state.copyTopology,
        onchange: (e: Event) => { state.copyTopology = (e.target as HTMLInputElement).checked; updateChanged(); }}),
      'Work on a copy of the topology')] : [];
      if (!groupFields.length && !extra.length) continue;
      sections.push(h('div', {class: 'oa-adv-group'},
        h('div', {class: 'oa-adv-title'}, label),
        others.length ? h('div', {class: 'oa-fields'}, ...others.map((f) => fieldEl(f, f.kind === 'file'))) : null,
        bools.length || extra.length ? h('div', {class: 'oa-checks'}, ...bools.map((f) => fieldEl(f)), ...extra) : null));
    }
    const details = h('details', {class: 'oa-details oa-advanced', open: state.advancedOpen},
      h('summary', null, 'Advanced options ', changedBadge), ...sections);
    details.addEventListener('toggle', () => { state.advancedOpen = details.open; });
    return details;
  }

  const essentials = h('div', {class: 'oa-fields oa-essentials'},
    ...ESSENTIAL.map((dest) => fieldEl(fields.get(dest)!, fields.get(dest)!.kind === 'file')),
    outputField());

  const commandBox = h('details', {class: 'oa-details', style: {marginTop: '12px'}},
    h('summary', null, 'Command'),
    h('div', {style: {marginTop: '8px'}}, commandHost),
    h('button', {class: 'oa-btn oa-btn-sm', type: 'button', style: {marginTop: '6px'}, onclick: () => preview &&
      navigator.clipboard.writeText(preview.command).then(() => toast('Command copied'))}, icon('copy'), 'Copy'));

  clear(root);
  root.append(h('div', {class: 'oa-run-layout'},
    h('div', {class: 'oa-card'}, dropZone(), uploads, detected, essentials, advanced()),
    h('div', {class: 'oa-run-side'},
      h('div', {class: 'oa-card'}, notes, h('div', {style: {marginTop: '12px'}}, runBtn), commandBox))));
  syncWaterEps();
  if (state.values.topology) onTopology(String(state.values.topology));
  refresh();
}
