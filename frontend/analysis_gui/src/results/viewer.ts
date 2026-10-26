import {api, SchemaData, Summary} from '../api';
import {clear, h, icon, toast} from '../dom';
import {ResultsContext} from './context';
import {FilesPanel} from './files';
import {InteractionsPanel} from './interactions';
import {MarkovPanel} from './markov';
import {ModesPanel} from './modes';
import {OverviewPanel} from './overview';
import {RmsdPanel} from './rmsd';
import {View3D} from './view3d';
import {WatersPanel} from './waters';

type Layout = 'docked' | 'full' | 'report';

interface Panel {
  el: HTMLElement;
  dispose(): void;
}

const TABS: {id: string; label: string}[] = [
  {id: 'overview', label: 'Overview'},
  {id: 'interactions', label: 'Interactions'},
  {id: 'modes', label: 'Binding modes'},
  {id: 'markov', label: 'Markov chain'},
  {id: 'rmsd', label: 'RMSD'},
  {id: 'waters', label: 'Waters'},
  {id: 'files', label: 'Files'},
];

const SCHEMA_LABELS: Record<string, string> = {residue: 'Residue', ligand: 'Ligand atom', default: 'Default'};
const LAYOUT_KEY = 'openmmdl-analysis-gui:layout';
const DOCK_KEY = 'openmmdl-analysis-gui:dock-width';

function stored(key: string, fallback: string): string {
  try {
    return localStorage.getItem(key) ?? fallback;
  } catch {
    return fallback;
  }
}

function store(key: string, value: string): void {
  try {
    localStorage.setItem(key, value);
  } catch {
    // Layout preferences are a convenience; ignore unavailable storage.
  }
}

export interface ViewerParams {
  root: string;
  schema?: string;
  tab?: string;
}

/** Results explorer for one analysis folder. */
export class ResultsViewer {
  readonly el = h('div');
  private ctx: ResultsContext | null = null;
  private view3d: View3D | null = null;
  private panels = new Map<string, Panel>();
  private layout: Layout = stored(LAYOUT_KEY, 'docked') as Layout;
  private tab: string;
  private schemaCache = new Map<string, SchemaData>();
  private summary: Summary | null = null;
  private body = h('div');
  private disposed = false;
  private onResize = () => this.fitHeight();

  constructor(private params: ViewerParams, private onRerun: (values: Record<string, unknown>) => void,
    private onState: (params: ViewerParams) => void) {
    this.tab = params.tab ?? 'overview';
    window.addEventListener('resize', this.onResize);
    this.el.append(h('div', {class: 'oa-card'}, h('div', {class: 'oa-empty', style: {border: '0'}},
      h('div', {class: 'oa-spinner'}), `Reading ${params.root}…`)));
    this.load().catch((e: Error) => {
      clear(this.el);
      this.el.append(h('div', {class: 'oa-card'}, h('div', {class: 'oa-note error'}, e.message)));
    });
  }

  private async load(): Promise<void> {
    this.summary = await api.summary(this.params.root);
    if (!this.summary.schemas.length) throw new Error('This folder has no binding-mode results (df_all.csv) yet.');
    const schema = this.summary.schemas.find((sc) => sc.id === this.params.schema) ?? this.summary.schemas[0];
    await this.openSchema(schema.id);
  }

  private async openSchema(id: string): Promise<void> {
    const summary = this.summary!;
    const schema = summary.schemas.find((sc) => sc.id === id)!;
    let data = this.schemaCache.get(id);
    if (!data) {
      data = await api.schema(summary.root, schema.dir);
      this.schemaCache.set(id, data);
    }
    if (this.disposed) return;
    const frame = this.ctx?.frame;
    this.teardown();
    this.ctx = new ResultsContext(summary, id, schema.dir, data);
    if (frame !== undefined) this.ctx.frame = frame;
    this.params = {...this.params, schema: id};
    this.view3d = new View3D(this.ctx);
    this.build();
    this.onState({...this.params, tab: this.tab});
    await this.view3d.mount();
  }

  private teardown(): void {
    this.panels.forEach((p) => p.dispose());
    this.panels.clear();
    this.view3d?.dispose();
    this.view3d = null;
  }

  private panel(id: string): Panel {
    let panel = this.panels.get(id);
    if (panel) return panel;
    const ctx = this.ctx!;
    const values = ctx.summary.run?.values;
    switch (id) {
    case 'overview': panel = new OverviewPanel(ctx, (tab) => this.showTab(tab), values ? () => this.onRerun(values) : null); break;
    case 'interactions': panel = new InteractionsPanel(ctx); break;
    case 'modes': panel = new ModesPanel(ctx); break;
    case 'markov': panel = new MarkovPanel(ctx); break;
    case 'rmsd': panel = new RmsdPanel(ctx); break;
    case 'waters': panel = new WatersPanel(ctx, (path) => this.view3d?.showExtra(path)); break;
    default: panel = new FilesPanel(ctx);
    }
    this.panels.set(id, panel);
    return panel;
  }

  private build(): void {
    const ctx = this.ctx!;
    const summary = ctx.summary;
    clear(this.el);

    const schemaSeg = h('div', {class: 'oa-seg', title: 'Residue: interactions keyed by protein residue and ligand atom. Ligand atom: keyed by ligand atom only, whichever residue it binds.'});
    for (const sc of summary.schemas)
      schemaSeg.append(h('button', {class: sc.id === ctx.schemaId ? 'active' : '',
        onclick: () => sc.id !== ctx.schemaId && this.openSchema(sc.id).catch((e: Error) => toast(e.message, 'error'))},
      SCHEMA_LABELS[sc.id] ?? sc.id));
    const layoutSeg = h('div', {class: 'oa-seg'});
    for (const [id, label] of [['docked', 'Docked 3D'], ['full', 'Full width'], ['report', 'Report']] as [Layout, string][])
      layoutSeg.append(h('button', {class: id === this.layout ? 'active' : '', onclick: () => this.setLayout(id)}, label));

    const meta = [summary.ligand ? `Ligand ${summary.ligand}` : summary.peptide ? `Peptide chain ${summary.peptide}` : null,
      `${ctx.nFrames} analysed frames`, summary.root].filter(Boolean).join(' · ');
    this.el.append(h('div', {class: 'oa-results-bar'},
      h('div', {class: 'grow'}, h('div', {class: 'oa-results-name'}, summary.name), h('div', {class: 'oa-results-meta'}, meta)),
      summary.schemas.length > 1 ? h('div', {class: 'oa-bar-group'}, h('span', {class: 'oa-bar-label'}, 'Partners'), schemaSeg) : null,
      h('div', {class: 'oa-bar-group'}, h('span', {class: 'oa-bar-label'}, 'Layout'), layoutSeg),
      h('button', {class: 'oa-btn oa-btn-sm', title: 'Read the folder again', onclick: () => this.reload()}, icon('reset'), 'Reload')));
    this.body = h('div');
    this.el.append(this.body);
    this.arrange();
  }

  private arrange(): void {
    const ctx = this.ctx!;
    clear(this.body);
    const grid = h('div', {class: `oa-results ${this.layout}`});
    grid.style.setProperty('--oa-dock-w', stored(DOCK_KEY, '42%'));
    const main = h('div', {class: 'oa-card', style: {minWidth: '0'}});
    grid.append(main);
    const view3d = this.view3d!;
    const modes = this.panels.get('modes') as ModesPanel | undefined;
    const embed = this.layout === 'full' && this.tab === 'modes';
    if (!embed) modes?.setEmbedded3D(null);

    if (this.layout === 'report') {
      for (const tab of TABS) {
        main.append(h('section', {class: 'oa-report-section'}, h('h2', null, tab.label), this.panel(tab.id).el));
        if (tab.id === 'overview')
          main.append(h('section', {class: 'oa-report-section'}, h('h2', null, '3D view'), view3d.el));
      }
    } else {
      const tabs = this.layout === 'full' ? [...TABS.slice(0, 1), {id: '3d', label: '3D view'}, ...TABS.slice(1)] : TABS;
      if (this.layout === 'docked' && this.tab === '3d') this.tab = 'overview';
      const strip = h('div', {class: 'oa-tabs', role: 'tablist'});
      const host = h('div');
      for (const tab of tabs) {
        const count = tab.id === 'interactions' ? ctx.features.length : tab.id === 'modes' ? ctx.modes.modes.length : null;
        strip.append(h('button', {class: `oa-tab ${tab.id === this.tab ? 'active' : ''}`, role: 'tab', dataset: {tab: tab.id},
          onclick: () => this.showTab(tab.id)}, tab.label, count !== null ? h('span', {class: 'count'}, String(count)) : null));
      }
      main.append(strip, host);
      host.append(this.tab === '3d' ? view3d.el : this.panel(this.tab).el);
      if (embed) (this.panel('modes') as ModesPanel).setEmbedded3D(view3d.el);
      if (this.layout === 'docked') {
        const splitter = h('div', {class: 'oa-splitter', title: 'Drag to resize'});
        this.wireSplitter(splitter, grid);
        grid.append(splitter, view3d.el);
      }
    }
    this.body.append(grid);
    this.fitHeight();
  }

  /** In the docked layout the results fill the window below the header, so the 3D view stays whole. */
  private fitHeight(): void {
    const grid = this.body.querySelector<HTMLElement>('.oa-results.docked');
    if (!grid) return;
    const top = grid.getBoundingClientRect().top + window.scrollY;
    grid.style.height = `${Math.max(520, window.innerHeight - top - 16)}px`;
  }

  private wireSplitter(splitter: HTMLElement, grid: HTMLElement): void {
    splitter.addEventListener('mousedown', (down) => {
      down.preventDefault();
      splitter.classList.add('dragging');
      const rect = grid.getBoundingClientRect();
      const move = (e: MouseEvent) => {
        const share = Math.min(70, Math.max(22, ((rect.right - e.clientX) / rect.width) * 100));
        grid.style.setProperty('--oa-dock-w', `${share.toFixed(1)}%`);
      };
      const up = () => {
        splitter.classList.remove('dragging');
        store(DOCK_KEY, grid.style.getPropertyValue('--oa-dock-w'));
        window.removeEventListener('mousemove', move);
        window.removeEventListener('mouseup', up);
      };
      window.addEventListener('mousemove', move);
      window.addEventListener('mouseup', up);
    });
  }

  showTab(id: string): void {
    if (this.layout === 'report') {
      this.panel(id).el.scrollIntoView({behavior: 'smooth', block: 'start'});
      return;
    }
    if (this.layout === 'docked' && id === '3d') id = 'overview';
    this.tab = id;
    this.arrange();
    this.onState({...this.params, tab: id});
  }

  private setLayout(layout: Layout): void {
    this.layout = layout;
    store(LAYOUT_KEY, layout);
    this.build();
  }

  private async reload(): Promise<void> {
    this.schemaCache.clear();
    await this.load().catch((e: Error) => toast(e.message, 'error'));
  }

  dispose(): void {
    this.disposed = true;
    window.removeEventListener('resize', this.onResize);
    this.teardown();
  }
}
