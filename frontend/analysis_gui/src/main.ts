import './styles.css';
import {clear, h} from './dom';
import {homePage, rememberResults} from './home';
import {pickPath} from './picker';
import {ResultsViewer, ViewerParams} from './results/viewer';
import {runFormPage} from './run-form';
import {runMonitorPage} from './run-monitor';

declare global {
  interface Window {
    OPENMMDL_GUI: {version: string; logo: string; initial: {results?: string}};
  }
}

const DOCS_URL = 'https://openmmdl.readthedocs.io/en/latest/openmmdl_analysis.html';

interface Page {
  el: HTMLElement;
  dispose(): void;
}

const config = window.OPENMMDL_GUI;
const page = h('div', {class: 'oa-page'});
const navLinks = new Map<string, HTMLButtonElement>();
let current: Page | null = null;
let rerunValues: Record<string, unknown> | null = null;

function nav(id: string, label: string, go: () => void): HTMLButtonElement {
  const btn = h('button', {class: 'oa-navlink', onclick: go}, label);
  navLinks.set(id, btn);
  return btn;
}

async function openResultsPicker(): Promise<void> {
  const root = await pickPath({title: 'Open an analysis results folder', mode: 'folder'});
  if (root) location.hash = resultsHash({root});
}

function resultsHash(params: ViewerParams): string {
  const q = new URLSearchParams({root: params.root});
  if (params.schema) q.set('schema', params.schema);
  if (params.tab) q.set('tab', params.tab);
  return `#/results?${q.toString()}`;
}

const header = h('header', {class: 'oa-header'},
  h('nav', {class: 'oa-header-nav'},
    nav('home', 'Home', () => (location.hash = '#/')),
    nav('run', 'New analysis', () => (location.hash = '#/run')),
    nav('results', 'Open results', openResultsPicker)),
  h('div', {class: 'oa-header-center'}, h('img', {class: 'oa-logo', src: config.logo, alt: 'OpenMMDL Analysis'})),
  h('div', {class: 'oa-header-actions'},
    h('a', {class: 'oa-btn oa-btn-sm', href: DOCS_URL, target: '_blank', rel: 'noopener'}, 'Documentation')));

const title = h('div', {class: 'oa-pagehead'});

function setTitle(text: string, sub = '', visible = true): void {
  title.replaceChildren(h('div', {class: 'oa-title'}, text), sub ? h('div', {class: 'oa-subtitle'}, sub) : '');
  title.style.display = visible ? '' : 'none';
  document.title = `OpenMMDL Analysis: ${text}`;
}

function show(next: Page, active: string): void {
  current?.dispose();
  current = next;
  clear(page);
  page.append(title, next.el);
  navLinks.forEach((btn, id) => btn.classList.toggle('active', id === active));
  window.scrollTo(0, 0);
}

function route(): void {
  const hash = location.hash.replace(/^#/, '') || '/';
  const [path, queryString] = hash.split('?');
  const query = new URLSearchParams(queryString ?? '');

  if (path === '/run') {
    setTitle('New analysis', 'Set up openmmdl analysis: inputs, ligand, options and the output folder.');
    const prefill = rerunValues;
    rerunValues = null;
    show(runFormPage((run) => (location.hash = `#/run/${run.id}`), prefill), 'run');
    return;
  }
  const runMatch = path.match(/^\/run\/([\w-]+)$/);
  if (runMatch) {
    setTitle('Running analysis');
    show(runMonitorPage(runMatch[1], {
      openResults: (root) => (location.hash = resultsHash({root})),
      edit: () => (location.hash = '#/run'),
    }), 'run');
    return;
  }
  if (path === '/results' && query.get('root')) {
    const params: ViewerParams = {root: query.get('root')!, schema: query.get('schema') ?? undefined, tab: query.get('tab') ?? undefined};
    const name = params.root.split(/[\\/]/).filter(Boolean).pop() ?? params.root;
    setTitle(name, '', false);
    rememberResults(params.root, name);
    const viewer = new ResultsViewer(params,
      (values) => {
        rerunValues = values;
        location.hash = '#/run';
      },
      (state) => history.replaceState(null, '', resultsHash(state)));
    show(viewer, 'results');
    return;
  }
  setTitle('Welcome');
  show({el: homePage({
    newRun: () => (location.hash = '#/run'),
    openResults: openResultsPicker,
    open: (root) => (location.hash = resultsHash({root})),
  }), dispose: () => {}}, 'home');
}

document.getElementById('app')!.append(h('div', {class: 'oa-shell'}, header, page));
if (!location.hash && config.initial?.results) location.hash = resultsHash({root: config.initial.results});
window.addEventListener('hashchange', route);
route();
