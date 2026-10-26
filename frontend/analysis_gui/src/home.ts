import {h, icon} from './dom';

export interface RecentEntry {
  root: string;
  name: string;
  opened: number;
}

const RECENT_KEY = 'openmmdl-analysis-gui:recent';

export function recentResults(): RecentEntry[] {
  try {
    return JSON.parse(localStorage.getItem(RECENT_KEY) ?? '[]') as RecentEntry[];
  } catch {
    return [];
  }
}

export function rememberResults(root: string, name: string): void {
  const list = [{root, name, opened: Date.now()}, ...recentResults().filter((r) => r.root !== root)].slice(0, 8);
  try {
    localStorage.setItem(RECENT_KEY, JSON.stringify(list));
  } catch {
    // Recent folders are a convenience; ignore unavailable storage.
  }
}

export function homePage(actions: {newRun: () => void; openResults: () => void; open: (root: string) => void}): HTMLElement {
  const recent = recentResults();
  return h('div', {class: 'oa-card'},
    h('div', {class: 'oa-hero'},
      h('h1', null, 'Welcome to OpenMMDL Analysis'),
      h('p', null, 'Analyze protein-ligand interactions along an MD trajectory: run a new analysis with a live log, or open the results folder of a finished one to explore interactions, binding modes, their Markov chain, RMSD and waters next to the 3D structure.')),
    h('div', {class: 'oa-choice-grid'},
      h('button', {class: 'oa-choice', onclick: actions.newRun},
        h('div', {class: 'oa-choice-icon'}, icon('run')),
        h('div', null, h('div', {class: 'oa-choice-title'}, 'Run a new analysis'),
          h('div', {class: 'oa-choice-desc'}, 'Pick the topology, trajectory and ligand (or an OpenMMDL simulation folder), set the options and start openmmdl analysis.'))),
      h('button', {class: 'oa-choice', onclick: actions.openResults},
        h('div', {class: 'oa-choice-icon'}, icon('chart')),
        h('div', null, h('div', {class: 'oa-choice-title'}, 'Open analysis results'),
          h('div', {class: 'oa-choice-desc'}, 'Select a folder written by openmmdl analysis (it contains BindingModes_residue, RMSD, …) and explore it.')))),
    recent.length ? h('div', {class: 'oa-recent'}, h('h3', null, 'Recently opened'),
      h('div', {class: 'oa-recent-list'}, ...recent.map((r) => h('button', {class: 'oa-recent-item', onclick: () => actions.open(r.root)},
        h('b', null, r.name), h('span', {class: 'oa-recent-path'}, r.root))))) : null);
}
