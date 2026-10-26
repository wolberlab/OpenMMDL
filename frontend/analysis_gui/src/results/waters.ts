import {api} from '../api';
import {h} from '../dom';
import {ResultsContext} from './context';
import {InteractionsPanel} from './interactions';

/** Water bridges from the interaction analysis and the stable-water clusters (-w). */
export class WatersPanel {
  readonly el: HTMLElement;
  private bridges: InteractionsPanel | null = null;

  constructor(ctx: ResultsContext, showExtra: (path: string | null) => void) {
    const hasBridges = ctx.features.some((f) => f.type === 'waterbridge');
    if (hasBridges) this.bridges = new InteractionsPanel(ctx, (type) => type === 'waterbridge');
    const stable = ctx.summary.stableWaters;

    const stableSection = stable.length ? h('div', null, ...stable.map((run) => h('div', {style: {marginBottom: '14px'}},
      h('div', {class: 'oa-small oa-muted', style: {marginBottom: '6px'}},
        `Cluster size ${run.eps.length > 1 ? `${run.eps[0]}.${run.eps.slice(1)}` : run.eps} Å · ${run.dir}`,
        run.csv ? h('span', null, ' · ', h('a', {href: api.fileUrl(ctx.root, run.csv), target: '_blank'}, 'stable_waters.csv')) : null),
      h('table', {class: 'oa-kv', style: {maxWidth: '760px'}},
        h('tr', null, ...['Present in at least', 'Clusters', 'Files', ''].map((t) => h('td', {style: {fontWeight: '700'}}, t))),
        ...run.levels.map((level) => h('tr', null,
          h('td', null, `${level.minFrames} frames`),
          h('td', null, String(level.clusters)),
          h('td', null,
            level.representative ? h('a', {href: api.fileUrl(ctx.root, level.representative), target: '_blank'}, 'representative_waters.pdb') : '–',
            level.residues ? h('span', null, ' · ', h('a', {href: api.fileUrl(ctx.root, level.residues), target: '_blank'}, 'interacting_residues.csv')) : null),
          h('td', null, level.representative ? h('button', {class: 'oa-btn oa-btn-sm',
            onclick: () => showExtra(level.representative)}, 'Show in 3D') : null)))),
      h('button', {class: 'oa-btn oa-btn-sm', style: {marginTop: '6px'}, onclick: () => showExtra(null)}, 'Hide stable waters')))) :
      h('div', {class: 'oa-empty'}, 'No stable water analysis in this folder. Run the analysis with “Stable water analysis” (-w) to cluster conserved water sites.');

    this.el = h('div', null,
      this.bridges ? this.bridges.el : h('div', {class: 'oa-empty'}, 'No water bridges were found.'),
      ctx.summary.interactingWaters.length ? h('div', {class: 'oa-small oa-muted', style: {margin: '8px 0 0'}},
        `${ctx.summary.interactingWaters.length} water molecules form bridges; the 3D view shows them (toggle “Waters”).`) : null,
      h('div', {style: {marginTop: '22px'}}, h('h3', {style: {margin: '0 0 10px', fontSize: '16px'}}, 'Stable waters'), stableSection));
  }

  dispose(): void {
    this.bridges?.dispose();
  }
}
