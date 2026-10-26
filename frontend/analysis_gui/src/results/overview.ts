import {typeInfo} from '../colors';
import {h, icon, pct, toast} from '../dom';
import {ResultsContext} from './context';
import {seriesLabel, stats} from './rmsd';

function tile(value: string, label: string, title?: string): HTMLElement {
  return h('div', {class: 'oa-tile', title}, h('div', {class: 'oa-tile-value'}, value), h('div', {class: 'oa-tile-label'}, label));
}

function formatValue(value: unknown): string {
  if (value === true) return 'yes';
  if (value === false || value === null || value === undefined || value === '') return '–';
  return String(value);
}

/** Summary of a results folder: key numbers, most frequent contacts, run settings, log messages. */
export class OverviewPanel {
  readonly el: HTMLElement;
  private tiles = h('div', {class: 'oa-tiles'});
  private disposers: (() => void)[] = [];

  constructor(private ctx: ResultsContext, private openTab: (tab: string) => void, private rerun: (() => void) | null) {
    const summary = ctx.summary;
    const top = [...ctx.features].sort((a, b) => b.count - a.count).slice(0, 12);
    const bars = h('div', {class: 'oa-hbars'});
    for (const f of top) {
      const share = (100 * f.count) / ctx.nFrames;
      const type = typeInfo(f.type);
      bars.append(h('div', {class: 'oa-hbar', title: `${f.label} · ${f.sublabel}`,
        onclick: () => { ctx.selectFeatures([f.id]); this.openTab('interactions'); }},
      h('span', {class: 'oa-hbar-label'}, f.label),
      h('div', {class: 'oa-hbar-track'}, h('div', {style: {width: `${share}%`, background: type.color}})),
      h('span', {class: 'oa-hbar-val'}, pct(share, 0))));
    }

    const values = summary.run?.values ?? {};
    const settings = h('table', {class: 'oa-kv'});
    const rows: [string, unknown][] = summary.run ? [
      ['Topology', values.topology], ['Trajectory', values.trajectory], ['Ligand', values.ligand_name ?? summary.ligand],
      ['Ligand SDF', values.ligand_sdf], ['Peptide chain', values.peptide], ['Engine', values.interaction_package ?? 'plip'],
      ['Mode threshold', `${summary.threshold}%`], ['Min transition', `${summary.minTransition}%`],
      ['Last frame', values.frames], ['Status', summary.run.status],
    ] : [
      ['Ligand', summary.ligand], ['Mode threshold', `${summary.threshold}% (CLI default; run file not found)`],
      ['Min transition', `${summary.minTransition}%`],
    ];
    for (const [key, value] of rows)
      settings.append(h('tr', null, h('td', null, key), h('td', {class: 'oa-mono', style: {overflowWrap: 'anywhere'}}, formatValue(value))));

    const messages = summary.log?.messages ?? [];
    this.el = h('div', null,
      this.tiles,
      h('div', {class: 'oa-two'},
        h('div', null, h('div', {class: 'oa-h4'}, 'Most frequent interactions'),
          top.length ? bars : h('div', {class: 'oa-empty'}, 'No interactions were found.')),
        h('div', null,
          h('div', {class: 'oa-h4'}, 'Analysis settings'), settings,
          summary.run?.command ? h('div', {style: {marginTop: '10px'}},
            h('div', {class: 'oa-command'}, summary.run.command),
            h('div', {style: {display: 'flex', gap: '6px', marginTop: '8px'}},
              h('button', {class: 'oa-btn oa-btn-sm', onclick: () => navigator.clipboard.writeText(summary.run!.command)
                .then(() => toast('Command copied'))}, icon('copy'), 'Copy command'),
              this.rerun ? h('button', {class: 'oa-btn oa-btn-sm', onclick: this.rerun}, icon('run'), 'Run again with changes') : null)) : null,
          messages.length ? h('div', {style: {marginTop: '16px'}},
            h('div', {class: 'oa-h4'}, `Log messages (${messages.length})`),
            h('div', {class: 'oa-notes'}, ...messages.slice(0, 12).map((m) =>
              h('div', {class: `oa-note ${m.level === 'warning' ? 'warn' : 'error'}`}, m.text)))) : null)));
    this.fillTiles();
    this.disposers.push(ctx.on('modes', () => this.fillTiles()));
  }

  private fillTiles(): void {
    const {ctx} = this;
    this.tiles.replaceChildren();
    const top = ctx.modes.modes[0];
    this.tiles.append(
      tile(String(ctx.nFrames), 'Analysed frames'),
      tile(String(ctx.features.length), 'Distinct interactions'),
      tile(String(ctx.modes.passing.length), `Above ${ctx.threshold}% threshold`),
      tile(String(ctx.modes.modes.length), 'Binding modes'),
      top ? tile(pct(top.percent), `Top mode (${top.label})`) : tile('–', 'Top mode'));
    const ligand = ctx.summary.rmsd?.series.find((sr) => sr.name.startsWith('resname ') || sr.name.startsWith('chainID '));
    if (ligand) this.tiles.append(tile(`${stats(ligand.values).median.toFixed(2)} Å`, `Median ${seriesLabel(ligand.name).toLowerCase()} RMSD`));
    if (ctx.summary.interactingWaters.length)
      this.tiles.append(tile(String(ctx.summary.interactingWaters.length), 'Bridging waters'));
  }

  dispose(): void {
    this.disposers.forEach((d) => d());
  }
}
