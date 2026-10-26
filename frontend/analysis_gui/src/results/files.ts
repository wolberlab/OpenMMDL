import {api, FileEntry} from '../api';
import {formatBytes, h, icon} from '../dom';
import {ResultsContext} from './context';

function lightbox(src: string, name: string): void {
  const box = h('div', {class: 'oa-lightbox', onclick: () => box.remove()},
    h('img', {src, alt: name}),
    h('a', {class: 'oa-btn oa-lightbox-download', href: src, download: name,
      onclick: (e: Event) => e.stopPropagation()}, icon('download'), `Download ${name}`));
  document.body.appendChild(box);
}

function zipLink(href: string, label: string): HTMLElement {
  return h('a', {class: 'oa-btn oa-btn-sm', href, download: '', onclick: (e: Event) => e.stopPropagation()},
    icon('download'), label);
}

function groupName(path: string): string {
  const i = path.indexOf('/');
  return i < 0 ? 'Run folder' : path.slice(0, i);
}

/** Everything the run wrote, with figures as a gallery. */
export class FilesPanel {
  readonly el: HTMLElement;

  constructor(ctx: ResultsContext) {
    const groups = new Map<string, FileEntry[]>();
    for (const file of ctx.summary.files) {
      const name = groupName(file.path);
      if (!groups.has(name)) groups.set(name, []);
      groups.get(name)!.push(file);
    }
    this.el = h('div', null,
      h('div', {class: 'oa-panel-head'}, h('div', null, h('h3', null, 'Files'),
        h('div', {class: 'oa-small oa-muted'}, `${ctx.summary.files.length} files in ${ctx.root}`)),
      h('div', {class: 'oa-panel-tools'},
        zipLink(api.zipUrl(ctx.root, '', true), 'All figures (.zip)'),
        zipLink(api.zipUrl(ctx.root), 'Everything (.zip)'))));
    const order = [...groups.keys()].sort((a, b) =>
      Number(b === 'Run folder') - Number(a === 'Run folder') || a.localeCompare(b));
    for (const name of order) {
      const files = groups.get(name)!;
      const images = files.filter((f) => f.image);
      const others = files.filter((f) => !f.image);
      const gallery = h('div', {class: 'oa-gallery'});
      for (const file of images) {
        const url = api.fileUrl(ctx.root, file.path);
        const label = file.path.slice(name === 'Run folder' ? 0 : name.length + 1);
        gallery.append(file.inline === false ?
          h('a', {class: 'oa-thumb', href: url, target: '_blank', title: 'Large image, opens in a new tab'},
            h('div', {class: 'oa-thumb-big'}, `Large image${file.width ? ` (${file.width}×${file.height} px)` : ''}, open in a new tab`),
            h('span', null, label)) :
          h('button', {class: 'oa-thumb', type: 'button', onclick: () => lightbox(url, file.path.split('/').pop()!)},
            h('img', {src: url, alt: label, loading: 'lazy'}), h('span', null, label)));
      }
      const table = h('table', {class: 'oa-file-table'});
      for (const file of others)
        table.append(h('tr', null,
          h('td', null, h('a', {href: api.fileUrl(ctx.root, file.path), target: '_blank'}, file.path)),
          h('td', null, formatBytes(file.size))));
      this.el.append(h('details', {class: 'oa-file-group oa-details', open: name === 'Run folder' || name === 'RMSD'},
        h('summary', null, `${name} (${files.length}) `,
          name === 'Run folder' ? null : zipLink(api.zipUrl(ctx.root, name), 'Folder (.zip)')),
        h('div', {style: {marginTop: '8px'}}, images.length ? gallery : null,
          others.length ? h('div', {style: {marginTop: images.length ? '10px' : '0'}}, table) : null)));
    }
  }

  dispose(): void {}
}
