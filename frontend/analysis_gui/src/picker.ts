import {api, DirListing} from './api';
import {clear, formatBytes, h, icon} from './dom';

export interface PickOptions {
  title: string;
  mode: 'file' | 'folder';
  ext?: string[];
  start?: string;
}

const LAST_DIR_KEY = 'openmmdl-analysis-gui:last-dir';

function lastDir(): string {
  try {
    return localStorage.getItem(LAST_DIR_KEY) ?? '';
  } catch {
    return '';
  }
}

function rememberDir(path: string): void {
  try {
    localStorage.setItem(LAST_DIR_KEY, path);
  } catch {
    // Storage can be unavailable (private window); the picker just starts at home then.
  }
}

function parentOf(path: string): string {
  const trimmed = path.replace(/[\\/]+$/, '');
  const i = Math.max(trimmed.lastIndexOf('/'), trimmed.lastIndexOf('\\'));
  return i > 0 ? trimmed.slice(0, i) : trimmed;
}

/** Browse the local file system through the GUI server. Resolves to a path or null. */
export function pickPath(options: PickOptions): Promise<string | null> {
  return new Promise((resolve) => {
    let listing: DirListing | null = null;
    let selected: string | null = null;
    const matches = (name: string) =>
      !options.ext?.length || options.ext.some((e) => name.toLowerCase().endsWith(e));

    const pathInput = h('input', {class: 'oa-input oa-input-mono', spellcheck: false,
      onkeydown: (e: KeyboardEvent) => { if (e.key === 'Enter') open(pathInput.value); }});
    const places = h('div', {class: 'oa-picker-places'});
    const list = h('div', {class: 'oa-picker-list'});
    const status = h('div', {class: 'oa-small oa-muted'});
    const selectBtn = h('button', {class: 'oa-btn oa-btn-primary', onclick: () => done(currentChoice())},
      options.mode === 'folder' ? 'Select this folder' : 'Select file');

    const close = h('button', {class: 'oa-btn oa-btn-icon', title: 'Close', onclick: () => done(null)}, icon('x'));
    const modal = h('div', {class: 'oa-modal', role: 'dialog', 'aria-label': options.title},
      h('div', {class: 'oa-modal-head'}, h('h3', null, options.title), close),
      h('div', {class: 'oa-modal-body'},
        places,
        h('div', {class: 'oa-picker-path'},
          h('button', {class: 'oa-btn oa-btn-icon', title: 'Parent folder',
            onclick: () => listing?.parent && open(listing.parent)}, icon('up')),
          pathInput,
          h('button', {class: 'oa-btn', onclick: () => open(pathInput.value)}, 'Go')),
        list),
      h('div', {class: 'oa-modal-foot'}, status,
        h('div', {style: {display: 'flex', gap: '8px'}},
          h('button', {class: 'oa-btn', onclick: () => done(null)}, 'Cancel'), selectBtn)));
    const back = h('div', {class: 'oa-modal-back', onclick: (e: MouseEvent) => { if (e.target === back) done(null); }}, modal);
    const onKey = (e: KeyboardEvent) => { if (e.key === 'Escape') done(null); };
    document.addEventListener('keydown', onKey);
    document.body.appendChild(back);

    function currentChoice(): string | null {
      if (options.mode === 'folder') return listing?.path ?? null;
      return selected;
    }

    function done(value: string | null) {
      document.removeEventListener('keydown', onKey);
      back.remove();
      if (value) rememberDir(options.mode === 'folder' ? value : parentOf(value));
      resolve(value);
    }

    function updateStatus() {
      const choice = currentChoice();
      selectBtn.disabled = !choice;
      status.textContent = choice ? choice : options.mode === 'file' ?
        `Pick a ${options.ext?.join(' / ') ?? ''} file` : '';
      if (options.mode === 'folder' && listing?.isResults)
        status.prepend(h('span', {class: 'oa-badge', style: {marginRight: '8px'}}, 'analysis results'));
    }

    async function open(path: string) {
      try {
        listing = await api.list(path);
      } catch (e) {
        status.textContent = (e as Error).message;
        return;
      }
      pathInput.value = listing.path;
      selected = null;
      clear(list);
      if (!listing.entries.length) list.append(h('div', {class: 'oa-picker-item dim'}, 'Empty folder'));
      for (const entry of listing.entries) {
        const usable = entry.dir || (options.mode === 'file' && matches(entry.name));
        const item = h('button', {class: `oa-picker-item ${usable ? '' : 'dim'}`, type: 'button'},
          icon(entry.dir ? 'folder' : 'file'),
          h('span', {class: 'name'}, entry.name),
          h('span', {class: 'size'}, entry.size === null ? '' : formatBytes(entry.size)));
        if (entry.dir) {
          item.addEventListener('click', () => open(entry.path));
        } else if (usable) {
          item.addEventListener('click', () => {
            list.querySelectorAll('.selected').forEach((el) => el.classList.remove('selected'));
            item.classList.add('selected');
            selected = entry.path;
            updateStatus();
          });
          item.addEventListener('dblclick', () => done(entry.path));
        }
        list.append(item);
      }
      updateStatus();
    }

    api.places().then((items) => {
      for (const place of items)
        places.append(h('button', {class: 'oa-btn oa-btn-sm', onclick: () => open(place.path)}, place.label));
    });
    const start = options.start ? (options.mode === 'file' ? parentOf(options.start) : options.start) : lastDir();
    open(start).then(() => {
      if (!listing && start) open('');
    });
  });
}
