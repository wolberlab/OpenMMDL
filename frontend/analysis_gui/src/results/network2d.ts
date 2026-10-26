import {api} from '../api';
import {clear, h} from '../dom';
import {FeatureView, ResultsContext} from './context';

const cache = new Map<string, Promise<string>>();
let visScript: Promise<string> | null = null;

/** vis-network is inlined: a sandboxed srcdoc frame has a null origin, and browsers may refuse its script loads. */
function withInlineVis(html: string): Promise<string> {
  visScript ??= fetch('/static/vis-network.min.js').then((r) => r.text());
  return visScript.then((code) =>
    html.replace(/<script[^>]*src="[^"]*vis-network[^"]*"[^>]*><\/script>/, () => `<script>${code}</script>`));
}

/** ProLIF LigNetwork of a set of interactions, toggled in place through postMessage. */
export class LigandNetwork {
  readonly el = h('div', {class: 'oa-lignet'});
  private iframe: HTMLIFrameElement | null = null;
  private ready = false;
  private key = '';
  private active: string[] = [];
  private emphasis: string[] = [];
  private onMessage = (event: MessageEvent) => {
    if (!this.iframe || event.source !== this.iframe.contentWindow) return;
    const data = event.data ?? {};
    if (data.type === 'oa-network-ready') {
      this.ready = true;
      this.post();
    } else if (data.type === 'oa-network-toggle') {
      this.onToggle(data.features as string[]);
    } else if (data.type === 'oa-network-image') {
      this.imageWaiter?.(data.url as string | null, data.error as string | null);
      this.imageWaiter = null;
    }
  };
  private imageWaiter: ((url: string | null, error: string | null) => void) | null = null;

  /** PNG data URL of the drawn network (its canvas lives in the sandboxed frame). */
  image(): Promise<string> {
    const frame = this.iframe?.contentWindow;
    if (!this.ready || !frame) return Promise.reject(new Error('the network is not drawn yet'));
    return new Promise((resolve, reject) => {
      // The frame is busy while vis lays the network out; give it time before giving up.
      const timer = window.setTimeout(() => reject(new Error('the network did not answer')), 15000);
      this.imageWaiter = (url, error) => {
        clearTimeout(timer);
        if (url) resolve(url);
        else reject(new Error(error ?? 'the network is not drawn yet'));
      };
      frame.postMessage({type: 'oa-network-export'}, '*');
    });
  }

  constructor(private ctx: ResultsContext, private onToggle: (ids: string[]) => void) {
    window.addEventListener('message', this.onMessage);
  }

  show(features: FeatureView[], hidden: Set<string>, emphasis: string[] = []): void {
    this.active = features.filter((f) => !hidden.has(f.id)).map((f) => f.id);
    this.emphasis = emphasis;
    const ids = features.map((f) => f.id);
    const key = `${this.ctx.schemaDir}|${ids.join(',')}`;
    if (key === this.key) {
      this.post();
      return;
    }
    this.key = key;
    this.ready = false;
    clear(this.el);
    if (!ids.length || !this.ctx.summary.ligand) {
      this.el.append(h('div', {class: 'oa-empty'}, 'Nothing to draw for this mode.'));
      return;
    }
    this.el.append(h('div', {class: 'oa-empty', style: {border: '0'}}, h('div', {class: 'oa-spinner'}), 'Building the ProLIF network…'));
    if (!cache.has(key)) {
      const promise = api.lignetwork(this.ctx.root, this.ctx.schemaDir, this.ctx.summary.ligand, ids).then(withInlineVis);
      promise.catch(() => cache.delete(key));
      cache.set(key, promise);
    }
    cache.get(key)!.then((html) => {
      if (key !== this.key) return;
      clear(this.el);
      this.iframe = h('iframe', {sandbox: 'allow-scripts', title: 'ProLIF ligand network'});
      this.iframe.srcdoc = html;
      this.el.append(this.iframe,
        h('button', {class: 'oa-btn oa-btn-sm oa-lignet-fit', title: 'Fit to view',
          onclick: () => this.iframe?.contentWindow?.postMessage({type: 'oa-network-fit'}, '*')}, 'Fit'));
    }).catch((e: Error) => {
      if (key !== this.key) return;
      clear(this.el);
      this.el.append(h('div', {class: 'oa-empty'}, `Could not build the ProLIF network: ${e.message}`));
    });
  }

  private post(): void {
    if (!this.ready || !this.iframe?.contentWindow) return;
    this.iframe.contentWindow.postMessage({type: 'oa-network-state', active: this.active, emphasis: this.emphasis}, '*');
  }

  dispose(): void {
    window.removeEventListener('message', this.onMessage);
  }
}
