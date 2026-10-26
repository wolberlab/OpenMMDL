type Child = Node | string | number | null | undefined | false | Child[];
type Attrs = Record<string, unknown> & {
  class?: string;
  style?: Partial<CSSStyleDeclaration>;
  dataset?: Record<string, string>;
};

/** Create an element: h('div', {class: 'x', onclick: fn}, child, ...). */
export function h<K extends keyof HTMLElementTagNameMap>(
  tag: K, attrs: Attrs | null = null, ...children: Child[]): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  if (attrs) {
    for (const [key, value] of Object.entries(attrs)) {
      if (value === undefined || value === null || value === false) continue;
      if (key === 'class') el.className = String(value);
      else if (key === 'style') Object.assign(el.style, value);
      else if (key === 'dataset') Object.assign(el.dataset, value);
      else if (key.startsWith('on') && typeof value === 'function')
        el.addEventListener(key.slice(2), value as EventListener);
      else if (value === true) el.setAttribute(key, '');
      else if (key in el && !key.includes('-')) (el as unknown as Record<string, unknown>)[key] = value;
      else el.setAttribute(key, String(value));
    }
  }
  append(el, children);
  return el;
}

export function append(parent: Node, children: Child[]): void {
  for (const child of children) {
    if (child === null || child === undefined || child === false) continue;
    if (Array.isArray(child)) append(parent, child);
    else parent.appendChild(child instanceof Node ? child : document.createTextNode(String(child)));
  }
}

export function clear(el: Element): void {
  while (el.firstChild) el.removeChild(el.firstChild);
}

const SVG_NS = 'http://www.w3.org/2000/svg';

export function s<K extends keyof SVGElementTagNameMap>(
  tag: K, attrs: Record<string, string | number | undefined> = {}, ...children: (SVGElement | string)[]): SVGElementTagNameMap[K] {
  const el = document.createElementNS(SVG_NS, tag);
  for (const [key, value] of Object.entries(attrs))
    if (value !== undefined) el.setAttribute(key, String(value));
  for (const child of children)
    el.appendChild(typeof child === 'string' ? document.createTextNode(child) : child);
  return el;
}

/** Inline SVG icon from a path list (24x24 viewBox, stroke style). */
export function icon(name: keyof typeof ICONS): SVGSVGElement {
  const svg = s('svg', {viewBox: '0 0 24 24', fill: 'none', stroke: 'currentColor',
    'stroke-width': 2, 'stroke-linecap': 'round', 'stroke-linejoin': 'round'});
  for (const d of ICONS[name]) svg.appendChild(s('path', {d}));
  return svg;
}

const ICONS = {
  folder: ['M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z'],
  file: ['M14 3H7a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h10a2 2 0 0 0 2-2V8z', 'M14 3v5h5'],
  play: ['M7 4l13 8-13 8z'],
  pause: ['M7 4h3v16H7z', 'M14 4h3v16h-3z'],
  up: ['M12 19V5', 'M5 12l7-7 7 7'],
  run: ['M5 3l14 9-14 9z'],
  chart: ['M4 19V5', 'M4 19h16', 'M8 15l3-4 3 3 5-7'],
  copy: ['M9 9h11v11H9z', 'M5 15H4V4h11v1'],
  reset: ['M3 12a9 9 0 1 0 3-6.7', 'M3 4v5h5'],
  x: ['M6 6l12 12', 'M18 6L6 18'],
  layers: ['M12 3l9 5-9 5-9-5z', 'M3 13l9 5 9-5'],
  download: ['M12 4v11', 'M7 10l5 5 5-5', 'M5 20h14'],
  search: ['M11 18a7 7 0 1 0 0-14 7 7 0 0 0 0 14z', 'M20 20l-4-4'],
  prev: ['M15 6l-6 6 6 6'],
  next: ['M9 6l6 6-6 6'],
};

let tooltipEl: HTMLDivElement | null = null;

export function showTooltip(event: MouseEvent, content: Node | string): void {
  if (!tooltipEl) {
    tooltipEl = h('div', {class: 'oa-tooltip'});
    document.body.appendChild(tooltipEl);
  }
  clear(tooltipEl);
  append(tooltipEl, [content]);
  tooltipEl.style.display = 'block';
  const pad = 14;
  const rect = tooltipEl.getBoundingClientRect();
  let x = event.clientX + pad;
  let y = event.clientY + pad;
  if (x + rect.width > window.innerWidth - 8) x = event.clientX - rect.width - pad;
  if (y + rect.height > window.innerHeight - 8) y = event.clientY - rect.height - pad;
  tooltipEl.style.left = `${Math.max(8, x)}px`;
  tooltipEl.style.top = `${Math.max(8, y)}px`;
}

export function hideTooltip(): void {
  if (tooltipEl) tooltipEl.style.display = 'none';
}

export function toast(message: string, kind: 'info' | 'error' = 'info', ms = 4000): void {
  const el = h('div', {class: `oa-toast ${kind === 'error' ? 'error' : ''}`}, message);
  document.body.appendChild(el);
  setTimeout(() => el.remove(), ms);
}

export function formatBytes(n: number): string {
  if (n < 1024) return `${n} B`;
  const units = ['KB', 'MB', 'GB', 'TB'];
  let value = n / 1024;
  let i = 0;
  while (value >= 1024 && i < units.length - 1) {
    value /= 1024;
    i++;
  }
  return `${value.toFixed(value < 10 ? 1 : 0)} ${units[i]}`;
}

export function pct(value: number, digits = 1): string {
  return `${value.toFixed(digits)}%`;
}
