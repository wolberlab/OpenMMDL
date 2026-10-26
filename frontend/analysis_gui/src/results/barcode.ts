import {h} from '../dom';

/** Size a canvas to its CSS box at device resolution; returns the 2D context and size. */
export function fitCanvas(canvas: HTMLCanvasElement): {g: CanvasRenderingContext2D; w: number; h: number} | null {
  const dpr = window.devicePixelRatio || 1;
  const w = Math.round(canvas.clientWidth * dpr);
  const h = Math.round(canvas.clientHeight * dpr);
  if (!w || !h) return null;
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w;
    canvas.height = h;
  }
  const g = canvas.getContext('2d')!;
  g.clearRect(0, 0, w, h);
  return {g, w, h};
}

/** Draw presence runs of a 0/1 frame sequence, plus an optional frame cursor. */
export function drawBarcode(canvas: HTMLCanvasElement, flags: Uint8Array, color: string, cursor: number | null): void {
  const fit = fitCanvas(canvas);
  if (!fit) return;
  const {g, w, h} = fit;
  const n = flags.length;
  g.fillStyle = color;
  let start = -1;
  for (let i = 0; i <= n; i++) {
    if (i < n && flags[i]) {
      if (start < 0) start = i;
    } else if (start >= 0) {
      const x0 = Math.floor((start / n) * w);
      const x1 = Math.max(x0 + 1, Math.ceil((i / n) * w));
      g.fillRect(x0, 0, x1 - x0, h);
      start = -1;
    }
  }
  drawCursor(g, cursor, n, w, h);
}

/** Colour every frame position by a category index (mode ethogram). */
export function drawCategories(canvas: HTMLCanvasElement, labels: ArrayLike<number>, color: (label: number) => string,
  cursor: number | null): void {
  const fit = fitCanvas(canvas);
  if (!fit) return;
  const {g, w, h} = fit;
  const n = labels.length;
  let start = 0;
  for (let i = 1; i <= n; i++) {
    if (i === n || labels[i] !== labels[start]) {
      g.fillStyle = color(labels[start]);
      const x0 = Math.floor((start / n) * w);
      g.fillRect(x0, 0, Math.max(1, Math.ceil((i / n) * w) - x0), h);
      start = i;
    }
  }
  drawCursor(g, cursor, n, w, h);
}

function drawCursor(g: CanvasRenderingContext2D, cursor: number | null, n: number, w: number, h: number): void {
  if (cursor === null || cursor < 0) return;
  const x = ((cursor + 0.5) / n) * w;
  g.fillStyle = 'rgba(15, 23, 42, 0.85)';
  g.fillRect(Math.round(x) - 1, 0, Math.max(2, Math.round((window.devicePixelRatio || 1) * 1.5)), h);
}

/** Frame position under the mouse on a barcode canvas. */
export function positionAt(canvas: HTMLCanvasElement, event: MouseEvent, n: number): number {
  const rect = canvas.getBoundingClientRect();
  return Math.min(n - 1, Math.max(0, Math.floor(((event.clientX - rect.left) / rect.width) * n)));
}

/** Frame-number labels aligned with barcodes of the same width. */
export function frameAxis(frames: number[], ticks = 6): HTMLElement {
  const el = h('div', {class: 'oa-axis'});
  const n = frames.length;
  if (!n) return el;
  const seen = new Set<number>();
  for (let t = 0; t <= ticks; t++) {
    const p = Math.min(n - 1, Math.round((t / ticks) * (n - 1)));
    if (seen.has(p)) continue;
    seen.add(p);
    el.append(h('span', {style: {left: `${((p + 0.5) / n) * 100}%`}}, String(frames[p])));
  }
  return el;
}
