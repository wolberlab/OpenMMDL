import {h, icon, toast} from './dom';

const SVG_NS = 'http://www.w3.org/2000/svg';
// The page stylesheet does not travel with a downloaded SVG, so these are copied inline.
const STYLE_PROPS = ['fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-opacity', 'stroke-dasharray',
  'stroke-linecap', 'opacity', 'font-family', 'font-size', 'font-weight', 'paint-order', 'display', 'visibility'];

export function download(blob: Blob, name: string): void {
  const url = URL.createObjectURL(blob);
  const a = h('a', {href: url, download: name});
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 10000);
}

export function safeName(name: string): string {
  return name.replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '') || 'figure';
}

/** Self-contained SVG text of an on-screen SVG: computed styles inlined, white background added. */
export function standaloneSvg(svg: SVGSVGElement): {text: string; width: number; height: number} {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  const source = [svg, ...svg.querySelectorAll('*')];
  const target = [clone, ...clone.querySelectorAll('*')];
  source.forEach((el, i) => {
    // Embedded depictions (nested <svg>) carry their own style attributes.
    if (el !== svg && el.parentElement?.closest('svg') !== svg) return;
    const style = getComputedStyle(el);
    const out = (target[i] as SVGElement).style;
    for (const prop of STYLE_PROPS) out.setProperty(prop, style.getPropertyValue(prop));
  });
  const box = svg.viewBox.baseVal;
  const width = box && box.width ? box.width : svg.clientWidth;
  const height = box && box.height ? box.height : svg.clientHeight;
  clone.setAttribute('xmlns', SVG_NS);
  clone.setAttribute('width', String(width));
  clone.setAttribute('height', String(height));
  clone.removeAttribute('class');
  const background = document.createElementNS(SVG_NS, 'rect');
  for (const [k, v] of Object.entries({x: box?.x ?? 0, y: box?.y ?? 0, width, height, fill: '#ffffff'}))
    background.setAttribute(k, String(v));
  clone.insertBefore(background, clone.firstChild);
  return {text: new XMLSerializer().serializeToString(clone), width, height};
}

export async function svgToPng(svg: SVGSVGElement, scale = 3): Promise<Blob> {
  const {text, width, height} = standaloneSvg(svg);
  const img = new Image();
  img.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(text)}`;
  await img.decode();
  return canvasToPng(width * scale, height * scale, (g) => g.drawImage(img, 0, 0, width * scale, height * scale));
}

export function canvasToPng(width: number, height: number, draw: (g: CanvasRenderingContext2D) => void): Promise<Blob> {
  const canvas = h('canvas', {width: Math.round(width), height: Math.round(height)});
  const g = canvas.getContext('2d')!;
  g.fillStyle = '#ffffff';
  g.fillRect(0, 0, canvas.width, canvas.height);
  draw(g);
  return new Promise((resolve, reject) =>
    canvas.toBlob((blob) => (blob ? resolve(blob) : reject(new Error('Could not render the image'))), 'image/png'));
}

export function downloadCsv(rows: (string | number | null)[][], name: string): void {
  const cell = (v: string | number | null) => {
    const text = v === null ? '' : String(v);
    return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
  };
  download(new Blob([rows.map((r) => r.map(cell).join(',')).join('\n') + '\n'], {type: 'text/csv'}), `${safeName(name)}.csv`);
}

/** Small "SVG · PNG" (or custom) download buttons. */
export function exportButtons(actions: Record<string, () => Promise<void> | void>): HTMLElement {
  const wrap = h('div', {class: 'oa-export', title: 'Download'}, icon('download'));
  for (const [label, run] of Object.entries(actions))
    wrap.append(h('button', {class: 'oa-btn oa-btn-sm', type: 'button', onclick: async () => {
      try {
        await run();
      } catch (e) {
        toast(`Export failed: ${(e as Error).message}`, 'error');
      }
    }}, label));
  return wrap;
}

/** SVG and PNG downloads of the SVG returned by getSvg (null when nothing is drawn yet). */
export function svgExportButtons(getSvg: () => SVGSVGElement | null, name: () => string): HTMLElement {
  const svgOrFail = () => {
    const svg = getSvg();
    if (!svg) throw new Error('nothing to export yet');
    return svg;
  };
  return exportButtons({
    SVG: () => download(new Blob([standaloneSvg(svgOrFail()).text], {type: 'image/svg+xml'}), `${safeName(name())}.svg`),
    PNG: async () => download(await svgToPng(svgOrFail()), `${safeName(name())}.png`),
  });
}
