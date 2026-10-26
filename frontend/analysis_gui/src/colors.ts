/** Interaction types as named by OpenMMDL's barcode grouping, with display colours.
 * Hues follow OpenMMDL's 2D binding-mode figures; UI shades are darkened for contrast on white. */
export interface InteractionType {key: string; label: string; color: string; halo: string}

export const TYPES: InteractionType[] = [
  {key: 'hydrophobic', label: 'Hydrophobic', color: '#c9a400', halo: 'rgb(255,255,0)'},
  {key: 'acceptor', label: 'H-bond acceptor', color: '#e0626b', halo: 'rgb(255,153,153)'},
  {key: 'donor', label: 'H-bond donor', color: '#3f6fe8', halo: 'rgb(77,128,255)'},
  {key: 'waterbridge', label: 'Water bridge', color: '#00a6b5', halo: 'rgb(0,255,230)'},
  {key: 'pistacking', label: 'π-stacking', color: '#2537c8', halo: 'rgb(0,0,255)'},
  {key: 'pication', label: 'π-cation', color: '#5b49c9', halo: 'rgb(0,0,255)'},
  {key: 'saltbridge_ni', label: 'Salt bridge (NI)', color: '#e08600', halo: 'rgb(255,153,0)'},
  {key: 'saltbridge_pi', label: 'Salt bridge (PI)', color: '#1f9e8f', halo: 'rgb(77,230,204)'},
  {key: 'halogen', label: 'Halogen bond', color: '#d02fb4', halo: 'rgb(255,0,230)'},
  {key: 'metal', label: 'Metal complex', color: '#b86b00', halo: 'rgb(255,153,0)'},
];

const BY_KEY = new Map(TYPES.map((t) => [t.key, t]));

export function typeInfo(key: string): InteractionType {
  return BY_KEY.get(key) ?? {key, label: key, color: '#64748b', halo: 'rgb(148,163,184)'};
}

/** Node colours of OpenMMDL's Markov chain figures. */
export const MARKOV_COLORS = {
  notTop10: {color: '#87ceeb', label: 'Not among the 10 most frequent modes'},
  first: {color: '#22a447', label: 'Mostly in the first third of the frames'},
  second: {color: '#f59e0b', label: 'Mostly in the second third'},
  third: {color: '#e5484d', label: 'Mostly in the last third'},
  even: {color: '#facc15', label: 'Spread over the whole trajectory'},
};

/** Qualitative palette for binding-mode identity (ethogram strip, mode markers). */
export const MODE_PALETTE = [
  '#4f46e5', '#0ea5e9', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#14b8a6', '#f97316',
  '#ec4899', '#84cc16',
];
export const MODE_OTHER = '#cbd5e1';
