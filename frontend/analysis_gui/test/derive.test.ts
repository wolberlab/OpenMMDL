import assert from 'node:assert/strict';
import {test} from 'node:test';
import {decodeRuns, deriveModes, markov, passingFeatures} from '../src/derive.ts';

const flags = (bits: string) => Uint8Array.from(bits, (b) => Number(b));
const feature = (id: string, bits: string) => {
  const f = flags(bits);
  return {id, flags: f, count: f.reduce((a, b) => a + b, 0)};
};

test('decodeRuns expands run-length pairs', () => {
  assert.deepEqual([...decodeRuns([0, 2, 3, 1], 5)], [1, 1, 0, 1, 0]);
  assert.deepEqual([...decodeRuns([], 3)], [0, 0, 0]);
});

test('features pass when present in at least the threshold share of frames', () => {
  const features = [feature('a', '1100'), feature('b', '1000'), feature('c', '1111')];
  assert.deepEqual(passingFeatures(features, 4, 50).map((f) => f.id), ['a', 'c']);
  assert.deepEqual(passingFeatures(features, 4, 0).map((f) => f.id), ['a', 'b', 'c']);
});

test('modes are labelled by first appearance and ranked by frame count', () => {
  const features = [feature('a', '110011'), feature('b', '011110'), feature('rare', '100000')];
  const result = deriveModes(features, 6, 50);
  assert.deepEqual(result.passing, ['a', 'b']);
  assert.deepEqual([...result.labels], [0, 1, 2, 2, 1, 0]);
  const byLabel = Object.fromEntries(result.modes.map((m) => [m.label, m]));
  assert.deepEqual(byLabel.Binding_Mode_1.features, ['a']);
  assert.deepEqual(byLabel.Binding_Mode_2.features, ['a', 'b']);
  assert.deepEqual(byLabel.Binding_Mode_3.features, ['b']);
  assert.equal(result.modes[0].label, 'Binding_Mode_1');
  assert.equal(result.modes[0].rank, 0);
  assert.equal(byLabel.Binding_Mode_1.percent, (100 * 2) / 6);
});

test('the representative frame is closest to the mode average over all features', () => {
  const features = [feature('a', '1111'), feature('noise', '1000'), feature('extra', '0011')];
  const result = deriveModes(features, 4, 100);
  assert.equal(result.modes.length, 1);
  assert.equal(result.modes[0].representative, 1);
});

test('markov keeps transitions above the share of frames, with their reverse', () => {
  // 0 0 1 1 0 2: 0->1 once, 1->0 once, 0->2 once, self loops 0 and 1 once each.
  const labels = Int32Array.from([0, 0, 1, 1, 0, 2]);
  const all = markov(labels, 0);
  const edge = (a: number, b: number) => all.edges.find((e) => e.from === a && e.to === b);
  assert.equal(edge(0, 1)?.count, 1);
  assert.equal(edge(0, 2)?.probability, (100 * 1) / 3);
  assert.equal(edge(2, 0), undefined);
  assert.equal(all.nodes.find((n) => n.mode === 0)?.selfProbability, (100 * 1) / 3);

  const strict = markov(labels, 20);
  assert.equal(strict.edges.length, 0);
  assert.equal(strict.nodes.length, 0);
  assert.equal(strict.hiddenModes, 3);
});

test('markov colours top modes by the third of the trajectory they mostly occupy', () => {
  const labels = Int32Array.from([0, 0, 0, 1, 1, 1, 2, 2, 2]);
  const {nodes} = markov(labels, 0);
  const color = (m: number) => nodes.find((n) => n.mode === m)?.color;
  assert.equal(color(0), 'first');
  assert.equal(color(1), 'second');
  assert.equal(color(2), 'third');
});
