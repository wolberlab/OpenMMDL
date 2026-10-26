import assert from 'node:assert/strict';
import {test} from 'node:test';
import {analysisSteps, matchSteps} from '../src/run-steps.ts';

const LOG = [
  'INFO: Files are preprocessed',
  'INFO: Ligand ring data gathered',
  'INFO: RMSD calculated',
  'INFO: The whole trajectory consisting of 50 frames will be analyzed',
  'INFO: Processing protein-ligand trajectory',
  'INFO: Protein-ligand trajectory processed',
  'INFO: Running residue bindig mode analysis',
  'INFO: interaction partners generated',
  'INFO: residue bindig mode: markov state figure generated',
  'INFO: PNG image saved as: ligand_numbering.png',
  'INFO: residue bindig mode: found representative frame for each binding mode',
  'INFO:  residue bindig mode: Binding mode pdb files saved',
  'INFO: residue bindig mode: Barcodes generated',
  'INFO: Running ligand bindig mode analysis',
  'INFO: interaction partners generated',
];

test('steps follow the options of the run', () => {
  const plain = analysisSteps({}).map((s) => s.label);
  assert.ok(!plain.includes('Representative frames'));
  assert.ok(!plain.includes('Stable water analysis'));
  const full = analysisSteps({representative_frame: true, generate_pml: true, stable_water_analysis: true, frame_rmsd: true});
  assert.equal(full.filter((s) => s.label === 'Representative frames').length, 2);
  assert.equal(full.at(-1)!.label, 'Stable water analysis');
  assert.ok(!analysisSteps({peptide: 'B'}).some((s) => s.label === 'Ligand ring information'));
});

test('repeated log messages complete the steps of each pass in order', () => {
  const steps = analysisSteps({representative_frame: true});
  const at = matchSteps(steps, LOG);
  const lineOf = (label: string, group?: string) => at[steps.findIndex((s) => s.label === label && s.group === group)];
  assert.equal(lineOf('Interactions in every frame'), 5);
  assert.equal(lineOf('Interaction partners', 'residue'), 7);
  assert.equal(lineOf('Interaction partners', 'ligand'), 14);
  assert.equal(lineOf('Representative frames', 'residue'), 11);
  assert.equal(lineOf('Barcodes', 'residue'), 12);
  assert.equal(lineOf('Markov chain', 'ligand'), null);
});
