/** The steps of `openmmdl analysis` and how its log lines complete them. */

export interface Step {
  id?: 'interactions';
  label: string;
  /** Log line that marks the step as done; steps without one finish when a later step does. */
  done?: RegExp;
  /** Rough share of the run time, for the overall progress bar. */
  weight: number;
  group?: string;
}

export const GROUP_TITLES: Record<string, string> = {
  residue: 'Binding modes by residue',
  ligand: 'Binding modes by ligand atom',
};

/** Steps in log order for the options of a run. */
export function analysisSteps(values: Record<string, unknown>): Step[] {
  const peptide = Boolean(values.peptide);
  const steps: Step[] = [{label: 'Prepare input files', done: /Files are preprocessed/, weight: 2}];
  if (!peptide) steps.push({label: 'Ligand ring information', done: /Ligand ring data gathered/, weight: 1});
  steps.push({label: values.frame_rmsd ? 'RMSD over time and between frames' : 'RMSD over time',
    done: /will be analyzed/, weight: values.frame_rmsd ? 6 : 2});
  steps.push({id: 'interactions', label: 'Interactions in every frame', done: /trajectory processed/, weight: 55});
  for (const group of ['residue', 'ligand']) {
    const tag = `${group} bindig mode`;
    steps.push(
      {group, label: 'Interaction partners', done: /interaction partners generated/, weight: 2},
      {group, label: 'Markov chain', done: new RegExp(`${tag}: markov state figure generated`), weight: 2},
      {group, label: peptide ? 'Peptide figures' : 'Binding mode figures', done: /PNG image saved as|Binding mode figure generated/, weight: 3});
    if (values.representative_frame)
      steps.push({group, label: 'Representative frames', done: new RegExp(`${tag}: Binding mode pdb files saved`), weight: 12});
    steps.push({group, label: 'Barcodes', done: new RegExp(`${tag}: Barcodes generated`), weight: 4});
    if (values.generate_pml)
      steps.push({group, label: 'Pharmacophores', done: new RegExp(`${tag}: Pharmacophores generated`), weight: 2});
  }
  steps.push({label: 'Interaction analysis complete', done: /Analysis is Finished/, weight: 1});
  if (values.stable_water_analysis) steps.push({label: 'Stable water analysis', weight: 10});
  return steps;
}

/** Index of the log line that completes each step (null while open). Repeated messages are
 * matched in order, so the second "interaction partners generated" belongs to the second pass. */
export function matchSteps(steps: Step[], lines: string[]): (number | null)[] {
  const at: (number | null)[] = steps.map(() => null);
  let cursor = 0;
  steps.forEach((step, s) => {
    if (!step.done) return;
    for (let i = cursor; i < lines.length; i++) {
      if (step.done.test(lines[i])) {
        at[s] = i;
        cursor = i + 1;
        return;
      }
    }
  });
  return at;
}
