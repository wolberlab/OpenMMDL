import {api, RunSnapshot} from './api';
import {clear, h, icon, toast} from './dom';
import {analysisSteps, GROUP_TITLES, matchSteps} from './run-steps';

const BAR_RE = /^(.*?):?\s*(\d+)%\|[^|]*\|\s*(\d+)\/(\d+)\s*\[([^<\]]*)<([^,\]]*)(?:,\s*([^\]]*))?\]/;
const WARNING_START = /\.py:\d+: \w*Warning: /;
const WARNING_END = /^\s+warnings\.warn\(/;

interface Bar {label: string; percent: number; done: number; total: number; remaining: string; rate: string}

function parseBar(text: string): Bar | null {
  const m = text.match(BAR_RE);
  if (!m) return null;
  return {label: m[1].trim(), percent: Number(m[2]), done: Number(m[3]), total: Number(m[4]),
    remaining: m[6].trim(), rate: (m[7] ?? '').trim()};
}

function clock(seconds: number): string {
  const s = Math.max(0, Math.floor(seconds));
  const m = Math.floor(s / 60) % 60;
  const hr = Math.floor(s / 3600);
  return hr ? `${hr}:${String(m).padStart(2, '0')}:${String(s % 60).padStart(2, '0')}` : `${m}:${String(s % 60).padStart(2, '0')}`;
}

function roughly(seconds: number): string {
  if (seconds < 90) return 'under 2 min';
  if (seconds < 3600) return `about ${Math.round(seconds / 60)} min`;
  return `about ${(seconds / 3600).toFixed(1)} h`;
}

export function runMonitorPage(runId: string, actions: {openResults: (root: string) => void; edit: () => void}): {el: HTMLElement; dispose(): void} {
  const status = h('span', {class: 'oa-status running'}, 'running');
  const time = h('span', {class: 'oa-muted oa-small'});
  const overallBar = h('div', {style: {width: '0%'}});
  const overallPct = h('span', {class: 'oa-overall-pct'}, '0%');
  const overallEta = h('span', {class: 'oa-small oa-muted'});
  const now = h('div', {class: 'oa-now'});
  const timeline = h('div', {class: 'oa-timeline'});
  const log = h('div', {class: 'oa-log'});
  const command = h('div', {class: 'oa-command', style: {marginBottom: '10px'}});
  const out = h('div', {class: 'oa-small oa-muted', style: {overflowWrap: 'anywhere'}});
  const warningsBox = h('input', {type: 'checkbox', onchange: () => renderLog()});
  const warningsLabel = h('span');
  const cancelBtn = h('button', {class: 'oa-btn oa-btn-danger', onclick: async () => {
    cancelBtn.disabled = true;
    await api.cancelRun(runId).catch((e: Error) => toast(e.message, 'error'));
  }}, icon('x'), 'Cancel run');
  const openBtn = h('button', {class: 'oa-btn oa-btn-primary', disabled: true}, icon('chart'), 'Open results');
  const editBtn = h('button', {class: 'oa-btn', onclick: actions.edit}, 'Back to settings');

  const el = h('div', {class: 'oa-monitor'},
    h('div', {class: 'oa-card'},
      h('div', {style: {display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: '8px'}}, status, time),
      h('div', {class: 'oa-overall'},
        h('div', {class: 'oa-overall-head'}, overallPct, overallEta),
        h('div', {class: 'oa-progress oa-progress-lg'}, overallBar)),
      now, timeline,
      h('div', {style: {display: 'flex', flexDirection: 'column', gap: '8px', marginTop: '16px'}}, openBtn, cancelBtn, editBtn)),
    h('div', {class: 'oa-card'}, command, out,
      h('div', {class: 'oa-log-head'}, h('span', {class: 'oa-h4', style: {margin: '0'}}, 'Log'),
        h('label', {class: 'oa-check oa-small'}, warningsBox, warningsLabel)),
      log));

  let since = 0;
  let timer = 0;
  let disposed = false;
  let follow = true;
  let latest: RunSnapshot | null = null;
  let timed = true;
  const lines: string[] = [];
  const times: number[] = [];
  log.addEventListener('scroll', () => {
    follow = log.scrollTop + log.clientHeight >= log.scrollHeight - 20;
  });

  /** Python warnings (header line through its `warnings.warn(` line) are hidden unless asked for. */
  function renderLog(): void {
    const show = warningsBox.checked;
    clear(log);
    let hidden = 0;
    let inWarning = 0;
    for (const raw of lines) {
      const bar = parseBar(raw);
      const line = bar ? raw.replace(BAR_RE, '').trim() : raw;
      if (!line) continue;
      if (WARNING_START.test(line)) inWarning = 6;
      if (inWarning > 0) {
        inWarning = WARNING_END.test(line) ? 0 : inWarning - 1;
        hidden++;
        if (!show) continue;
        log.append(h('div', {class: 'warn'}, line));
        continue;
      }
      const cls = /^(ERROR|CRITICAL)|Traceback|Error:|Exception/.test(line) ? 'err' : /^WARNING/.test(line) ? 'warn' : 'info';
      log.append(h('div', {class: cls}, line));
    }
    warningsLabel.textContent = hidden ? `Show Python warnings (${hidden})` : 'Show Python warnings';
    if (latest?.status === 'failed')
      log.append(h('div', {class: 'err'}, `\nThe analysis stopped with exit code ${latest.returncode}. See the messages above.`));
    if (follow) log.scrollTop = log.scrollHeight;
  }

  function renderSteps(run: RunSnapshot): void {
    const steps = analysisSteps(run.values ?? {});
    const at = matchSteps(steps, lines);
    const finished = run.status === 'finished';
    const lastMatched = at.reduce<number>((acc, v, i) => (v !== null ? i : acc), -1);
    const done = steps.map((_, s) => finished || at[s] !== null || s < lastMatched);
    const current = run.status === 'running' ? done.indexOf(false) : -1;
    const failed = run.status === 'failed' || run.status === 'cancelled' ? done.indexOf(false) : -1;

    // When a step ended: its own line, or the next matched line for steps that print nothing.
    const endTime = steps.map((_, s) => {
      for (let t = s; t < steps.length; t++) if (at[t] !== null) return times[at[t]!] ?? null;
      return finished ? run.elapsed : null;
    });
    const duration = (s: number): number | null => {
      if (!timed) return null;
      const end = done[s] ? endTime[s] : s === current ? run.elapsed : null;
      if (end === null || end === undefined) return null;
      const start = s > 0 ? endTime[s - 1] ?? 0 : 0;
      return Math.max(0, end - start);
    };

    const bar = run.status === 'running' ? parseBar(run.progress) : null;
    const total = steps.reduce((sum, st) => sum + st.weight, 0);
    let fraction = steps.reduce((sum, st, s) => sum + (done[s] ? st.weight : 0), 0) / total;
    if (current >= 0 && bar) fraction += (steps[current].weight * bar.percent) / 100 / total;
    if (finished) fraction = 1;
    overallBar.style.width = `${(100 * fraction).toFixed(1)}%`;
    overallPct.textContent = `${Math.floor(100 * fraction)}%`;
    overallBar.parentElement!.classList.toggle('failed', failed >= 0);
    overallEta.textContent = finished ? `done in ${clock(run.elapsed)}` :
      run.status !== 'running' ? run.status :
      fraction > 0.03 ? `${roughly((run.elapsed * (1 - fraction)) / fraction)} left` : 'estimating…';

    clear(now);
    if (current >= 0) {
      const frames = lines.map((l) => l.match(/consisting of (\d+) frames|until frame (\d+)/)).find(Boolean);
      const cpus = lines.map((l) => l.match(/Using (\d+) CPUs/)).find(Boolean);
      const step = steps[current];
      const perFrame = step.id === 'interactions' || /frame/i.test(bar?.label ?? '');
      const parts = [step.group ? `${GROUP_TITLES[step.group]}: ${step.label}` : step.label];
      if (bar) parts.push(`${bar.done} / ${bar.total}${perFrame ? ' frames' : ''}`);
      else if (step.id === 'interactions' && frames) parts.push(`${frames[1] ?? frames[2]} frames`);
      if (bar?.rate) parts.push(bar.rate.replace('it/s', '/s').replace('s/it', 's each'));
      if (cpus && step.id === 'interactions') parts.push(`${cpus[1]} CPUs`);
      now.append(h('div', {class: 'oa-now-label'}, h('span', {class: 'oa-now-dot'}), parts.join(' · ')));
      if (bar) now.append(h('div', {class: 'oa-progress'}, h('div', {style: {width: `${bar.percent}%`}})),
        h('div', {class: 'oa-small oa-muted'}, `${bar.label || 'Progress'}: ${bar.percent}%${bar.remaining && bar.remaining !== '?' ? ` · ${bar.remaining} left` : ''}`));
    }

    clear(timeline);
    const row = (s: number, indent: boolean) => {
      const state = s === failed ? 'failed' : done[s] ? 'done' : s === current ? 'current' : 'pending';
      const d = duration(s);
      return h('div', {class: `oa-tl-step ${state} ${indent ? 'sub' : ''}`},
        h('span', {class: 'dot'}, state === 'done' ? '✓' : state === 'failed' ? '✕' : ''),
        h('span', {class: 'label'}, steps[s].label),
        h('span', {class: 'time'}, d !== null && state !== 'pending' ? clock(d) : ''));
    };
    for (let s = 0; s < steps.length; s++) {
      const group = steps[s].group;
      if (!group) {
        timeline.append(row(s, false));
        continue;
      }
      const members: number[] = [];
      while (s < steps.length && steps[s].group === group) members.push(s++);
      s--;
      const doneCount = members.filter((m) => done[m]).length;
      const state = members.includes(failed) ? 'failed' : doneCount === members.length ? 'done' :
        members.includes(current) ? 'current' : 'pending';
      const start = members[0] > 0 ? endTime[members[0] - 1] : 0;
      const end = state === 'done' ? endTime[members[members.length - 1]] : state === 'current' ? run.elapsed : null;
      timeline.append(h('div', {class: `oa-tl-step group ${state}`},
        h('span', {class: 'dot'}, state === 'done' ? '✓' : state === 'failed' ? '✕' : ''),
        h('span', {class: 'label'}, GROUP_TITLES[group], h('span', {class: 'count'}, ` ${doneCount}/${members.length}`)),
        h('span', {class: 'time'}, timed && end !== null && end !== undefined && start !== null && start !== undefined ? clock(end - start) : '')));
      if (state === 'current' || state === 'failed') for (const m of members) timeline.append(row(m, true));
    }
  }

  function render(run: RunSnapshot): void {
    latest = run;
    status.className = `oa-status ${run.status}`;
    status.textContent = run.status;
    time.textContent = clock(run.elapsed);
    if (!command.textContent) command.textContent = run.command;
    out.textContent = `Output folder: ${run.outputDir}`;
    lines.push(...run.lines);
    // Servers before line timestamps send no times; durations are then left out.
    timed = timed && run.times !== undefined;
    times.push(...(run.times ?? run.lines.map(() => run.elapsed)));
    since = run.next;
    if (run.lines.length || run.status !== 'running') renderLog();
    renderSteps(run);
    cancelBtn.style.display = run.status === 'running' ? '' : 'none';
    openBtn.disabled = run.status !== 'finished';
    openBtn.onclick = () => actions.openResults(run.outputDir);
  }

  async function poll(): Promise<void> {
    if (disposed) return;
    try {
      const run = await api.run(runId, since);
      render(run);
      if (run.status === 'running') timer = window.setTimeout(poll, 1000);
      else if (run.status === 'finished') toast('Analysis finished');
    } catch (e) {
      clear(log);
      log.append(h('div', {class: 'err'}, (e as Error).message));
    }
  }
  poll();
  return {el, dispose: () => { disposed = true; clearTimeout(timer); }};
}
