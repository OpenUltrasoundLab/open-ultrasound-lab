#!/usr/bin/env python3
"""Run MATLAB scripts in GNU Octave (headless) and save every figure as PNG.
Usage: python3 mtest.py matlab/m02_svd_cyst.m [more.m ...]
Prints the script's stdout/stderr, wall time, and the PNG paths (shots/matlab/<script>-figN.png).
Exit code 1 if any script errored. Scripts must be base-MATLAB compatible (see CONTRIBUTING.md)."""
import sys, subprocess, pathlib, time
root = pathlib.Path(__file__).parent
out = root / 'shots' / 'matlab'; out.mkdir(parents=True, exist_ok=True)
bad = 0
for f in sys.argv[1:]:
    p = pathlib.Path(f).resolve(); stem = p.stem
    for old in out.glob(stem + '-fig*.png'): old.unlink()
    prelude = ("graphics_toolkit('gnuplot'); set(0,'defaultfigurevisible','off'); "
               "warning('off','all'); more off; ")
    post = (f"figs = sort(get(0,'children')); for k=1:numel(figs), "
            f"print(figs(k), sprintf('{out}/{stem}-fig%d.png', k), '-dpng', '-r80'); end")
    cmd = prelude + f"cd('{p.parent}'); try, run('{p.name}'); catch err, disp(['ERROR: ' err.message]); "\
          "for s = err.stack', fprintf('  at %s line %d\\n', s.name, s.line); end; exit(1); end; " + post
    t0 = time.time()
    r = subprocess.run(['octave', '--no-gui', '-q', '--eval', cmd], capture_output=True, text=True, timeout=600)
    dt = time.time() - t0
    noise = ('gnuplot', 'qt toolkit', 'iconv', 'internal list of open figure', 'be notified', 'one-directional pipe', 'Octave interpreter so most', 'will not be reflected', 'example, if the plot window', 'not update its internal')
    err = '\n'.join(l for l in r.stderr.splitlines() if not any(n in l for n in noise))
    print(f'===== {p.name}  ({dt:.1f} s, exit {r.returncode})')
    print(r.stdout.rstrip()[-6000:])
    if err.strip(): print('[stderr]', err[-3000:])
    for png in sorted(out.glob(stem + '-fig*.png')): print('  figure:', png)
    bad += r.returncode != 0 or 'ERROR:' in r.stdout
sys.exit(1 if bad else 0)
