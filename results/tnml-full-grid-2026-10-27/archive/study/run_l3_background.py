"""Run audited L3 points sequentially and refresh the figure after each point.

This is the local launch record for this frozen study, not a new decoder.
The tested full-grid driver controls gates, source identity, limits, and resume.
"""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
repository = root.parent.parent
scripts = repository / 'examples/tnml_repaired_seed'
selection = json.loads((root / 'scan-selection.json').read_text())
julia = '/Users/rongyuqing/.juliaup/bin/julia'
environment = os.environ.copy()
environment.update({
    'JULIA_DEPOT_PATH': '/tmp/toponoise-julia-depot:/Users/rongyuqing/.julia',
    'JULIA_CONDAPKG_BACKEND': 'Null',
    'JULIA_PYTHONCALL_EXE': str(repository / '.CondaPkg/.pixi/envs/default/bin/python'),
    'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1',
    'PYTHONDONTWRITEBYTECODE': '1', 'MPLCONFIGDIR': '/private/tmp/tnml-full-grid-mpl',
})
state = {'pid': os.getpid(), 'selected_sizes': [3], 'requested_sizes': selection['sizes'],
         'status': 'running', 'points': [], 'started_at_utc': datetime.now(timezone.utc).isoformat()}

def save():
    state['updated_at_utc'] = datetime.now(timezone.utc).isoformat()
    temporary = root / 'background-status.json.tmp'
    temporary.write_text(json.dumps(state, indent=2) + '\n')
    temporary.replace(root / 'background-status.json')

save()
for p in selection['p_values']:
    if p == 0:
        continue
    state['current_p'] = p
    save()
    command = [sys.executable, str(scripts/'run_full_grid.py'), '--root', str(root),
               '--phase', 'production', '--sizes', '3', '--p-values', str(p), '--julia', julia]
    result = subprocess.run(command, cwd=repository, env=environment)
    plot = subprocess.run([sys.executable, str(scripts/'plot_full_grid.py'), '--root', str(root)],
                          cwd=repository, env=environment)
    state['points'].append({'p': p, 'driver_returncode': result.returncode,
                            'plot_returncode': plot.returncode})
    if result.returncode != 0:
        state['status'] = 'stopped' if result.returncode == 3 else 'blocked'
        save()
        raise SystemExit(result.returncode)
    if plot.returncode != 0:
        state['status'] = 'report_error'
        save()
        raise SystemExit(plot.returncode)
    save()
state['status'] = 'completed_selected_L3_grid'
save()
