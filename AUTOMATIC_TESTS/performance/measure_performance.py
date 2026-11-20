#!/usr/bin/python3
################################################################################
#                Measure the run time and memory usage of stella               #
################################################################################
# This script records a baseline of the speed and memory usage of stella, so that
# optimisations can be quantified. It is not a pass/fail test, since timings are
# too noisy for that; the numerical tests check that the results did not change.
#
# Each case is one of the input files of the numerical reproducibility tests
# (test 3), with the resolution increased to make the timings meaningful. For each
# case we record:
#     - the wall time of the simulation (the minimum over <repeats> runs),
#     - the peak memory (maximum resident set size) of every MPI process,
#     - the timers which stella reports itself (print_extra_info_to_terminal),
#       e.g. the time spent in each term of the gyrokinetic equation.
#
# Usage (from the main stella directory):
#     python3 AUTOMATIC_TESTS/performance/measure_performance.py --output baseline.json
#     ... modify stella and recompile ...
#     python3 AUTOMATIC_TESTS/performance/measure_performance.py --compare baseline.json
#
# Run the baseline and the comparison on the same machine, with the same number of
# processors, and with as few other programs running as possible.
################################################################################

import os
import re
import sys
import json
import time
import shutil
import socket
import argparse
import pathlib
import platform
import tempfile
import statistics
import subprocess

performance_directory = pathlib.Path(__file__).absolute().parent
input_directory = performance_directory.parent / 'numerical_tests' / 'test_3_numerical_reproducibility'
default_stella_path = performance_directory.parent.parent / 'stella'

# Electrostatic and electromagnetic flux tube cases
default_cases = ['es_nonlinear', 'em_nonlinear', 'adiabatic_electrons',
                 'collisions_dougherty_implicit', 'collisions_fokker_planck_implicit']

# Resolutions: the 'small' resolution is the one of the numerical tests
resolutions = {
    'small':  {},
    'medium': {'nx': 16, 'ny': 16, 'nzed': 24, 'nmu': 6, 'nvgrid': 8, 'nstep': 100, 'nwrite': 10},
    'large':  {'nx': 32, 'ny': 32, 'nzed': 32, 'nmu': 8, 'nvgrid': 12, 'nstep': 100, 'nwrite': 10},
}

#-------------------------------------------------------------------------------
def write_input_file(case, resolution, folder):
    '''Copy the input file of <case> to <folder>, with the given resolution, and
    with print_extra_info_to_terminal = .true. so that stella prints its timers.'''
    text = open(input_directory / f'{case}.in', 'r').read()
    for key, value in resolutions[resolution].items():
        text, n = re.subn(rf'^(\s*{key}\s*=\s*)\S+', rf'\g<1>{value}', text, flags=re.M)
        if n != 1: sys.exit(f'ERROR: expected one "{key} = ..." in {case}.in, found {n}.')
    text, n = re.subn(r'print_extra_info_to_terminal\s*=\s*\S+', 'print_extra_info_to_terminal = .true.', text)
    if n != 1: sys.exit(f'ERROR: expected "print_extra_info_to_terminal" in {case}.in.')
    (folder / f'{case}.in').write_text(text)
    return f'{case}.in'

#-------------------------------------------------------------------------------
def run_stella(stella_path, input_file, folder, nproc):
    '''Run stella in <folder>. Each MPI process is started through this script
    (with --measure-rss), which records the peak memory of that stella process.'''
    for f in folder.glob('rss_*.txt'): f.unlink()
    command = ['mpirun', '--oversubscribe', '-np', str(nproc), sys.executable, str(pathlib.Path(__file__).absolute()),
               '--measure-rss', str(folder), str(stella_path), input_file]
    start = time.perf_counter()
    result = subprocess.run(command, cwd=folder, capture_output=True, text=True)
    wall_time = time.perf_counter() - start
    if result.returncode != 0:
        print(result.stdout[-3000:], result.stderr[-3000:])
        sys.exit(f'ERROR: stella failed for {input_file}.')
    peak_memory_mb = sorted(float(f.read_text()) for f in folder.glob('rss_*.txt'))
    return wall_time, peak_memory_mb, result.stdout

#-------------------------------------------------------------------------------
def measure_rss(folder, command):
    '''Run <command> and write its peak resident set size (in MB) to <folder>.'''
    import resource
    returncode = subprocess.run(command).returncode
    maxrss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    maxrss_mb = maxrss / 1024**2 if platform.system() == 'Darwin' else maxrss / 1024  # bytes on macOS, kB on Linux
    pathlib.Path(folder, f'rss_{os.getpid()}.txt').write_text(f'{maxrss_mb:.3f}')
    sys.exit(returncode)

#-------------------------------------------------------------------------------
def parse_stella_timers(stdout):
    '''Read the "ELAPSED TIME" report of stella, e.g. {'GYROKINETIC EQUATION/ExB nonlin': 0.12}.'''
    timers, section = {}, None
    report = stdout.split('ELAPSED TIME', 1)[-1] if 'ELAPSED TIME' in stdout else ''
    lines = report.split('\n')
    for i, line in enumerate(lines):
        if i + 1 < len(lines) and re.match(r'^\s*-{3,}\s*$', lines[i + 1]) and line.strip():
            section = line.strip()
        match = re.match(r'^\s*(.+?):\s+([0-9.]+)\s+(sec|min|hours)', line)
        if match and section:
            factor = {'sec': 1, 'min': 60, 'hours': 3600}[match.group(3)]
            key = f'{section}/{match.group(1)}'
            while key in timers: key += "'"  # some names appear twice in a section
            timers[key] = float(match.group(2)) * factor
    return timers

#-------------------------------------------------------------------------------
def git_info():
    try:
        commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True, cwd=performance_directory).stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', '../../STELLA_CODE'], capture_output=True, text=True, cwd=performance_directory).stdout.strip() != ''
        return commit + ('+modified' if dirty else '')
    except Exception:
        return 'unknown'

#-------------------------------------------------------------------------------
def measure(args):
    results = {'commit': git_info(), 'date': time.strftime('%Y-%m-%d %H:%M'), 'host': socket.gethostname(),
               'nproc': args.nproc, 'resolution': args.resolution, 'repeats': args.repeats, 'cases': {}}
    print(f'Measuring stella ({results["commit"]}) on {args.nproc} processors, resolution "{args.resolution}":')
    for case in args.cases:
        with tempfile.TemporaryDirectory() as folder:
            folder = pathlib.Path(folder)
            input_file = write_input_file(case, args.resolution, folder)
            runs = [run_stella(args.stella, input_file, folder, args.nproc) for _ in range(args.repeats)]
        wall_times = [run[0] for run in runs]
        best = min(range(len(runs)), key=lambda i: wall_times[i])
        results['cases'][case] = {
            'wall_time_s': min(wall_times),
            'wall_time_median_s': statistics.median(wall_times),
            'peak_memory_per_process_mb': max(runs[best][1]),
            'peak_memory_total_mb': sum(runs[best][1]),
            'stella_timers_s': parse_stella_timers(runs[best][2]),
        }
        r = results['cases'][case]
        print(f'  {case:36s} {r["wall_time_s"]:8.2f} s   {r["peak_memory_per_process_mb"]:8.1f} MB/process   {r["peak_memory_total_mb"]:8.1f} MB total')
    return results

#-------------------------------------------------------------------------------
def compare(baseline, current, threshold):
    '''Print the relative change of every quantity with respect to the baseline.'''
    print(f'\nComparison with the baseline ({baseline["commit"]}, {baseline["date"]}, {baseline["host"]}):')
    for key in ['nproc', 'resolution']:
        if baseline[key] != current[key]: print(f'  WARNING: {key} differs ({baseline[key]} vs {current[key]}), the comparison is not meaningful.')
    if baseline['host'] != current['host']: print('  WARNING: the baseline was measured on a different machine.')
    def change(old, new):
        if old == 0: return ''
        percent = 100 * (new - old) / old
        flag = '  <--' if abs(percent) > threshold else ''
        return f'{percent:+7.1f}%{flag}'
    for case, new in current['cases'].items():
        old = baseline['cases'].get(case)
        if old is None: print(f'\n  {case}: not in the baseline'); continue
        print(f'\n  {case}')
        print(f'    {"quantity":46s} {"baseline":>10s} {"current":>10s}   change')
        for key, label in [('wall_time_s', 'wall time [s]'), ('peak_memory_per_process_mb', 'peak memory per process [MB]'), ('peak_memory_total_mb', 'peak memory total [MB]')]:
            print(f'    {label:46s} {old[key]:10.2f} {new[key]:10.2f}  {change(old[key], new[key])}')
        for key in new['stella_timers_s']:
            if key not in old['stella_timers_s']: continue
            a, b = old['stella_timers_s'][key], new['stella_timers_s'][key]
            if max(a, b) < 0.01 * new['wall_time_s']: continue  # skip timers below 1% of the wall time
            print(f'    {key[:46]:46s} {a:10.2f} {b:10.2f}  {change(a, b)}')
    print(f'\n  Changes larger than {threshold}% are marked with "<--". Timings vary between runs, so')
    print(f'  check that a change is reproducible (e.g. with --repeats 5) before drawing conclusions.')

#-------------------------------------------------------------------------------
if __name__ == '__main__':

    # Internal mode: run a single stella process and record its peak memory
    if len(sys.argv) > 1 and sys.argv[1] == '--measure-rss':
        measure_rss(sys.argv[2], sys.argv[3:])

    parser = argparse.ArgumentParser(description='Measure the run time and memory usage of stella.')
    parser.add_argument('--stella', type=pathlib.Path, default=pathlib.Path(os.environ.get('STELLA_EXE_PATH', default_stella_path)), help='stella executable')
    parser.add_argument('--nproc', type=int, default=4, help='number of MPI processes (default: 4)')
    parser.add_argument('--resolution', choices=resolutions.keys(), default='medium', help='resolution of the cases (default: medium)')
    parser.add_argument('--repeats', type=int, default=3, help='number of runs per case, the fastest is kept (default: 3)')
    parser.add_argument('--cases', nargs='+', default=default_cases, help=f'input files of test 3 (default: {" ".join(default_cases)})')
    parser.add_argument('--output', type=pathlib.Path, help='write the results to this json file')
    parser.add_argument('--compare', type=pathlib.Path, help='compare the results with this json file')
    parser.add_argument('--threshold', type=float, default=5, help='mark changes larger than this percentage (default: 5)')
    args = parser.parse_args()
    args.stella = args.stella.absolute()
    if not args.stella.exists(): sys.exit(f'ERROR: the stella executable {args.stella} does not exist.')

    # When comparing, use the same settings as the baseline by default
    baseline = None
    if args.compare:
        baseline = json.loads(args.compare.read_text())
        if '--nproc' not in sys.argv: args.nproc = baseline['nproc']
        if '--resolution' not in sys.argv: args.resolution = baseline['resolution']
        if '--cases' not in sys.argv: args.cases = list(baseline['cases'].keys())

    results = measure(args)
    if args.output:
        args.output.write_text(json.dumps(results, indent=2))
        print(f'\nWrote the results to {args.output}')
    if baseline: compare(baseline, results, args.threshold)
