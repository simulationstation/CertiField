"""Build boost/check.cpp against two Boost.Math source trees and evaluate the same calls
with both (the commit SciPy bundles, and current develop).

    python boostcheck.py --old <math-old>/include --new <math-new>/include build
    python boostcheck.py run < calls.txt          (lines: <name> <args...>)

Windows/MSVC: vcvars64.bat is located under "C:/Program Files/Microsoft Visual Studio".
Binaries go to the local state folder (not the repository).
"""
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runner  # noqa: E402

OUT = os.path.join(runner.state_dir(), 'boost')


def vcvars():
    c = sorted(glob.glob(r'C:\Program Files\Microsoft Visual Studio\*\*\VC\Auxiliary\Build\vcvars64.bat'))
    if not c:
        raise SystemExit('vcvars64.bat not found')
    return c[-1]


def build(inc, exe):
    os.makedirs(OUT, exist_ok=True)
    src = os.path.join(HERE, 'boost', 'check.cpp')
    inner = ('"%s" >nul && cl /nologo /EHsc /O2 /std:c++17 /fp:precise /I "%s" "%s" /Fe:"%s" /Fo:"%s\\\\"'
             % (vcvars(), inc, src, exe, OUT))
    # one string, so that cmd sees the quotes unchanged (/s: strip only the outer pair)
    r = subprocess.run('cmd /s /c "%s"' % inner, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise SystemExit('build failed')


def run(lines):
    outs = {}
    for tag in ('old', 'new'):
        exe = os.path.join(OUT, 'check_%s.exe' % tag)
        r = subprocess.run([exe], input='\n'.join(lines) + '\n', capture_output=True, text=True)
        outs[tag] = r.stdout.splitlines()
    return outs


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    if len(sys.argv) > 1 and sys.argv[-1] == 'build':
        for tag in ('old', 'new'):
            if '--' + tag in sys.argv:
                build(sys.argv[sys.argv.index('--' + tag) + 1], os.path.join(OUT, 'check_%s.exe' % tag))
        print('built in', OUT)
        return
    lines = [l.strip() for l in sys.stdin if l.strip()]
    outs = run(lines)
    for a, b in zip(outs['old'], outs['new']):
        print('old: ' + a)
        print('new: ' + b.split(' -> ', 1)[-1].join(['    -> ', '']) if ' -> ' in b else 'new: ' + b)


if __name__ == '__main__':
    main()
