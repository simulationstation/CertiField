"""Resource limits of the owner's machine: import this module FIRST in every entry script (before numpy).

At most 2 threads per process for the numerical libraries; worker pools of at most BENCHPACK_NPROC (default 2)
processes (see p2rig.evaluate).  Values already set in the environment are kept.
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '2')
os.environ.setdefault('BENCHPACK_NPROC', '2')
os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')

WORK_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'work')
