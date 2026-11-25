"""In-process task runner with a hard per-task timeout, for running in one process.

A C-level watchdog (faulthandler.dump_traceback_later with exit=True) ends the process if
a task runs longer than `timeout` seconds, even inside a long C call that holds the GIL.
Before each task its key is written to a small local "in-flight" file; the wrapper script
(run.sh) restarts the program, which then records that task as a timeout and skips it.
"""
import faulthandler
import json
import os
import tempfile

# resource limits requested by the owner: two threads per process
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '2'


def state_dir():
    """Local (not synced) folder for the in-flight file and the watchdog log."""
    d = os.environ.get('AUDIT_STATE_DIR') or os.path.join(
        os.environ.get('LOCALAPPDATA', tempfile.gettempdir()), 'certifield_audit')
    os.makedirs(d, exist_ok=True)
    return d


class Runner:
    def __init__(self, name, timeout):
        d = state_dir()
        self.inflight = os.path.join(d, name + '.inflight.json')
        self.log = open(os.path.join(d, name + '.watchdog.log'), 'a')
        self.timeout = timeout

    def take_timed_out(self):
        """Key of the task that was running when a previous process was ended, or None."""
        try:
            with open(self.inflight, encoding='utf-8') as f:
                key = json.load(f)
        except (FileNotFoundError, ValueError):
            return None
        os.remove(self.inflight)
        return key

    def run(self, key, fn, timeout=None):
        with open(self.inflight, 'w', encoding='utf-8') as f:
            json.dump(key, f)
        faulthandler.dump_traceback_later(timeout or self.timeout, exit=True, file=self.log)
        try:
            return fn()
        finally:
            faulthandler.cancel_dump_traceback_later()
            try:
                os.remove(self.inflight)
            except FileNotFoundError:
                pass
