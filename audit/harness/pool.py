"""A process pool with a per-task timeout: a worker that exceeds it is killed and replaced,
and the task is recorded as {'flag': 'timeout'} (never counted as a defect)."""
import multiprocessing as mpc
import os
import time
from multiprocessing.connection import wait


def _worker(conn):
    os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')
    import tasks
    while True:
        try:
            item = conn.recv()
        except EOFError:
            return
        if item is None:
            return
        tid, task = item
        try:
            res = tasks.run(task)
        except Exception as e:
            res = {'flag': 'noref', 'err': repr(e)[:300]}
        conn.send((tid, res))


class Pool:
    def __init__(self, nproc, timeout):
        self.nproc, self.timeout = nproc, timeout
        self.ctx = mpc.get_context('spawn')
        self.workers = [self._spawn() for _ in range(nproc)]

    def _spawn(self):
        a, b = self.ctx.Pipe()
        p = self.ctx.Process(target=_worker, args=(b,), daemon=True)
        p.start()
        return {'p': p, 'c': a, 'task': None, 't0': None}

    def map(self, tasks, timeout=None, progress=None):
        timeout = timeout or self.timeout
        results = [None] * len(tasks)
        todo = list(range(len(tasks)))[::-1]
        done = 0
        last = time.time()
        while done < len(tasks):
            for w in self.workers:
                if w['task'] is None and todo:
                    tid = todo.pop()
                    w['c'].send((tid, tasks[tid]))
                    w['task'], w['t0'] = tid, time.time()
            busy = [w['c'] for w in self.workers if w['task'] is not None]
            ready = wait(busy, timeout=0.5)
            for c in ready:
                w = next(w for w in self.workers if w['c'] is c)
                try:
                    tid, res = c.recv()
                except (EOFError, OSError):
                    tid, res = w['task'], {'flag': 'noref', 'err': 'worker died'}
                    self._replace(w)
                results[tid] = res
                w['task'] = None
                done += 1
            now = time.time()
            for i, w in enumerate(self.workers):
                if w['task'] is not None and now - w['t0'] > timeout:
                    results[w['task']] = {'flag': 'timeout'}
                    done += 1
                    w['p'].kill()
                    self.workers[i] = self._spawn()
            if progress and now - last > 30:
                progress(done, len(tasks))
                last = now
        return results

    def _replace(self, w):
        try:
            w['p'].kill()
        except Exception:
            pass
        i = self.workers.index(w)
        self.workers[i] = self._spawn()
        w.update(self.workers[i])

    def close(self):
        for w in self.workers:
            try:
                w['c'].send(None)
            except Exception:
                pass
        for w in self.workers:
            w['p'].join(timeout=2)
            if w['p'].is_alive():
                w['p'].kill()
