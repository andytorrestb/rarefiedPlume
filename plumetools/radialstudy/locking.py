"""Operating-system locks released automatically when a controller dies."""
from contextlib import contextmanager
import os


@contextmanager
def exclusive_lock(path):
    with path.open('a+b') as stream:
        if os.name=='nt':
            import msvcrt
            if stream.seek(0,2)==0:
                stream.write(b'0')
                stream.flush()
            stream.seek(0)
            lock = lambda:msvcrt.locking(stream.fileno(),msvcrt.LK_NBLCK,1)
            unlock = lambda:msvcrt.locking(stream.fileno(),msvcrt.LK_UNLCK,1)
        else:
            import fcntl
            lock = lambda:fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            unlock = lambda:fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
        try:
            lock()
        except OSError as exc:
            raise RuntimeError(f'a controller or native segment is already running (locked: {path})') from exc
        try:
            yield
        finally:
            unlock()
