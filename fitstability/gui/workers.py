"""Background execution of the (slow) experiments so the window stays responsive."""
from __future__ import annotations

import traceback

from PySide6.QtCore import QThread, Signal

from ..experiments import Cancelled


class Worker(QThread):
    """Runs ``job(progress, cancel)`` off the GUI thread.

    ``progress(done, total)`` and ``cancel()`` are the callbacks the library experiments accept.
    Exactly one of ``succeeded(result)``, ``cancelled()`` or ``failed(message)`` is emitted.
    """

    progressed = Signal(int, int)
    succeeded = Signal(object)
    cancelled = Signal()
    failed = Signal(str)

    def __init__(self, job, parent=None):
        super().__init__(parent)
        self._job = job
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            out = self._job(lambda d, t: self.progressed.emit(int(d), int(t)), lambda: self._cancel)
        except Cancelled:
            self.cancelled.emit()
        except Exception as e:  # noqa: BLE001 - reported to the user, not swallowed
            self.failed.emit(f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}")
        else:
            if self._cancel:
                self.cancelled.emit()
            else:
                self.succeeded.emit(out)
