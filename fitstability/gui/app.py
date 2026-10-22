"""Entry point: ``fitstability-gui [file]`` or ``python -m fitstability.gui [file]``."""
from __future__ import annotations

import argparse
import multiprocessing
import sys


def main(argv=None) -> int:
    multiprocessing.freeze_support()
    ap = argparse.ArgumentParser(prog="fitstability-gui", description=__doc__)
    ap.add_argument("file", nargs="?", help="VAMAS or CSV spectrum to open")
    ap.add_argument("--theme", choices=["light", "dark"], help="colour theme (remembered)")
    args = ap.parse_args(argv)

    from PySide6.QtWidgets import QApplication
    from .main_window import MainWindow
    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("fitstability")
    win = MainWindow(args.theme)
    win.show()
    if args.file:
        win.open_file(args.file)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
