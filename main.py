# main.py
"""
Main entry point for launching the PyBranch GUI application.
"""

import sys
from PyQt6.QtWidgets import QApplication
from pybranch.gui import MainWindow

def main() -> None:
    """Initializes and runs the PyQt application."""
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()