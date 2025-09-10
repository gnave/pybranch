
# pybranch/gui.py
"""
Defines the PyQt6 graphical user interface for the interactive branching
fraction calculator.
"""
import os
import sys
import logging
from typing import Optional, List, Dict, Any
from .analysis import (
    BranchingRatioAnalysis, format_intensity_table_as_string,
    format_final_results_as_string
)

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QComboBox, QLineEdit, QTextEdit, QLabel, QGridLayout,
    QInputDialog, QMessageBox, QFileDialog
)
from PyQt6.QtGui import QFont

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from .analysis import BranchingRatioAnalysis, format_intensity_table_as_string, calculate_final_results
from .utils import load_config
from .data_io import read_lifetimes, glob, Spectrum
from .visualization import plot_spectrum_line


class MainWindow(QMainWindow):
    """The main application window for the Branching Fraction Calculator."""

    def __init__(self):
        super().__init__()
        self.analysis: Optional[BranchingRatioAnalysis] = None
        self.all_levels: List[str] = []
        self.all_spectra_files: List[str] = []
        self.logging_configured: bool = False
        
        # --- FIX 1: Resolve all paths in the config to be absolute ---
        self.config: Dict[str, Any] = self._load_and_resolve_paths()

        self.initUI()
        self.load_initial_data()

    def _load_and_resolve_paths(self) -> Dict[str, Any]:
        """Loads config and resolves all file paths to be absolute."""
        config = load_config()
        # Get the directory where this gui.py script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Assume the project root (where config.yaml is) is one level up
        project_root = os.path.dirname(script_dir)
        
        files_config = config.get('files', {})
        for key, value in files_config.items():
            if not os.path.isabs(value):
                # Construct an absolute path from the project root
                absolute_path = os.path.join(project_root, value)
                files_config[key] = absolute_path
        return config

    def initUI(self) -> None:
        """Initializes the user interface, widgets, and layout."""
        # (This method is largely the same, but the connections are confirmed correct)
        self.setWindowTitle("Branching Fraction Calculator (PyQt Edition)")
        self.setGeometry(100, 100, 1200, 750)

        central_widget = QWidget()
        main_layout = QHBoxLayout(central_widget)
        self.setCentralWidget(central_widget)

        controls_panel = QWidget()
        controls_layout = QVBoxLayout(controls_panel)
        controls_panel.setFixedWidth(250)

        top_button_layout = QHBoxLayout()
        self.quit_button = QPushButton("Quit")
        self.results_button = QPushButton("Results")
        top_button_layout.addWidget(self.quit_button)
        top_button_layout.addWidget(self.results_button)
        controls_layout.addLayout(top_button_layout)

        grid_layout = QGridLayout()
        grid_layout.setSpacing(10)

        # Log File
        self.log_file_edit = self._add_grid_widget(grid_layout, 0, "Log file name:")
        
        # Upper Level (Label + ComboBox)
        self.upper_level_combo = self._add_grid_widget(grid_layout, 1, "Upper level:", QComboBox())
        
        # Other controls...
        self.ref_level_button = QPushButton("Ref. Level:")
        self.ref_level_combo = QComboBox()
        grid_layout.addWidget(self.ref_level_button, 2, 0)
        grid_layout.addWidget(self.ref_level_combo, 2, 1)

        self.rescale_file_button = QPushButton("Rescale File:")
        self.rescale_file_combo = QComboBox()
        grid_layout.addWidget(self.rescale_file_button, 3, 0)
        grid_layout.addWidget(self.rescale_file_combo, 3, 1)

        self.rescale_level_button = QPushButton("Rescale Level:")
        self.rescale_level_combo = QComboBox()
        grid_layout.addWidget(self.rescale_level_button, 4, 0)
        grid_layout.addWidget(self.rescale_level_combo, 4, 1)

        self.delete_file_button = QPushButton("Delete file:")
        self.delete_file_combo = QComboBox()
        grid_layout.addWidget(self.delete_file_button, 5, 0)
        grid_layout.addWidget(self.delete_file_combo, 5, 1)

        self.delete_level_button = QPushButton("Delete level:")
        self.delete_level_combo = QComboBox()
        grid_layout.addWidget(self.delete_level_button, 6, 0)
        grid_layout.addWidget(self.delete_level_combo, 6, 1)

        controls_layout.addLayout(grid_layout)
        controls_layout.addStretch(1)

        # Plotting Controls
        plot_layout = QGridLayout()
        plot_label = QLabel("Plot Line:")
        self.plot_level_combo = QComboBox()
        self.plot_file_button = QPushButton("Select File & Plot")
        plot_layout.addWidget(plot_label, 0, 0)
        plot_layout.addWidget(self.plot_level_combo, 0, 1)
        plot_layout.addWidget(self.plot_file_button, 1, 1)
        controls_layout.addLayout(plot_layout)

        self.comment_button = QPushButton("Add comment")
        controls_layout.addWidget(self.comment_button)

        main_layout.addWidget(controls_panel)

        self.display_text = QTextEdit()
        self.display_text.setReadOnly(True)
        self.display_text.setFont(QFont("Courier New", 9))
        main_layout.addWidget(self.display_text)

        # --- Connect Signals to Slots ---
        self.quit_button.clicked.connect(self.close)
        self.log_file_edit.returnPressed.connect(self._on_set_log_file)
        self.upper_level_combo.currentTextChanged.connect(self._on_upper_level_selected)
        self.results_button.clicked.connect(self._on_calculate_results)
        self.comment_button.clicked.connect(self._on_add_comment)
        self.plot_file_button.clicked.connect(self._on_plot_line)
        self.ref_level_button.clicked.connect(self._on_normalize)
        self.rescale_level_button.clicked.connect(self._on_rescale)
        self.delete_level_button.clicked.connect(self._on_delete)

    def _add_grid_widget(self, grid, row, label_text, widget=None):
        label = QLabel(label_text)
        if widget is None: widget = QLineEdit()
        grid.addWidget(label, row, 0)
        grid.addWidget(widget, row, 1)
        return widget

    def load_initial_data(self) -> None:
        """Loads data needed to populate the GUI for the first time."""
        default_log = self.config.get('defaults', {}).get('log_filename', 'analysis')
        self.log_file_edit.setText(default_log)
        self._on_set_log_file()

        files_config = self.config.get('files', {})
        # Because paths are now absolute, glob will work reliably
        _, _, self.all_levels = read_lifetimes(files_config.get('levels_glob', '*.lev'))
        self.all_spectra_files = glob(files_config.get('spectrum_glob', '*.II'))
        
        # Block signals while we set up the combo box to prevent premature firing
        self.upper_level_combo.blockSignals(True)
        self.upper_level_combo.clear()
        self.upper_level_combo.addItems([""] + sorted(self.all_levels))
        self.upper_level_combo.blockSignals(False)
        
        # --- FIX 2: Manually set the default and trigger the first analysis ---
        default_upper = self.config.get('defaults', {}).get('upper_level')
        if default_upper in self.all_levels:
            self.upper_level_combo.setCurrentText(default_upper)
        else:
            # If no default, manually trigger with a blank to clear the menus
            self._on_upper_level_selected("")

    def _on_set_log_file(self) -> None:
        """Configures the logging system based on the QLineEdit content."""
        log_name = self.log_file_edit.text()
        if not log_name: return
        
        logging.basicConfig(level=logging.INFO, 
                            format='%(message)s', # Simpler format for log file
                            filename=f"{log_name}.log",
                            filemode='w', force=True)
        
        if not self.logging_configured:
            self._log(f"Logging initialized to file: {log_name}.log")
            self.logging_configured = True
        else:
            self._log(f"\nLog file changed to: {log_name}.log")

    def _on_upper_level_selected(self, level: str) -> None:
        """The single entry point to start or reset an analysis for a given level."""
        self.display_text.clear()
        if not level:
            self.analysis = None
            self._log("Please select an upper level to begin analysis.")
            self._populate_menus() # Call with no analysis to clear menus
            return
            
        self._log(f"--- Starting analysis for Upper Level: {level} ---")
        
        self.analysis = BranchingRatioAnalysis(upper_level=level, config=self.config)
        self.analysis.load_atomic_data()
        self.analysis.aggregate_observed_data()

        self._update_display("Data After Aggregation (Raw Values)")
        self._populate_menus()

    def _update_display(self, title: str) -> None:
        """Appends a formatted data table to the main display."""
        # --- FIX 3: This will now display the real table ---
        if self.analysis:
            table_string = format_intensity_table_as_string(self.analysis, title)
            self._log(table_string)
        else:
            self._log(f"\n--- {title} ---\n(No analysis data to display)")

    def _populate_menus(self) -> None:
        """Populates all combo boxes based on the current analysis data."""
        # This function will now be called correctly and have data to work with.
        # Block signals to prevent changes from triggering events
        for combo in [self.ref_level_combo, self.rescale_level_combo, self.delete_level_combo, self.plot_level_combo]:
            combo.blockSignals(True)
            combo.clear()

        # Populate with data if an analysis is active
        if self.analysis and self.analysis.transition_ids:
            lower_levels = sorted(list(set(self.analysis.transition_ids.values())))
            wavenumber_items = [f"{wn:.3f} ({ll})" for wn, ll in sorted(self.analysis.transition_ids.items(), reverse=True)]

            self.ref_level_combo.addItems(lower_levels)
            self.rescale_level_combo.addItems(lower_levels)
            self.delete_level_combo.addItems(lower_levels)
            self.plot_level_combo.addItems(wavenumber_items)

        # File-based menus are populated regardless of analysis state
        self.rescale_file_combo.clear()
        self.delete_file_combo.clear()
        file_basenames = sorted([os.path.basename(f) for f in self.all_spectra_files])
        self.rescale_file_combo.addItems(file_basenames)
        self.delete_file_combo.addItems(file_basenames)

        for combo in [self.ref_level_combo, self.rescale_level_combo, self.delete_level_combo, self.plot_level_combo]:
            combo.blockSignals(False)

    def _log(self, text: str) -> None:
        """Appends text to the display and writes it to the configured log file."""
        self.display_text.append(text)
        self.display_text.verticalScrollBar().setValue(self.display_text.verticalScrollBar().maximum())
        if self.logging_configured:
            logging.info(text)
            
    # --- Other methods (_on_normalize, _on_rescale, etc.) are unchanged ---
    # They should now work correctly because the menus will be populated.

    def _on_normalize(self) -> None:
        if not self.analysis: return
        ref_level = self.ref_level_combo.currentText()
        if not ref_level: return
        self.analysis.normalize_by_reference_line(ref_level)
        self._update_display(f"Data After Normalizing by '{ref_level}'")

    def _on_rescale(self) -> None:
        if not self.analysis: return
        rescale_file = self.rescale_file_combo.currentText()
        rescale_level = self.rescale_level_combo.currentText()
        ref_level = self.ref_level_combo.currentText()
        if not all([rescale_file, rescale_level, ref_level]): return
        
        # Find the full filename to match the key
        full_rescale_file = next((f for f in self.all_spectra_files if os.path.basename(f) == rescale_file), None)
        if not full_rescale_file:
            self._log(f"Error: Could not find full path for {rescale_file}")
            return

        self.analysis.rescale_spectrum_by_transfer_line(
            transfer_level=rescale_level,
            reference_file=full_rescale_file,
            initial_reference_level=ref_level
        )
        self._update_display(f"Data After Rescaling '{os.path.basename(full_rescale_file)}'")

    
    def _on_delete(self) -> None:
        """
        Handles deleting a single data point and automatically removes the
        entire transition if it was the last remaining measurement.
        """
        if not self.analysis:
            return
        file_to_delete_from = self.delete_file_combo.currentText()
        level_to_delete = self.delete_level_combo.currentText()
        if not all([file_to_delete_from, level_to_delete]):
            return

        # Find the full, absolute file path from the selected basename
        full_file_path = next((f for f in self.all_spectra_files if os.path.basename(f) == file_to_delete_from), None)
        if not full_file_path:
            self._log(f"Error: Could not find full path for {file_to_delete_from}")
            return

        # Step 1: Perform the simple deletion of the single data point
        self.analysis.delete_line(full_file_path, level_to_delete)

        # Step 2: Check if any measurements of this transition remain across all files
        remaining_count = sum(1 for key in self.analysis.intensities if key[1] == level_to_delete)

        display_title = "Data After Deleting Line"

        # Step 3: If no measurements remain, promote the deletion to a full transition removal
        if remaining_count == 0:
            self._log(f"Info: Last measurement of '{level_to_delete}' was deleted. Removing transition entirely from analysis.")
            
            # We must find the wavenumber that corresponds to this lower level to remove it
            wavenumber_to_remove = next((wn for wn, ll in self.analysis.transition_ids.items() if ll == level_to_delete), None)
            
            if wavenumber_to_remove:
                self.analysis.remove_transition_entirely(wavenumber_to_remove)
                display_title = "Data After Deleting Line (Transition Removed)"
                # Crucially, we must update the dropdown menus now that a line has vanished
                self._populate_menus()
            else:
                self._log(f"Warning: Could not find a wavenumber for '{level_to_delete}' to remove it from the transition list.")

        # Step 4: Update the display to show the result of the operation
        self._update_display(display_title)

    def _on_calculate_results(self) -> None:
        """Calculates and displays the final branching fraction results."""
        if not self.analysis:
            self._log("\nError: Please select an upper level and load data first.")
            return
            
        # Run the calculation
        results = self.analysis.calculate_branching_ratios()

        # --- THIS IS THE FIX ---
        # Format the entire results dictionary into a string
        results_string = format_final_results_as_string(results)
        
        # Log the formatted string to the display and the log file
        self._log(results_string)

    def _on_add_comment(self) -> None:
        """Opens a dialog to add a comment to the log."""
        text, ok = QInputDialog.getText(self, 'Add Comment', 'Enter your comment:')
        if ok and text:
            self._log(f"\n# --- User Comment ---")
            self._log(f"# {text}")
            self._log(f"# --------------------")

  # Replace the existing _on_plot_line method with this one.

    def _on_plot_line(self) -> None:
        """Opens a file dialog and plots the selected line, applying corrections."""
        if not self.analysis:
            QMessageBox.warning(self, "Warning", "Please select an upper level first.")
            return
            
        wavenumber_to_plot_str = self.plot_level_combo.currentText()
        if not wavenumber_to_plot_str:
            QMessageBox.warning(self, "Warning", "Please select a line to plot.")
            return
        
        # This wavenumber is from the dropdown, so it's already in the corrected scale
        wavenumber = float(wavenumber_to_plot_str.split('(')[0].strip())
        
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Spectrum File", "", "Data Files (*.dat)")
        if not file_path:
            return
        
        base_path = file_path[:-4]
        try:
            spectrum = Spectrum(base_path)
            spectrum.load()
            
            # --- THIS IS THE FIX ---
            # Get the wavenumber correction factor from the loaded header
            wavcorr = float(spectrum.header.get('wavcorr', 0.0))
            # --- END OF FIX ---

            window_size = self.config.get('plotting',{}).get('window_length', 32)
            
            # Pass the wavcorr value to the plotting function
            fig = plot_spectrum_line(spectrum, wavenumber, window_size, wavcorr_applied=wavcorr)
            
            if fig:
                self.plot_window = PlotWindow(fig)
                self.plot_window.show()
            else:
                QMessageBox.critical(self, "Error", "Failed to generate plot.")

        except FileNotFoundError:
            QMessageBox.critical(self, "Error", f"Could not find associated files for {base_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred during plotting: {e}")

class PlotWindow(QMainWindow):
    """A simple window for displaying a Matplotlib figure."""
    def __init__(self, fig: Figure):
        super().__init__()
        self.setWindowTitle("Spectrum Plot")
        canvas = FigureCanvas(fig)
        self.setCentralWidget(canvas)