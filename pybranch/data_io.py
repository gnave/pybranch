# pybranch/data_io.py
"""Handles all data input and output operations for the PyBranch project.

This module is the data access layer for the application. It contains functions
for reading various specialized data formats (e.g., .dat, .hdr, .lev) and
encapsulates the data for a single observation within the `Spectrum` class.
"""

import numpy as np
from struct import unpack
from os import path
from glob import glob
from typing import List, Dict, Any, Tuple, Optional

# --- Main Data Class ---

class Spectrum:
    """Represents a single atomic spectrum and all its associated data.

    This class acts as a container for the data loaded from a set of related
    files (.hdr, .dat, .lin). Upon initialization, it holds the path to the
    data. The `load()` method must be called to populate the data attributes.

    Attributes:
        path_base (str): The base path to the spectrum files, without extension.
        header (Dict[str, Any]): A dictionary of metadata loaded from the .hdr file.
        raw_data (Optional[np.ndarray]): The raw, unprocessed data array loaded
            from the .dat file.
        linelist (List[Dict[str, Any]]): A list of dictionaries, where each
            dictionary contains fitted line parameters from the .lin file.
        wavenumbers (Optional[np.ndarray]): The calculated wavenumber (x-axis)
            of the spectrum. Populated after `load()` is called.
        intensities (Optional[np.ndarray]): The calculated intensity (y-axis)
            of the spectrum. Populated after `load()` is called.
    """

    def __init__(self, specfile_path_base: str):
        """Initializes the Spectrum object.

        Args:
            specfile_path_base (str): The base path to the spectrum files,
                e.g., 'data/Cr102700.003.I' (without extension).

        Raises:
            FileNotFoundError: If the required .dat file for the given base
                path does not exist.
        """

        if not path.exists(specfile_path_base + ".dat"):
            raise FileNotFoundError(f"Required .dat file not found for base path: {specfile_path_base}")
            
        self.path_base: str = specfile_path_base
        self.header: Dict[str, Any] = {}
        self.raw_data: Optional[np.ndarray] = None
        self.linelist: List[Dict[str, Any]] = []
        self.wavenumbers: Optional[np.ndarray] = None
        self.intensities: Optional[np.ndarray] = None

    def load(self):
        """Loads all associated data for this spectrum (.hdr, .dat, .lin).

        This method populates the `header`, `raw_data`, `linelist`,
        `wavenumbers`, and `intensities` attributes by calling the
        specialized file-reading functions.
        """

        self._load_header()
        self._load_data()
        self._load_linelist()
        self._calculate_axes()
        print(f"Successfully loaded spectrum: {path.basename(self.path_base)}")

    def _load_header(self):
        """Private method to load metadata from the .hdr file."""
        self.header = read_header(self.path_base)

    def _load_data(self):
        """Private method to load intensity data from the .dat file."""
        if not self.header:
            self._load_header()
        self.raw_data = read_spectrum_data(self.path_base, self.header)

    def _load_linelist(self):
        """Private method to load fitted line parameters from the .lin file."""
        self.linelist = read_linelist(self.path_base)
        
    def _calculate_axes(self):
        """Private method to calculate the final wavenumber and intensity axes."""
        if self.raw_data is None:
            return
            
        # Extract parameters from header with defaults
        is_complex = 'Complex' in self.header.get('data_is', '')
        cmplx = 2 if is_complex else 1
        
        wstart = float(self.header.get('wstart', 0.0))
        delw = float(self.header.get('delw', 1.0))
        wavcorr = float(self.header.get('wavcorr', 0.0))
        rdsclfct = float(self.header.get('rdsclfct', 1.0))
        npts = len(self.raw_data) // cmplx
        
        # Apply wavenumber correction
        if wavcorr != 0:
            wstart = wstart * (1 + wavcorr)
            delw = delw * (1 + wavcorr)

        self.wavenumbers = np.linspace(wstart, wstart + (npts - 1) * delw, npts)
        self.intensities = self.raw_data[::cmplx] * rdsclfct


# --- File Reading Helper Functions ---

def read_header(specfile_path_base: str) -> Dict[str, str]:
    """
    Reads a .hdr file and creates a metadata dictionary.
    """
    header = {}
    hdr_path = specfile_path_base + ".hdr"
    val = "" # Initialize val to handle 'continue' lines safely
    if not path.exists(hdr_path):
        print(f"Warning: Header file not found: {hdr_path}")
        return header

    with open(hdr_path, 'r') as hdr:
        for line in hdr:
            if line.startswith("/") or line.startswith("END"):
                continue
            if line.startswith("continue"):
                val += line[9:32].strip()
            else:
                key = line[0:8].rstrip()
                val = line[9:32].strip()
            if line.startswith("id"):
                val = line[9:80].strip()
            header[key] = val
            header[key + "_comment"] = line[34:80].strip()
    return header

def read_spectrum_data(specfile_path_base: str, header: Dict[str, Any]) -> np.ndarray:
    """
    Reads binary spectral data from a .dat file.
    """
    dat_path = specfile_path_base + ".dat"
    if not path.exists(dat_path):
        print(f"Warning: Data file not found: {dat_path}")
        return np.array([])
        
    with open(dat_path, "rb") as f:
        spec = np.fromfile(f, np.float32)
        
    npo = int(header.get("npo", 0))
    if npo > 0 and len(spec) != npo:
        print(f'Warning: No. of points in {path.basename(dat_path)} does not match header: npo = {npo}, length = {len(spec)}')
    
    return spec

def read_linelist(specfile_path_base: str) -> List[Dict[str, Any]]:
    """
    Reads a binary .lin file from Xgremlin into a list of dictionaries.
    """
    linel = []
    lin_path = specfile_path_base + ".lin"
    if not path.exists(lin_path):
        return linel # This is a common case, so no warning is printed.
        
    with open(lin_path, "rb") as flin:
        # Read the header information
        try:
            nlin = unpack("i", flin.read(4))[0]
            flin.read(4)  # Skip linlen
            flin.read(312) # Skip rest of header

            # Read in all the lines
            for _ in range(nlin):
                sp = {}
                sp['sig'], sp['xint'], sp['width'], sp['dmping'], sp['itn'], sp['ihold'] = unpack("dfffhh", flin.read(24))
                sp['tags'] = flin.read(4)
                sp['epstot'], sp['epsevn'], sp['epsodd'], sp['epsran'], sp['spare'] = unpack("fffff", flin.read(20))
                # Decode and clean the identifier string
                sp['ident'] = flin.read(32).decode('utf-8', errors='ignore').strip('\x00')
                linel.append(sp)
        except Exception as e:
            print(f"Error reading linelist file {lin_path}: {e}")

    return linel

def read_lifetimes(lifetime_files_glob: str) -> Tuple[Dict[str, float], Dict[str, str], List[str]]:
    """
    Imports level lifetimes from .lev files matching a glob pattern.
    """
    lifetimes = {}
    upper_values = {}
    levels = []
    lifetime_files = glob(lifetime_files_glob)

    for lifetime_file in lifetime_files:
        with open(lifetime_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) < 7:
                    continue
                upper_level = parts[6]
                upper_values[upper_level] = parts[2]
                try:
                    lifetime = float(parts[5])
                    levels.append(upper_level)
                except (ValueError, IndexError):
                    lifetime = None
                lifetimes[upper_level] = lifetime
    
    return lifetimes, upper_values, levels

def read_theoretical_A_values(calc_file: str) -> Dict[Tuple[str, str], float]:
    """
    Reads calculated theoretical A-values from a file.
    """
    E1 = {}
    if not path.exists(calc_file):
        print(f"Warning: Theoretical calculation file not found: {calc_file}")
        return E1

    with open(calc_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            waveno_str = parts[0]
            upper_level = parts[-1]
            try:
                E1[(upper_level, waveno_str)] = float(parts[1])
            except ValueError:
                continue
    return E1

def read_identified_lines(id_lines_file: str, energy_levels_file: str, upper_level: str) -> List[str]:
    """
    Finds all previously identified lines from a given upper level.
    """
    upper_energy_key = ""
    try:
        with open(energy_levels_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) > 6 and parts[6] == upper_level:
                    upper_energy_key = parts[6]
                    break
    except FileNotFoundError:
        print(f"Warning: Energy levels file not found: {energy_levels_file}")
        return []

    if not upper_energy_key:
        return []
    return grep_open(upper_energy_key, id_lines_file)

def grep_open(grep_key: str, open_file: str) -> List[str]:
    """
    Opens a file and returns all lines containing the grep_key.
    """
    try:
        with open(open_file, 'r') as f:
            read_lines = f.readlines()
        return [line for line in read_lines if grep_key in line]
    except FileNotFoundError:
        print(f"Warning: Grep target file not found: {open_file}")
        return []
