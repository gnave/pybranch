# pybranch/visualization.py
"""Handles all plotting and visualization for the PyBranch project."""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
from typing import Optional

from .data_io import Spectrum
# --- FIX 1: Import the analysis module to access its functions ---
from . import analysis

def plot_spectrum_line(
    spectrum_obj: Spectrum,
    wavenumber: float,
    window_length: int,
    wavcorr_applied: float = 0.0
) -> Optional[Figure]:
    """
    Generates a Matplotlib figure of a specific spectral line from a Spectrum object.

    Args:
        spectrum_obj: The loaded Spectrum object containing the data.
        wavenumber: The central wavenumber of the line to plot (in corrected scale).
        window_length: The number of data points to include in the plot window.
        wavcorr_applied: The wavenumber correction factor that was applied to the
            experimental data's x-axis.

    Returns:
        A Matplotlib Figure object if successful, otherwise None.
    """
    if spectrum_obj.wavenumbers is None or spectrum_obj.intensities is None:
        print("Error: Spectrum data is not loaded.")
        return None

    center_index = np.argmin(np.abs(spectrum_obj.wavenumbers - wavenumber))
    start_index = max(0, center_index - window_length // 2)
    end_index = min(len(spectrum_obj.wavenumbers), start_index + window_length)

    # The x_data is already corrected from when the Spectrum object was loaded
    x_data = spectrum_obj.wavenumbers[start_index:end_index]
    y_data = spectrum_obj.intensities[start_index:end_index]

    fig, ax = plt.subplots()
    fig.suptitle(f"Line near {wavenumber:.3f} cm⁻¹ in {os.path.basename(spectrum_obj.path_base)}")
    ax.set_xlabel('Wavenumber / cm⁻¹')
    ax.set_ylabel('Intensity / arb. units')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%9.3f'))
    ax.plot(x_data, y_data, label='Experimental Data')

    # --- THIS IS THE CORRECTION LOGIC ---
    # Temporarily "de-correct" the wavenumber to match the uncorrected .lin file data
    wavenumber_to_fit_uncorrected = wavenumber
    if wavcorr_applied != 0:
        # Use division, as it's the exact inverse of original = corrected * (1 + corr)
        wavenumber_to_fit_uncorrected = wavenumber / (1 + wavcorr_applied)

    # Perform the fit in the uncorrected coordinate system
    fit_result = analysis.fit_voigt_profile(wavenumber_to_fit_uncorrected, x_data, spectrum_obj.linelist)
    
    if fit_result:
        # The returned fit_x is initially uncorrected
        fit_x_uncorrected, fit_y, linefit = fit_result
        
        # "Re-correct" the fit's x-axis so it aligns with the corrected data for plotting
        fit_x_corrected = fit_x_uncorrected
        if wavcorr_applied != 0:
            fit_x_corrected = fit_x_uncorrected * (1 + wavcorr_applied)

        ax.plot(fit_x_corrected, fit_y, ls='dotted', label='Voigt Fit')
        
        # Also correct the line center from the fit for display purposes
        corrected_sig = linefit['sig'] * (1 + wavcorr_applied)
        damp = (linefit['dmping'] - 1) / 25
        fwhm = linefit['width'] / 1000.
        params = f"Fit: σ = {corrected_sig:.3f}, Int = {linefit['xint']:.0f}, FWHM = {fwhm:.3f}, Damp = {damp:.3f}"
        fig.text(0.5, 0.02, params, ha='center', fontsize=8)
    # --- END OF CORRECTION LOGIC ---
    
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    return fig