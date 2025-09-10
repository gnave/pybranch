# pybranch/analysis.py

import math
import numpy as np
from scipy.special import voigt_profile
from typing import List, Dict, Any, Tuple, Optional
import os

# Import the data reading functions that the analysis will depend on
from . import data_io

# --- Main Analysis Class ---

class BranchingRatioAnalysis:
    """Manages the end-to-end analysis of branching fractions for an upper level."""
    # (Docstrings are unchanged)

    def __init__(self, upper_level: str, config: Dict[str, Any]):
        """Initializes the analysis for a specific upper level."""
        # This is the correct, simple initialization.
        # self.upper_level will contain the asterisk, which is the correct key.
        self.upper_level: str = upper_level
        
        self.config: Dict[str, Any] = config
        self.lifetimes: Dict[str, float] = {}
        self.theoretical_A: Dict[Tuple[str, str], float] = {}
        self.known_lines: List[str] = []
        self.snrs: Dict[Tuple[str, str], float] = {}
        self.intensities: Dict[Tuple[str, str], float] = {}
        self.uncertainties: Dict[Tuple[str, str], float] = {}
        self.transition_ids: Dict[float, str] = {}
        self.wavenumbers: List[float] = []
        self.results: Optional[Dict[str, Any]] = None

    def run_full_analysis(self) -> None:
        """Executes all the major steps of the analysis in sequence."""
        print("1. Loading prerequisite atomic data...")
        self.load_atomic_data()
        print("2. Aggregating data from observed experimental files...")
        self.aggregate_observed_data()
        print("3. Calculating final branching fractions and A-values...")
        self.calculate_branching_ratios()
        print(f"Analysis complete for upper level: {self.upper_level}")

    def load_atomic_data(self) -> None:
        """Loads all necessary atomic data from files specified in the config."""
        conf_files = self.config.get('files', {})
        self.lifetimes, _, _ = data_io.read_lifetimes(conf_files.get('levels_glob', '*.lev'))
        self.theoretical_A = data_io.read_theoretical_A_values(conf_files.get('theoretical_lines', ''))
        master_lev_file = conf_files.get('master_level_file')
        if not master_lev_file or not os.path.exists(master_lev_file):
            raise FileNotFoundError(
                f"Configuration Error: 'master_level_file' ({master_lev_file}) not found or not specified in config.yaml."
            )
        self.known_lines = data_io.read_identified_lines(
            conf_files.get('identified_lines', ''),
            master_lev_file,
            self.upper_level
        )

    def aggregate_observed_data(self) -> None:
        """Aggregates and processes data from all experimental summary files."""
        files_config = self.config.get('files', {})
        unc_config = self.config.get('uncertainties', {})
        spectrum_files_glob = files_config.get('spectrum_glob', '*.II')
        all_spectra_files = data_io.glob(spectrum_files_glob)
        master_lev_file = files_config.get('master_level_file')
        if not master_lev_file or not os.path.exists(master_lev_file):
             raise FileNotFoundError(
                f"Configuration Error: 'master_level_file' ({master_lev_file}) not found or not specified in config.yaml."
            )
        with open(master_lev_file) as f:
            energy_levels = f.readlines()
        w_maxI = {}
        calunc_per_1000 = {}
        for spectrum_file in all_spectra_files:
            with open(spectrum_file) as f:
                params = f.readline().split()
            if len(params) < 3: continue
            resoln, lowE, hiE = float(params[0]), float(params[1]), float(params[2])
            cal_factor = unc_config.get('calibration_factor', 70.0)
            calunc_per_1000[spectrum_file] = cal_factor / (hiE - lowE) if (hiE - lowE) != 0 else 0
            
            transitions = data_io.grep_open("- " + self.upper_level, spectrum_file)
            
            w_maxI[spectrum_file], maxI = 0, 0
            for t in transitions:
                parts = t.split()
                if len(parts) < 7: continue
                # Stripping the asterisk from the LOWER level is correct, as per original code.
                lower_level, intensity = parts[6].strip('*'), float(parts[1])
                self.snrs[spectrum_file, lower_level] = float(parts[0])
                self.intensities[spectrum_file, lower_level] = intensity
                if intensity > maxI:
                    maxI, w_maxI[spectrum_file] = intensity, float(parts[4])
                root_npts = (0.001 * float(parts[2]) / resoln) if resoln > 0 else 0
                snr_val = self.snrs.get((spectrum_file, lower_level), 0)
                variance = 2.25 / (snr_val**2 * root_npts) if snr_val > 0 and root_npts > 0 else float('inf')
                self.uncertainties[spectrum_file, lower_level] = math.sqrt(variance)
                for level in energy_levels:
                    level_parts = level.split()
                    if len(level_parts) > 6 and level_parts[6] == lower_level:
                        for line in self.known_lines:
                            line_parts = line.split()
                            if len(line_parts) > 3 and line_parts[3].strip('*') == level_parts[6].strip('*'):
                                self.transition_ids[float(line_parts[0])] = lower_level
        self.wavenumbers = sorted(list(set(self.transition_ids.keys())), reverse=True)
        for (spectrum, lower_level), unc in self.uncertainties.items():
            wnum = next((k for k, v in self.transition_ids.items() if v == lower_level), 0)
            calunc = calunc_per_1000.get(spectrum, 0) * (wnum - w_maxI.get(spectrum, 0)) / 1000.0
            if unc != float('inf'):
                self.uncertainties[spectrum, lower_level] = math.sqrt(unc**2 + calunc**2)

    def calculate_branching_ratios(self) -> Dict[str, Any]:
        """Performs the final calculation and stores the results."""
        analysis_params = self.config.get('analysis', {})
        unc_params = self.config.get('uncertainties', {})
        self.results = calculate_final_results(
            upper_level=self.upper_level, # Pass the key with the asterisk
            wavenumbers=self.wavenumbers,
            transition_ids=self.transition_ids,
            intensities=self.intensities,
            uncertainties=self.uncertainties,
            lifetimes=self.lifetimes,
            lifetime_unc=unc_params.get('lifetime', 0.10),
            theoretical_A_values=self.theoretical_A,
            discrim=analysis_params.get('wavenumber_discriminator', 0.2)
        )
        return self.results
    
    def normalize_by_reference_line(self, reference_level: str) -> None:
        """
        Rescales the intensities within each individual spectrum.

        For each spectrum file, the intensities of all lines are rescaled
        such that the specified reference line has an intensity of 1000. This
        puts all spectra on a common intensity scale before averaging.

        Args:
            reference_level: The lower level key of the line to use
                             for normalization.
        """
        if not self.intensities:
            print("Warning: Cannot normalize. Intensity data has not been loaded yet.")
            return

        self.intensities = normalize_intensities_by_reference_line(
            self.intensities,
            reference_level
        )
        print(f"Intensities for all spectra normalized using '{reference_level}' as the reference (value = 1000).")

    def rescale_spectrum_by_transfer_line(
        self,
        transfer_level: str,
        reference_file: str,
        initial_reference_level: str
    ) -> None:
        """
        Rescales a single spectrum against the weighted average of all others.
        """
        if not self.intensities or not self.uncertainties:
            print("Warning: Cannot rescale. Intensity data has not been loaded yet.")
            return

        # --- NEW ROBUSTNESS CHECK ---
        # Get a list of all unique spectrum filenames from the data keys
        all_spectra_keys = sorted(list(set(key[0] for key in self.intensities.keys())))
        
        correct_reference_file = None
        # First, check for an exact match
        if reference_file in all_spectra_keys:
            correct_reference_file = reference_file
        else:
            # If no exact match, try to find a key that starts with the user's input.
            # This handles cases where the user provides "file" instead of "file.II".
            matches = [key for key in all_spectra_keys if os.path.basename(key).startswith(reference_file)]
            if len(matches) == 1:
                correct_reference_file = matches[0]
                print(f"Info: Matched partial filename '{reference_file}' to full key '{correct_reference_file}'.")
            elif len(matches) > 1:
                print(f"Warning: Ambiguous reference file '{reference_file}'. Found multiple matches: {matches}. No rescaling performed.")
                return
            else:
                print(f"Warning: Reference file '{reference_file}' not found in the loaded dataset. No rescaling performed.")
                return
        # --- END OF CHECK ---

        weighting_cap = self.config.get('analysis', {}).get('weighting_cap', 555)
        
        self.intensities, self.uncertainties = rescale_intensities_by_transfer_line(
            self.intensities,
            self.uncertainties,
            transfer_level,
            correct_reference_file, # <-- Use the validated, correct key
            initial_reference_level,
            weighting_cap
        )
        print(f"Spectrum '{correct_reference_file}' rescaled using '{transfer_level}' as the transfer line.")

    def log_intensity_table(self, title: str) -> None:
        """
        Prints a formatted table of the current SNRs and intensities.

        This is a debugging and logging tool, similar to the original program's
        'show_values' function, to inspect the state of the data at any
        point in the analysis pipeline.

        Args:
            title: A title to print above the table to identify the context.
        """
        if not self.intensities:
            print(f"\n--- {title} ---")
            print("No intensity data loaded to display.")
            return

        all_spectra = sorted(list(set(key[0] for key in self.intensities.keys())))
        file_num = len(all_spectra)
        width = 46 + 18 * file_num + 18 # Add space for mean/stdev

        print("\n" + "="*80)
        print(f"--- {title} ---")
        print("="*80)
        
        # --- Print Header ---
        print('-'*width)
        file_name_line = f"| {'File Name:':<23s} | "
        for spectrum in all_spectra:
            # Shorten the filename for display
            display_name = os.path.basename(spectrum)
            file_name_line += f'{display_name:>14s}  | '
        file_name_line += f'{"Mean":>7s} | {"Stdev":>7s} |'
        print(file_name_line)
        print('-'*width)
        print(f"| {'Wavenumber':<10s} | {'L Level':<10s} | " + f'{"SNR":>5s} | {"Intensity":>8s} | '*file_num + f'{"":>7s} | {"(frac)":>7s} |')
        print('-'*width)

        # --- Print Data Rows ---
        for wavenumber in self.wavenumbers:
            lower_level = self.transition_ids[wavenumber]
            data_line = f"| {wavenumber:>10.3f} | {lower_level:<10s} | "

            # Get data for just this line
            values_for_level = {s: i for (s, l), i in self.intensities.items() if l == lower_level}
            uncs_for_level = {s: u for (s, l), u in self.uncertainties.items() if l == lower_level}
            
            # Calculate the current weighted mean
            mean_result = calculate_weighted_mean(values_for_level, uncs_for_level)
            
            for spectrum in all_spectra:
                snr = self.snrs.get((spectrum, lower_level))
                intensity = self.intensities.get((spectrum, lower_level))
                
                if snr is not None and intensity is not None:
                    data_line += f"{round(snr):>5d} | {round(intensity):>8d} | "
                else:
                    data_line += f"  --- |      --- | "

            if mean_result:
                mean_I, stdev_frac = mean_result
                data_line += f"{round(mean_I):>7d} | {stdev_frac:>7.3f} |"
            else:
                data_line += f"  --- |     --- |"
            
            print(data_line)
            
        print('-'*width + "\n")

    def delete_line(self, spectrum_file: str, lower_level: str) -> None:
        """
        Deletes a single data point (intensity, SNR, uncertainty) for a
        specific line in a specific spectrum.

        Args:
            spectrum_file: The full filename of the spectrum to modify.
            lower_level: The lower level key for the line to be deleted.
        """
        key = (spectrum_file, lower_level)

        # Use .pop() with a default value of None. This safely removes the key
        # if it exists and does nothing if it's already gone, preventing errors.
        intensity_popped = self.intensities.pop(key, None)
        self.snrs.pop(key, None)
        self.uncertainties.pop(key, None)

        if intensity_popped is not None:
            print(f"Info: Deleted line '{lower_level}' from spectrum '{spectrum_file}'.")
        else:
            print(f"Warning: Could not delete line '{lower_level}' from '{spectrum_file}' (was not found).")

    def remove_transition_entirely(self, wavenumber_to_remove: float) -> None:
        """
        Removes a transition completely from the analysis.

        This is a more powerful deletion tool. It removes all associated data
        points (SNR, intensity) for a given transition AND removes the transition
        itself from the lists of wavenumbers and transition IDs. This ensures
        the line is treated as "unobserved" in the final calculation.

        Args:
            wavenumber_to_remove: The precise wavenumber of the transition to remove.
        """
        # 1. Find the corresponding lower level key
        lower_level_to_remove = self.transition_ids.get(wavenumber_to_remove)

        if not lower_level_to_remove:
            print(f"Warning: Wavenumber {wavenumber_to_remove} not found in transition list. Cannot remove.")
            return

        # 2. Remove all data points associated with this lower level
        # Create a list of keys to pop to avoid modifying dict while iterating
        keys_to_pop = [key for key in self.intensities if key[1] == lower_level_to_remove]
        for key in keys_to_pop:
            self.intensities.pop(key, None)
            self.snrs.pop(key, None)
            self.uncertainties.pop(key, None)

        # 3. Remove the transition from the core analysis structure
        self.transition_ids.pop(wavenumber_to_remove, None)
        if wavenumber_to_remove in self.wavenumbers:
            self.wavenumbers.remove(wavenumber_to_remove)
        
        print(f"Info: Completely removed transition at {wavenumber_to_remove} ({lower_level_to_remove}) from analysis.")

 

# --- Pure, Stateless Scientific Functions ---

def calc_gaussian_width_from_voigt(voigt_width: float, damping: float) -> float:
    """
    Calculates the Gaussian width of a Voigt profile from the total Voigt width
    and the damping parameter, using the approximation from Kielkopf (1973).

    Args:
        voigt_width: The Full Width at Half Maximum (FWHM) of the Voigt profile.
        damping: The damping parameter (ratio of Lorentzian to Gaussian width).

    Returns:
        The FWHM of the Gaussian component.
    """
    # This is equation 8 from Kielkopf, JOSA 63, 987 (1973)
    eta = 0.099
    A = 1 + eta * np.log(2)
    B = eta * np.log(2)

    # Ensure the term inside the square root is non-negative
    inner_term = 1 - A * damping + B * damping * damping
    if inner_term < 0:
        return 0.0  # Or handle as an error
    
    gauss_width = voigt_width * np.sqrt(inner_term)
    return gauss_width



def fit_voigt_profile(
    wavenumber_to_fit: float,
    x_grid: np.ndarray,
    linelist: List[Dict[str, Any]],
    wavenumber_window: float = 0.1
) -> Optional[Tuple[np.ndarray, np.ndarray, Dict[str, Any]]]:
    """
    Calculates a Voigt profile for a line if it exists in a given linelist.

    Args:
        wavenumber_to_fit: The central wavenumber of the line to fit.
        x_grid: The numpy array of wavenumbers to calculate the profile on.
        linelist: A list of line parameter dictionaries (from data_io.read_linelist).
        wavenumber_window: The tolerance for matching the wavenumber.

    Returns:
        A tuple containing (x_grid, y_profile, line_parameters) if a match is found,
        otherwise None.
    """
    linefit = None
    for line in linelist:
        if abs(line['sig'] - wavenumber_to_fit) < wavenumber_window:
            linefit = line
            break

    if linefit is None:
        return None

    # Xgremlin widths are in 0.001 cm-1, so convert to cm-1.
    # voigt_profile uses HWHM, not FWHM, so divide by 2.
    width_hwhm = linefit['width'] / 2000.0

    # Xgremlin damping parameter runs from 1 to 26, convert to a 0-1 scale.
    damping = (linefit['dmping'] - 1) / 25.0

    # Calculate Gaussian width (FWHM) and convert to standard deviation for voigt_profile.
    gauss_fwhm = calc_gaussian_width_from_voigt(width_hwhm * 2, damping)
    gauss_sigma = gauss_fwhm / (2 * np.sqrt(2 * np.log(2)))
    
    # Lorentzian width (HWHM)
    lorentz_hwhm = width_hwhm * damping
   
    # Calculate the profile and normalize by the peak intensity.
    y = voigt_profile(x_grid - linefit['sig'], gauss_sigma, lorentz_hwhm)
    peak_value = voigt_profile(0, gauss_sigma, lorentz_hwhm)
    y_scaled = (y / peak_value) * linefit['xint']

    return x_grid, y_scaled, linefit

def calculate_weighted_mean(
    values: Dict[str, float],
    uncertainties: Dict[str, float]
) -> Optional[Tuple[float, float]]:
    # (Implementation is unchanged)
    avg_int, sum_weight, max_weight = 0.0, 0.0, 555
    for key, value in values.items():
        if key in uncertainties and uncertainties[key] > 0 and uncertainties[key] != float('inf'):
            weight = 1 / (uncertainties[key] ** 2)
            if weight > max_weight: weight = max_weight
            avg_int += value * weight
            sum_weight += weight
    if sum_weight > 0:
        return (avg_int / sum_weight), (1 / math.sqrt(sum_weight))
    return None

def calculate_final_results(
    upper_level: str,
    wavenumbers: List[float],
    transition_ids: Dict[float, str],
    intensities: Dict[Tuple[str, str], float],
    uncertainties: Dict[Tuple[str, str], float],
    lifetimes: Dict[str, float],
    lifetime_unc: float,
    theoretical_A_values: Dict[Tuple[str, str], float],
    discrim: float
) -> Dict[str, Any]:
    """Calculates the final branching fractions, transition probabilities, and uncertainties."""
    
    results_output = {'summary': {}, 'table': [], 'residuals_log': []}
    level_keys = list(set(transition_ids.values()))
    mean_intensities, stdevs, total_intensity = {}, {}, 0
    for lk in level_keys:
        vals = {s: i for (s, l), i in intensities.items() if l == lk}
        uncs = {s: u for (s, l), u in uncertainties.items() if l == lk}
        mean_res = calculate_weighted_mean(vals, uncs)
        if mean_res: 
            mean_intensities[lk], stdevs[lk] = mean_res
            total_intensity += mean_intensities[lk]
        else: 
            mean_intensities[lk], stdevs[lk] = 0, 0
    
    sum_unobserved_A, theoval = 0, {}
    for (ul_from_file, wno_str), A_val in theoretical_A_values.items():
        if ul_from_file == upper_level:
            is_obs = any(abs(obs_wno - float(wno_str)) < discrim for obs_wno in wavenumbers)
            if is_obs:
                 obs_wno = next(w for w in wavenumbers if abs(w - float(wno_str)) < discrim)
                 theoval[(ul_from_file, obs_wno)] = A_val
            else:
                sum_unobserved_A += A_val
                results_output['residuals_log'].append(f'{wno_str:>10s} | {A_val:>7.4f}')

    upper_level_lifetime = lifetimes.get(upper_level, 0)
    frac_resid = (sum_unobserved_A * upper_level_lifetime / 1000.0) if upper_level_lifetime > 0 else 0
    total_intensity *= (1 + frac_resid)
    
    results_output['summary'] = {'level': upper_level, 'lifetime': upper_level_lifetime, 'residual_percent': frac_resid * 100.0}

    BFsq = sum((mean_intensities[lk] / total_intensity)**2 * stdevs[lk]**2 for lk in level_keys if total_intensity > 0 and lk in mean_intensities and lk in stdevs)
    BFsq += (frac_resid**2) * 0.25

    for wavenumber in wavenumbers:
        level_key = transition_ids[wavenumber]
        mean_I = mean_intensities.get(level_key, 0)
        
        branching_fraction_percent = (mean_I / total_intensity * 100) if total_intensity > 0 else 0
        
        # --- THIS IS THE CORRECTED FORMULA ---
        if upper_level_lifetime > 0:
            # A [10^8 s^-1] = BF_percent / (10 * tau_ns)
            A_val = branching_fraction_percent / (10 * upper_level_lifetime)
        else:
            A_val = 0
        # --- END OF CORRECTION ---
        
        variance = (stdevs.get(level_key, 0)**2) * (1 - 2 * (branching_fraction_percent/100)) + BFsq
        bf_unc_pct = 100 * math.sqrt(variance) if variance > 0 else 0
        aval_unc_pct = 100 * math.sqrt(variance + lifetime_unc**2) if (variance + lifetime_unc**2) > 0 else 0
        
        results_output['table'].append({
            'wavenumber': wavenumber, 'level_key': level_key,
            'branching_fraction_percent': branching_fraction_percent,
            'bf_uncertainty_percent': bf_unc_pct, 'A_value': A_val,
            'A_value_uncertainty_percent': aval_unc_pct,
            'theoretical_A_value': theoval.get((upper_level, wavenumber), 0)
        })
    return results_output

def normalize_intensities_by_reference_line(
    intensities: Dict[Tuple[str, str], float],
    reference_level: str
) -> Dict[Tuple[str, str], float]:
    """
    Normalizes each spectrum's intensities so the reference line is 1000.

    Args:
        intensities: The dictionary of intensity data, keyed by (spectrum_file, lower_level).
        reference_level: The lower level key of the line to normalize to.

    Returns:
        A new dictionary with the normalized intensities.
    """
    normalized_intensities = intensities.copy()
    # Get a unique list of all spectra present in the data
    all_spectra = sorted(list(set(key[0] for key in intensities.keys())))

    for spectrum_file in all_spectra:
        reference_key = (spectrum_file, reference_level)
        
        # Check if the reference line exists in this specific spectrum
        if reference_key in intensities:
            norm_factor = intensities[reference_key]
            
            if norm_factor > 0:
                # Find all lines belonging to this spectrum and rescale them
                for (s, lower_level), intensity in intensities.items():
                    if s == spectrum_file:
                        normalized_intensities[(s, lower_level)] = 1000.0 * intensity / norm_factor
        else:
            print(f"Warning: Reference level '{reference_level}' not found in spectrum '{spectrum_file}'. Skipping normalization for this file.")
            
    return normalized_intensities

def rescale_intensities_by_transfer_line(
    intensities: Dict[Tuple[str, str], float],
    uncertainties: Dict[Tuple[str, str], float],
    transfer_level: str,
    reference_file: str,
    initial_reference_level: str,
    weighting_cap: int = 555
) -> Tuple[Dict[Tuple[str, str], float], Dict[Tuple[str, str], float]]:
    """
    Rescales a single spectrum using a weighted average from all other spectra.

    Args:
        intensities: Dictionary of current intensity data.
        uncertainties: Dictionary of current fractional uncertainty data.
        transfer_level: The key for the transfer line used for rescaling.
        reference_file: The path of the single spectrum to be rescaled.
        initial_reference_level: The original reference level, used to filter
            which spectra contribute to the average.
        weighting_cap: The cap on statistical weights.

    Returns:
        A tuple containing the new, rescaled intensities and uncertainties dictionaries.
    """
    new_intensities = intensities.copy()
    new_uncertainties = uncertainties.copy()

    # 1. Calculate the weighted average intensity of the transfer line
    #    across all OTHER spectra.
    values_for_transfer_level = {}
    uncs_for_transfer_level = {}
    for (spectrum, lower_level), intensity in intensities.items():
        # --- THIS IS THE CRITICAL FIX ---
        # Exclude the file we are trying to rescale from the average calculation.
        if spectrum == reference_file:
            continue
        # --- END OF FIX ---

        # Per original logic, only include spectra that have BOTH the initial
        # reference level and the new transfer level.
        if lower_level == transfer_level and (spectrum, initial_reference_level) in intensities:
            values_for_transfer_level[spectrum] = intensity
            uncs_for_transfer_level[spectrum] = uncertainties[(spectrum, lower_level)]

    mean_result = calculate_weighted_mean(values_for_transfer_level, uncs_for_transfer_level)
    if not mean_result:
        print(f"Warning: Could not calculate weighted average for transfer level '{transfer_level}'. No valid other spectra found. No rescaling performed.")
        return new_intensities, new_uncertainties
    
    avg_int, fractional_unc_of_avg = mean_result

    # 2. Calculate the rescaling factor
    intensity_in_ref_file = intensities.get((reference_file, transfer_level))
    if not intensity_in_ref_file or intensity_in_ref_file == 0:
        print(f"Warning: Transfer level '{transfer_level}' has zero intensity in reference file '{reference_file}'. No rescaling performed.")
        return new_intensities, new_uncertainties
        
    rescale_factor = avg_int / intensity_in_ref_file

    # 3. Apply the factor to all lines within the reference file and propagate uncertainty
    for (spectrum, lower_level), intensity in intensities.items():
        if spectrum == reference_file:
            # Rescale the intensity
            new_intensities[(spectrum, lower_level)] = intensity * rescale_factor

            # Propagate the uncertainty
            old_frac_unc_sq = uncertainties.get((spectrum, lower_level), 0)**2
            frac_unc_of_factor_sq = fractional_unc_of_avg**2
            
            new_uncertainties[(spectrum, lower_level)] = math.sqrt(old_frac_unc_sq + frac_unc_of_factor_sq)

    return new_intensities, new_uncertainties

def format_intensity_table_as_string(analysis_obj: 'BranchingRatioAnalysis', title: str) -> str:
    """
    Formats the current SNRs and intensities into a multi-line string table.

    Args:
        analysis_obj: The active BranchingRatioAnalysis instance.
        title: A title for the table.

    Returns:
        A formatted string representing the data table for display or logging.
    """
    if not analysis_obj.intensities:
        return f"\n--- {title} ---\nNo intensity data loaded to display."

    lines = []
    all_spectra = sorted(list(set(key[0] for key in analysis_obj.intensities.keys())))
    file_num = len(all_spectra)
    width = 46 + 18 * file_num + 18  # Add space for mean/stdev

    lines.append("\n" + "="*80)
    lines.append(f"--- {title} ---")
    lines.append("="*80)
    
    # --- Header ---
    lines.append('-'*width)
    file_name_line = f"| {'File Name:':<23s} | "
    for spectrum in all_spectra:
        display_name = os.path.basename(spectrum)
        file_name_line += f'{display_name:>14s}  | '
    file_name_line += f'{"Mean":>7s} | {"Stdev":>7s} |'
    lines.append(file_name_line)
    lines.append('-'*width)
    lines.append(f"| {'Wavenumber':<10s} | {'L Level':<10s} | " + f'{"SNR":>5s} | {"Intensity":>8s} | '*file_num + f'{"":>7s} | {"(frac)":>7s} |')
    lines.append('-'*width)

    # --- Data Rows ---
    for wavenumber in analysis_obj.wavenumbers:
        lower_level = analysis_obj.transition_ids[wavenumber]
        data_line = f"| {wavenumber:>10.3f} | {lower_level:<10s} | "

        values_for_level = {s: i for (s, l), i in analysis_obj.intensities.items() if l == lower_level}
        uncs_for_level = {s: u for (s, l), u in analysis_obj.uncertainties.items() if l == lower_level}
        mean_result = calculate_weighted_mean(values_for_level, uncs_for_level)
        
        for spectrum in all_spectra:
            snr = analysis_obj.snrs.get((spectrum, lower_level))
            intensity = analysis_obj.intensities.get((spectrum, lower_level))
            if snr is not None and intensity is not None:
                data_line += f"{round(snr):>5d} | {round(intensity):>8d} | "
            else:
                data_line += f"  --- |      --- | "

        if mean_result:
            mean_I, stdev_frac = mean_result
            data_line += f"{round(mean_I):>7d} | {stdev_frac:>7.3f} |"
        else:
            data_line += f"  --- |     --- |"
        lines.append(data_line)
            
    lines.append('-'*width + "\n")
    return "\n".join(lines)

def format_final_results_as_string(results: Dict[str, Any]) -> str:
    """
    Formats the final results dictionary into a multi-line string for display.

    Args:
        results: The final results dictionary from calculate_final_results.

    Returns:
        A formatted string of the summary, results table, and residuals.
    """
    lines = []
    if not results or not results.get('table'):
        return "\nAnalysis finished, but no results were generated."

    lines.append("\n" + "="*110)
    lines.append("                              Branching Fraction Results")
    lines.append("="*110)
    
    summary = results['summary']
    lines.append(f"Level: {summary.get('level', 'N/A'):>10s} | Lifetime: {summary.get('lifetime', 0.0):.3f} ns | Unobserved Residual: {summary.get('residual_percent', 0.0):>6.3f} %")
    
    lines.append("-"*110)
    lines.append(f"| {'Wavenumber':^12s} | {'Lower Level':^12s} | {'Branching Fr.':^15s} | {'BF Unc.':^12s} | {'Trans. Prob.':^15s} | {'A-Val Unc.':^12s} | {'Theoretical A':^15s} |")
    lines.append(f"| {'(cm-1)':^12s} | {'':^12s} | {'(%)':^15s} | {'(%)':^12s} | {'(10^8 s^-1)':^15s} | {'(%)':^12s} | {'(10^8 s^-1)':^15s} |")
    lines.append("-"*110)

    for row in results['table']:
        aval_scaled = row['A_value']
        theo_aval_scaled = row['theoretical_A_value']
        lines.append(f"| {row['wavenumber']:>12.3f} | {row['level_key']:<10s} | {row['branching_fraction_percent']:>15.3f} | {row['bf_uncertainty_percent']:>12.1f} | {aval_scaled:>15.4f} | {row['A_value_uncertainty_percent']:>12.1f} | {theo_aval_scaled:>15.4f} |")
    lines.append("-"*110)

    if results.get('residuals_log'):
        lines.append("\nUnobserved transitions contributing to residual:")
        lines.append("---------------------------------------------")
        for line in results['residuals_log']:
            lines.append(line)
            
    return "\n".join(lines)