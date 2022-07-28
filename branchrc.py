# This file is the resource file for the python program 'branch.py'
# All entries must be valid python code

def set_parameters():

    upper_level = "a3F4p_y4F4*"
    file_name = "y4F4"
    reference_level = "3F4s__a4F4"
    transfer_file = "Cr102700.03.II"
    transfer_level = "d5____a4G5"
	
    identified_lines_file = "CrII.CS"
    lifetime_file = "cr_ii.lev"
    delfile = "Cr102700.003.II"
    delval = "d5____a4P1"
    lifetime_unc = 0.10
    all_spectrum_files = '*.II'
    calc_file = "CrII_waveno.E1"
    discrim = 0.2            # Wavelength discriminator to match observations with calculations
    unc_cal = 0.07            # Calibration uncertainty estimated from 5% discrepancy in measurements
                             # coverts to 2x5/sqrt(3) = 5.8 %. Add to lamp uncertainty of 3.5% = 6.8%
                             # Round up to 7 %

    return (file_name, upper_level, reference_level, transfer_file,
			transfer_level, identified_lines_file, lifetime_file, lifetime_unc, all_spectrum_files,
            delfile,delval,calc_file,discrim, unc_cal)
