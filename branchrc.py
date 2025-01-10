# This file is the resource file for the python program 'branch.py'
# All entries must be valid python code

def set_parameters():

    file_name = "y4F4"                     # 0
    upper_level = "a3F4p_y4F4*"            # 1
    reference_level = "3F4s__a4F4"         # 2
    reference_file = "Cr102700.003.I"      # 3
    normal_level = "d5____a4G5"            # 4
    identified_lines_file = "CrII.CS"      # 5
    level_file = "cr_ii.lev"               # 6
    lifetime_unc = 0.10                    # 7
    spectrum_files = '*.II'                # 8
    calc_file = "CrII_waveno.E1"           # 9
    discrim = 0.2            # 10 : Wavenumber discriminator to match observations with calculations
    plotwin_length = 32                    # 11        


    return (file_name, upper_level, reference_level, reference_file,
			normal_level, identified_lines_file, level_file, lifetime_unc, spectrum_files,
            calc_file,discrim,plotwin_length)
