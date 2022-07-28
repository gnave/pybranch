#  
#  A script for calculating branching fractions from sets of intensity ratios
#  measured in various emission spectra.
#  
import sys
import math
from branchgrep import grep_open
from branchrc import set_parameters
from glob import glob
from tkinter import *
from tkinter.scrolledtext import ScrolledText

file_name= set_parameters()[0]
upper_level= set_parameters()[1]
reference_level= set_parameters()[2]
transfer_file= set_parameters()[3]
normal_level= set_parameters()[4]
identified_lines_file= set_parameters()[5]
lifetime_file= set_parameters()[6]
lifetime_unc= set_parameters()[7]
all_spectrum_files= glob(set_parameters()[8])
delfile= set_parameters()[9]
delval= set_parameters()[10]
calc_file= set_parameters()[11]
discrim = set_parameters()[12]
unc_cal = set_parameters()[13]


def get_calculations(E1):
#
# This function gets the calculated A values for calculating a residual
#

    for line in open(calc_file, 'r').readlines():
        line = line.split()
        waveno = line[0]             # Wavenumber must be first item in list
        upper_level = line[-1]       # Key is the last item in the list
        E1[upper_level,waveno] = float( line[1])   # Aval is 2nd item, units 1e6 s-1.
    return(E1)
     
class BranchingFractionCalc(Frame):
#  
#   The frame the user interacts with.
#
#   Set other initial variables


    def __init__(self, parent=None):
#
#       Description:
#           Initializes frame, sets title and resolution.
#
#       Inputs:
#           parent - Frame containing Results window. Default is None.
###
        Frame.__init__(self, parent=None)
        self.master.title("Branching Fraction Calculator")
        self.master.geometry("1200x700")
        self.master.option_add('TkFixedFont','mono 12')
        self.pack(expand=YES, fill=BOTH)
        self.make_widgets()
        self.get_lifetimes()
        self.firstlog = 0
       
    def get_lifetimes(self):
#     
#       Description:
#           Imports level lifetimes from .lev files.
#      
        self.lifetimes = {}
        self.upper_value = {}

        for line in open(lifetime_file, 'r').readlines():
            line = line.split()
            upper_level = line[0]
            self.upper_value[upper_level] = line[2]
            try: lifetime = eval(line[4])
            except: lifetime = None
            self.lifetimes[upper_level] = lifetime
            self.upper_value[upper_level] = line[2]

    def log(self, text):
#      
#       Description:
#            Writes text to log file and displays line in ScrolledText on main display.
#
#       Inputs:
#           text - text to be written to log and put on main display.
#      
        self.logfile.write(text+'\n')
        self.display_text.insert('end', text+'\n')

    def set_output(self, event):
#       
#       Description:
#           Sets log_name to entry contents after pressing Enter key (event).
#     
        log_name = self.file_name_entry.get()+'.log'
        try:
        	close(self.logfile)                            # Check if file open already - if so, close it
        	                                               # and re-open with new filename.
        	self.logfile = open(log_name,'a',buffering=1)   # Line buffering forces immediate write
        except:
        	self.logfile = open(log_name,'a', buffering=1)    # Line buffering forces immediate write
       	self.log(f'Output file name set to {log_name}\n')
        event.widget["foreground"] = "Black"

    def get_upper_level(self, event):
#     
#       Description:
#           Sets self.upper_level to entry contents after pressing Enter key (event).
#           First displays table of previously identified transitions and lifetime
#           for given upper level, then displays the wavenumbers of transitions,
#           lower levels, SNRs, Intensities, LS configurations, and notes for observed lines.
#      
        self.upper_level = self.upper_level_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f'Upper level: {self.upper_level}')
        self.show_identified_lines()
        self.get_snrs_and_ints()
        self.show_values()

    def show_identified_lines(self):
#     
#       Description
#           Finds and displays all previously identified lines associated with
#           upper level in the .sort file.
#      

        self.known_lines = []

        # Check level file for given level key.
        for line in open(lifetime_file, 'r').readlines():
            line = line.split()

            if len(line) > 0:

                level_file_key = line[0]

                # If level Key is found -> return upper level conf key
                if(level_file_key == self.upper_level):
                    upper_energy_key = line[5]
                    self.known_lines = grep_open(upper_energy_key, self.id_lines_file)   # JJL: use grep_open for Windows compatibility
                    break

        self.log(f'Lifetime of upper level: {self.lifetimes[self.upper_level]} ns\n')
        self.log(f'Identified transitions from upper level {self.upper_level}:')
        self.log(''.join(['-']*86))
        self.log(f'| Wavenumber | L Energy   | U Energy   | L Conf      | U Conf      | Notes           |')
        self.log(''.join(['-']*86))

        for transition in self.known_lines:

            transition = transition.split()
            transition_key = transition[5]
            transition_key = transition_key.strip('*')

            if (transition_key == upper_energy_key):

                wavenumber = eval(transition[0])
                l_energy   = eval(transition[1])
                u_energy   = eval(transition[2])
                l_conf     = transition[3]
                u_conf     = transition[5]
                try: note  = transition[6]
                #try: print(transition[6])
                except IndexError: note = ''
                self.log(f'| {wavenumber:>10.3f} | {l_energy:>10.3f} | '
                         f'{u_energy:10.3f} | {l_conf:>11s} | '
                         f'{u_conf:>11s} | {note:15s} |')

        self.log(''.join(['-']*86)+'\n')

    def get_snrs_and_ints(self):
#  
#       Description:
#           
#           Creates wavenumbers list, and SNRs, Intensities, and transition IDs
#           dictionaries for upper level.
#      
        self.snrs = {} # {(spectrum, lower wavenumber) : snr}
        self.intensities = {}  # {(spectrum, lower wavenumber) : intensity}
        self.unc={}
        self.transition_ids = {}  # {transition wavenumber : lower wavenumber}
        root_npts={}
        self.w_maxI={}            # {(spectrum): Wavenumber of transition with max intensity in a spectrum}
        calunc_per_1000={}

        # Get all lines in files containing this level and cross reference against known linelist for wavenumbers.

        # Check level file for given lower level key.
        levels = open(lifetime_file, 'r').readlines()

        for spectrum in all_spectrum_files:

            with open(spectrum) as f:
                params = f.readline().split()
            resoln=float(params[0])  # Used for no. points/fwhm
            lowE=float(params[1])    # Used to get range for 7% unc.
            hiE=float(params[2])     
                                     # 7% from +- 5% discrepancy in D2 measurements
                                     # converts to unc of 2x5/sqrt(3) = 5.8%
                                     # added onto lamp unc. of 3.5% = 6.8% - round to 7 %
                                     # Calunc now set in branchrc.py
            calunc_per_1000[spectrum] = unc_cal/(hiE-lowE) # Calibration uncertainty per 1000 cm-1.
                                              # Estimated from 7%/(sig_hi - sig_low)

            transitions = grep_open("- "+self.upper_level, spectrum)

            self.w_maxI[spectrum] = 0
            maxI = 0
            
            for transition in transitions:

                transition = transition.split()
                lower_level = transition[6]
                lower_level = lower_level.strip('*')
                self.snrs[spectrum, lower_level] = eval(transition[0])
                self.intensities[spectrum, lower_level] = eval(transition[1])
                if float(transition[1])> maxI :        # Find the wavenumber of the strongest line
                    maxI = float(transition[1])        # Used for calibration uncertainty
                    self.w_maxI[spectrum] = float(transition[4])
                    
                root_npts[spectrum,lower_level]=(0.001*float(transition[2])/resoln)    # no. points/FWHM.
                self.unc[spectrum,lower_level]=2.25/(self.snrs[spectrum,lower_level]**2 *root_npts[spectrum,lower_level])
                                                                             # Variance
                for level in levels:
                    level = level.split()

                    if (level[0] == lower_level):
                        lower_level_key = level[5].strip('*')

                for line in self.known_lines:
                    line = line.split()
                    if (line[3].strip('*') == lower_level_key):
                        self.transition_ids[eval(line[0])] = lower_level  
        self.wavenumbers = sorted(self.transition_ids.keys(), reverse=True)
   
        for spectrum,lower_level in self.unc :
            wnum =  (list(self.transition_ids.keys())[list(self.transition_ids.values()).index(lower_level)])
            calunc = calunc_per_1000[spectrum]*(wnum-self.w_maxI[spectrum])/1000
            self.unc[spectrum,lower_level] = math.sqrt(calunc*calunc+self.unc[spectrum,lower_level])

    def show_values(self):
#   
#       Description:
#           Outputs SNRs and Intensities for upper level in each .I file
#   
        file_num = len(all_spectrum_files)
        width = 27+18*file_num

        self.log(f'SNRs and intensities of observed transitions from '
                 f'upper level {self.upper_level}:')
        self.log(''.join(['-']*width))

        file_name_line = f'| File Name:              | '
        for spectrum in all_spectrum_files: file_name_line += f'{spectrum}  | '

        self.log(file_name_line)
        self.log(''.join(['-']*width))
        self.log(f'| Wavenumber | L Level    | '+f'SNR | Intensity | '*file_num)
        self.log(''.join(['-']*width))

        for wavenumber in self.wavenumbers:

            data_line = f"| {wavenumber:>10.3f} | {self.transition_ids[wavenumber]:>10s} | "

            for spectrum in all_spectrum_files:

                try:
                     data_line += (f"{round(self.snrs[spectrum, self.transition_ids[wavenumber]]):>3d} | "
                                  f"{round(self.intensities[spectrum, self.transition_ids[wavenumber]]):9d} | ")

                except: data_line += f"--- | --------- | "

            self.log(data_line)

        self.log(''.join(['-']*width)+'\n')

    def get_reference_level(self, event):
#      
#       Description:
#           Sets self.reference_level to entry contents after pressing Enter key (event).
#           Normalizes the Intensities of each line with respect to the reference level,
#           with its intensity set to 1000, then displays the wavenumbers of transitions,
#           lower levels, SNRs, Intensities, LS configurations, and notes for observed lines.
#     
        self.reference_level = self.reference_level_entry.get()
        event.widget["foreground"] = "Black"
        self.normalize_spectrum()
        self.log(f"Normalizing each spectrum with respect to level {self.reference_level} = 1000.\n")
        self.show_values()

    def normalize_spectrum(self):
#     
#       Normalizes each emission spectrum file's intensities, with self.reference_level
#       having an intensity of 1000.
#      
        for spectrum in all_spectrum_files:

            try:

                normal = self.intensities[spectrum, self.reference_level]

            except:

                self.log(f"Spectrum {spectrum} does not contain level {self.reference_level}\n")
                continue

            for lower_level in self.transition_ids.values():

                try: self.intensities[spectrum, lower_level] = 1000*self.intensities[spectrum, lower_level]/normal
                except KeyError: pass

    def get_transfer_file(self, event):
#     
#       Description:
#           Sets self.transfer_file to entry contents after pressing Enter key (event).
#    
        self.transfer_file = self.transfer_file_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f"Transfer spectrum file selected: {self.transfer_file}\n")

    def get_transfer_level(self, event):

        self.normal_level = self.normal_level_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f"Normalized levels in spectrum {self.transfer_file} using weighted average of level {self.normal_level} in other spectra")
        self.normalize_from_transfer_level()
        self.show_values()

    def normalize_from_transfer_level(self):
# 
#       Description:
#           Finds the weighted average intensity for the transfer level for all files
#           containing it, and uses this to renormalize all other intensities in the transfer file.
#     
        weighted_int = 0
        sum_weight = 0
        sqsum_int = 0

        for spectrum in all_spectrum_files:

            if (((spectrum,self.normal_level) in self.intensities.keys()) and ((spectrum,self.reference_level)) in self.intensities.keys()):
                
                # Intensity weighted
                weight = 1/(self.unc[spectrum,self.normal_level])**2
                if weight > 555:           # 555 = 1/(0.06**2/2)
                                           # Cap on unc. from calibration, as line ratio.
                                           # So all lines with SNR>24 get same weight.
                                           # Omit lamp unc. here as it's the same for all spectra.
                    weight = 555
                weighted_int += self.intensities[spectrum, self.normal_level] * weight
                sum_weight += weight
                sqsum_int += self.intensities[spectrum, self.normal_level]**2

        # Take weighted average/weight as effective SNR for transfer level
        # Renormalize SNRs by 1/sqrt(1/SNR(meas)^2 + 1/SNR(calc)^2)
        # -- this is root sum of squares of uncertainties

        avg_int = weighted_int / sum_weight
        sqsum_int = sqsum_int**(0.5) / sum_weight

        unc_renorm = (1/(self.unc[self.transfer_file,self.normal_level]**2) + sum_weight)**(-0.5)
        int_renorm = avg_int / self.intensities[self.transfer_file, self.normal_level]
        
        for lower_level in self.transition_ids.values():

            try:

                self.intensities[self.transfer_file, lower_level] *= int_renorm

                # If level is transfer level, then SNR = unc_renorm
                if lower_level == self.normal_level:
                    self.unc[self.transfer_file, lower_level] = unc_renorm
                
                # Otherwise add uncertainties to unc. of ref. level
                else:
                    self.unc[self.transfer_file, lower_level] = (self.unc[spectrum, lower_level]**(2) + unc_renorm**(2))**(0.5)

            except Exception as e: pass
        

    def Delfil(self,event):        # Put back ability to delete lines (GN, July 21)

        self.delfile = self.Dfilenentry.get()
        event.widget["foreground"] = "Black"
        self.log(f"Deleting value from file {self.delfile}")
        
    def DelLine(self,event):       # Put back ability to delete lines (GN, July 21)
        self.delval =self.Delentry.get()
        del self.snrs[self.delfile,self.delval]
        del self.intensities[self.delfile,self.delval]
        self.log(f'{self.delfile}, {self.delval} deleted from list' )
        self.show_values()       

    def get_id_lines(self, event):

        self.id_lines_file = self.id_lines_entry.get()
        event.widget["foreground"] = "Black"
        self.log(f"ID'd lines file selected: {self.id_lines_file}\n")

    def results(self):
#  
#       Description:
#           Routine is called when the 'Results' button is pressed.
#           Calculates the branching fractions and transition probabilities
#           Intensities, LS configurations, and notes for observed lines.
#  
        total_int = 0
        level_keys = self.transition_ids.values()
        avg_int = {}
        sum_weight = {}
        sqsum = {}
        variance = {}
        stdev = {}
        waveno = {}
        upper_level = {}
        E1 = {}
        jupp = {}

        BFsq = 0

        for level_key in level_keys:

            avg_int[level_key] = 0
            sum_weight[level_key] = 0
            sqsum[level_key] = 0

            for spectrum in all_spectrum_files:

                if (spectrum, level_key) in self.intensities.keys():
                    weight = 1/(self.unc[spectrum,level_key])**2
                    if weight > 555:           # 555 = 1/(0.06**2/2)
                                               # Cap on unc. from calibration, as line ratio.
                                               # So all lines with SNR>24 get same weight.
                                               # Omit lamp unc. here as it's the same for all spectra.
                        weight = 555
                        
                    avg_int[level_key] += self.intensities[spectrum, level_key] * weight
                    sum_weight[level_key]  += weight
                    sqsum[level_key]   += self.intensities[spectrum, level_key]**(2)

            avg_int[level_key] /= sum_weight[level_key]

            # This line is fractional weighted standard deviation of intensities (u(I)/I in Sikstrom eqn 6)
            # Calibration uncertainty has already been taken into account in calculation of the individual uncs.
            stdev[level_key] = 1/math.sqrt(sum_weight[level_key])

            total_int += avg_int[level_key]

#
# Calculate and print out the residuals (lines for which calcluated A values exist but are not seen in experiment)
# 

        self.log('Residuals:')
        
        get_calculations(E1)
        sumAresid=0
        Aval_seen=[]
        theoval={}
        for (x,y) in E1.keys():  
            if x==self.upper_level:
                for wno in self.wavenumbers: 
                    if -discrim  < wno - float(y) < discrim:    # Save all the theo. values that are in experiment.
                        Aval_seen.append((x,y))  
                        theoval[x,wno]=E1[x,y]
        for (x,y) in E1.keys():                         # Sum up all the theo. values that are not seen for a residual.
            if x==self.upper_level and (x,y) not in Aval_seen:
                sumAresid += E1[x,y]
                self.log(f'{y:>10s} | {E1[x,y]:>7.4f}')   # Print these lines to verify they're really not seen.
        frac_resid = sumAresid*self.lifetimes[self.upper_level]/1000.    # Convert to a fractional residual
        total_int *= (1+frac_resid)                       # Use residual as correction to the total intensity.
       
#
# Print the header for the results
#
        self.log(''.join(['-']*108)+'\n')
        self.log('Summary:')
        self.log(f'Level = {self.upper_level:>10s}, Lifetime = {self.lifetimes[self.upper_level]} ns, residual =, {frac_resid*100.:>6.3f} %')
        self.log(''.join(['-']*108))
        self.log(f'| Wavenumber | L Level    | Weight | B. Fraction | Uncertainty | Trans. Prob. | Uncertainty |  Theoretical |')
        self.log(''.join(['-']*108))

        for level_key in level_keys:            # This is 2nd term in eqn 7 of Sikstrom. 
            BFsq += (avg_int[level_key]/total_int)**2 * stdev[level_key]**2

# Correction for the residual (50 % uncertainty estimate on residual)
        BFsq += frac_resid**2 * 0.25

        for wavenumber in self.wavenumbers:

            level_key = self.transition_ids[wavenumber]
                                                       # First two terms of eqn 7 of Sikstrom (BF uncertainty)
            variance[level_key] = stdev[level_key]**2 * (1 - 2*avg_int[level_key]/total_int)+ BFsq

            avg_int[level_key] /= total_int
            aval = (1000 * avg_int[level_key]) / (self.lifetimes[self.upper_level])
            unc_aval = 100 * math.sqrt(variance[level_key] + lifetime_unc**2)    # Combine BF unc and lifetime unc in quadrature

            try: 
                theo_val = theoval[self.upper_level,wavenumber]
            except:
                theo_val = 0

            self.log(f'| {wavenumber:>10.3f} | {level_key:>10s} | {round(sum_weight[level_key]):>6d} | '
                     f'{100*avg_int[level_key]:>10.3f}% | {100*math.sqrt(variance[level_key]):>10.1f}% | '
                     f'{aval:>12.3f} | {unc_aval:>10.0f}% | {theo_val:>12.3f} |')


        self.log(''.join(['-']*108)+'\n')

    def quit(self):
    	self.logfile.close()
    	sys.exit()

    def make_widgets(self):
#     
#       Description:
#           Creates the widgets the user interacts with.
#     

        self.display_frame = Frame(self)
        self.display_frame.pack(side=RIGHT, fill=BOTH, expand=YES)
        self.display_scroll_y = Scrollbar(self.display_frame, orient=VERTICAL)
        self.display_scroll_y.pack(side=RIGHT, fill=Y)
        self.display_text = Text(self.display_frame, wrap=NONE, width=0)
        self.display_text.pack(side=TOP, fill=BOTH, expand=YES)
        self.display_scroll_y.config(command=self.display_text.yview)
        self.display_scroll_x = Scrollbar(self.display_frame, orient=HORIZONTAL, command=self.display_text.xview)
        self.display_scroll_x.pack(side=TOP, fill=X)

        self.display_text.config(xscrollcommand=self.display_scroll_x.set,
                                 yscrollcommand=self.display_scroll_y.set)

        self.widgets_frame = Frame(self)
        self.widgets_frame.pack(side=LEFT, fill=BOTH)
        self.toolbar_frame = Frame(self.widgets_frame, height=30)
        self.toolbar_frame.pack(side=TOP, fill=X)
        Button(self.toolbar_frame, text='Quit', command=self.quit).pack(side=LEFT, fill=X, expand=YES)
        Button(self.toolbar_frame, text='Results', command=self.results).pack(side=LEFT, fill=X, expand=YES)

        self.label_frame = Frame(self.widgets_frame)
        self.label_frame.pack(side=LEFT, fill=BOTH)
        self.entry_frame = Frame(self.widgets_frame)
        self.entry_frame.pack(side=RIGHT, fill=BOTH)

        Label(self.label_frame, text="Save As:").pack(side=TOP, fill=X, pady=4.5)
        self.file_name_entry = Entry(self.entry_frame, foreground="Blue")
        self.file_name_entry.pack(side=TOP, fill=X, pady=2)
        self.file_name_entry.insert(0, file_name)
        self.file_name_entry.bind("<Key-Return>", self.set_output)

        Label(self.label_frame, text="ID'd lines file:").pack(side=TOP, fill=X, pady=4.5)
        self.id_lines_entry = Entry(self.entry_frame, foreground="Blue")
        self.id_lines_entry.pack(side=TOP, fill=X, pady=2)
        self.id_lines_entry.insert(0, identified_lines_file)
        self.id_lines_entry.bind("<Key-Return>", self.get_id_lines)

        Label(self.label_frame, text="Upper Level:").pack(side=TOP, fill=X, pady=4.5)
        self.upper_level_entry = Entry(self.entry_frame, foreground="Blue")
        self.upper_level_entry.pack(side=TOP, fill=X, pady=2)
        self.upper_level_entry.insert(0, upper_level)
        self.upper_level_entry.bind("<Key-Return>", self.get_upper_level)

        Label(self.label_frame, text="Reference Level:").pack(side=TOP, fill=X, pady=4.5)
        self.reference_level_entry = Entry(self.entry_frame, foreground="Blue")
        self.reference_level_entry.pack(side=TOP, fill=X, pady=2)
        self.reference_level_entry.insert(0, reference_level)
        self.reference_level_entry.bind("<Key-Return>", self.get_reference_level)

        Label(self.label_frame, text="Transfer File:").pack(side=TOP, fill=X, pady=4.5)
        self.transfer_file_entry = Entry(self.entry_frame, foreground="Blue")
        self.transfer_file_entry.pack(side=TOP, fill=X, pady=2)
        self.transfer_file_entry.insert(0, transfer_file)
        self.transfer_file_entry.bind("<Key-Return>", self.get_transfer_file)

        Label(self.label_frame, text="Transfer Level:").pack(side=TOP, fill=X, pady=4.5)
        self.normal_level_entry = Entry(self.entry_frame, foreground="Blue")
        self.normal_level_entry.pack(side=TOP, fill=X, pady=2)
        self.normal_level_entry.insert(0, normal_level)
        self.normal_level_entry.bind("<Key-Return>", self.get_transfer_level)
     
        Label(self.label_frame, text="Delete file:").pack(side=TOP, fill=X, pady=4.5)
        self.Dfilenentry = Entry(self.entry_frame, foreground="Blue")
        self.Dfilenentry.pack(side=TOP, fill=X, pady=2)
        self.Dfilenentry.insert(0, delfile)
        self.Dfilenentry.bind("<Key-Return>", self.Delfil)
        
        Label(self.label_frame, text="Delete level:").pack(side=TOP, fill=X, pady=4.5)
        self.Delentry = Entry(self.entry_frame, foreground="Blue")
        self.Delentry.pack(side=TOP, fill=X, pady=2)
        self.Delentry.insert(0, delval)
        self.Delentry.bind("<Key-Return>", self.DelLine)
     

if __name__ == '__main__':

    BranchingFractionCalc().mainloop()
