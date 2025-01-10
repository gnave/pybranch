"""
A script for calculating branching fractions from sets of intensity ratios
measured in various emission spectra.
"""
import logging
import math
from branchgrep import grep_open
from branchrc import set_parameters
from glob import glob
from itertools import product
from os import popen, path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from scipy.special import voigt_profile
from tkinter import *
from tkinter import messagebox
from tkinter.scrolledtext import ScrolledText
from tkinter import simpledialog
from tkinter import filedialog
from nso_to_hdf import *

class BranchingFractionCalc(Frame):
    """
    The frame the user interacts with.
    """
    # Set other initial variables
    spectrum_extension=set_parameters()[8]
    all_spectrum_files = glob(spectrum_extension)
    id_lines_file = set_parameters()[5]
    calc_file = set_parameters()[9]
    plotwin_length = set_parameters()[11]
    life_unc = set_parameters()[7]
    discrim = set_parameters()[10]

    OldMenu = 0

    def __init__(self, parent=None):
        """
        Description:
            Initializes frame, sets title and resolution.

        Inputs:
            parent - Frame containing Results window. Default is None.
        """
        Frame.__init__(self, parent=None)
        self.master.title("Branching Fraction Calculator")
        self.master.geometry("1200x700")
        self.master.option_add('*Font','mono 8')
        self.pack(expand=YES, fill=BOTH)
        self.get_lifetimes()
        self.make_widgets()

    def get_calculations(self,E1):
	  #
	  # This function gets the calculated A values for calculating a residual
	  #
	  
        for line in open(self.calc_file, 'r').readlines():
            line = line.split()
            waveno = line[0]
            upper_level = line[-1]
            E1[upper_level,waveno] = float( line[1])
        return(E1)
      
    def get_lifetimes(self):
        """
        Description:
            Imports level lifetimes from .lev files.
        """
        self.lifetimes = {}
        self.upper_value = {}
        self.levels = []
        lifetime_files = glob('*.lev')

        for lifetime_file in lifetime_files:

            for line in open(lifetime_file, 'r').readlines():
                line = line.split()
                upper_level = line[6]
                self.upper_value[upper_level] = line[2]
                try: 
                    lifetime = eval(line[5])
                    self.levels.append(line[6]) 
                except: lifetime = None
                self.lifetimes[upper_level] = lifetime
                self.upper_value[upper_level] = line[2]
                
    def log(self, text):
        """
        Description:
            Writes text to log file and displays line in ScrolledText on main display.

        Inputs:
            text - text to be written to log and put on main display.
        """
        logging.info(text)
        self.display_text.insert('end', text+'\n')

    def set_output(self, event):
        """
        Description:
            Sets log_name to entry contents after pressing Enter key (event).
        """
        log_name = self.file_name_entry.get()
        logging.basicConfig(filemode='a', filename=log_name+'.log', level=logging.INFO, format='%(message)s')
        self.log(f'Log file name set to {log_name}.log\n')
        event.widget["foreground"] = "Black"

    def get_upper_level(self):
        """
        Description:
            Sets self.upper_level to entry contents after pressing Enter key (event).
            First displays table of previously identified transitions and lifetime
            for given upper level, then displays the wavenumbers of transitions,
            lower levels, SNRs, Intensities, LS configurations, and notes for observed lines.
        """
        
        self.upper_level = self.Upper.get()
        self.log(f'Upper level: {self.upper_level}')
        self.show_identified_lines()
        self.get_snrs_and_ints()
        self.show_values()
        Upper_lev = StringVar()

    def show_identified_lines(self):
        """
        Description
            Finds and displays all previously identified lines associated with
            upper level in the .sort file.
        """

        self.known_lines = []
        energy_levels_file = set_parameters()[6]

        # Check level file for given level key.
        for line in open(energy_levels_file, 'r').readlines():
            line = line.split()

            if len(line) > 0:

                level_file_key = line[6]

                # If level Key is found -> return upper level conf key
                if(level_file_key == self.upper_level):
                    upper_energy_key = line[6]
                    # JJL
                    self.known_lines = grep_open(upper_energy_key, self.id_lines_file)
                    # self.known_lines = popen(f"grep {upper_energy_key} {self.id_lines_file}", 'r').readlines()
                    break
    

        self.log(f'Lifetime of upper level: {self.lifetimes[self.upper_level]} ns\n')
        self.log(f'Identified transitions from upper level {self.upper_level}:')
        self.log(''.join(['-']*86))
        self.log(f'| Wavenumber | L Conf       | Intensity and Notes                                    |')
        self.log(''.join(['-']*86))

        for tr_line in self.known_lines:

            transition = tr_line.split()
            transition_key = transition[5]

            if (transition_key == upper_energy_key):
                wavenumber = eval(transition[0])
                l_conf     = transition[3]
                u_conf     = transition[5]
                # Any intensities should be after the upper config. Notes will be after this
                try:
                    inten_pos = tr_line.index(u_conf)+len(u_conf)
                    note = (tr_line[inten_pos:].rstrip())
                except:
                    note = " "

                self.log(f'| {wavenumber:>10.3f} | {l_conf:>12s} | {note:<54s} |') 


        self.log(''.join(['-']*86)+'\n')

    def calc_gaussian_width_from_voigt(self,width,damping):    
        """
        This function uses equation 8 from Kielkopf, JOSA 63, 987 (1973)
        It's used to get the Gaussian width of a Voigt profile given the Voigt width
        and damping parameter. There's a factor of sqrt(log(2)) in this equation that I don't understand
        """
        eta = 0.099
        A = 1+eta*np.log(2)
        B = eta*np.log(2)

        gauss = width * np.sqrt(1-A*damping+B*damping*damping)

        return(gauss)

    def voigt_fit(self,specfile, wnum, x):

        """
        Find out if there are any lines at wnum in the linelist corresponding to specfile.
        If so, calculate a voigt profile using the grid in x
        """

        window = 0.1      # Lines within this window are matched. Units are cm-1
        # Find a line at wavenumber wnum in the linelist corresponding to specfile and calculate its profile
        linefit = NONE
        linelist = read_linelist(specfile)
        for j in linelist:
            if wnum - window < j['sig'] < wnum + window:
                linefit=j
                break
            #       linefit = (x for x in linelist if  wnum - 1. < x['sig']< wnum + 1. )

        # Do nothing if there is no line in the list       
        # and return NONE         
        if linefit == NONE:
            return

        # Xgremlin widths are in 0.001 cm-1, so convert to cm-1.
        # voigt_profile also uses HWHM, rather than FWHM, so divide by 2.
        width = linefit['width']/2000.           

        # Xgremlin damping parameter runs from 1 to 26, so convert to go from 0 to 1
        damping = (linefit['dmping']-1)/25

        # Calculate Gaussian width.
        gauss = self.calc_gaussian_width_from_voigt(width,damping)  
        # Convert to std. devn of normal distribution using 1/2sqrt(2ln2) - needed for voigt_profile
        gauss = gauss/(np.sqrt(2*np.log(2)))
        lorentz = width*damping
   
        y = voigt_profile(x-linefit['sig'],gauss,lorentz)
        y = y/voigt_profile(0,gauss,lorentz) * linefit['xint']

        return(x,y, linefit)

    def PlotLevel(self):
        # First, find the file and the level. Set the wavenumber
        file_path = filedialog.askopenfilename(title="Select file", filetypes=[("Spectrum file",('*.dat'))])
        specfile = file_path[:-4]
        lev = self.showlev.get()  
        wnum =  (list(self.transition_ids.keys())[list(self.transition_ids.values()).index(lev)])    

        # Read the header in and set the parameters 
        if exists(specfile + ".hdr"):
            header = read_header(specfile)
        else:
            messagebox.showinfo("Error","Header file not found")
            return()

        if 'Complex' in header['data_is'] :
            cmplx = int(2)           
        else:
            cmplx = int(1)
        wstart = float(header['wstart'])
        delw = float(header['delw'])
        wavcorr = float(header['wavcorr'])

        if wavcorr != 0:
            wstart = wstart*(1+wavcorr)
            delw = delw*(1+wavcorr)

        # Set the starting index to plot plotwin_length/2 points each side of the line
        # Make sure that you start on a real point if its a complex spectrum.
        startidx = cmplx * int((wnum - wstart)/delw - self.plotwin_length/2.)
        # Calculate the starting wavenumber from the actual starting index, rather than what was requested
        wref = wstart+startidx*delw/cmplx       
        npts = self.plotwin_length*cmplx

        # Open the file and read in the spectrum
        with open(specfile+".dat","rb") as fb:
            tmp = np.fromfile(fb,np.float32)
            spec = tmp[startidx:startidx+npts:cmplx]*float(header['rdsclfct'])
        
        # Plot the line

        x=np.linspace(wref,wref+(self.plotwin_length-1)*delw,self.plotwin_length)
        plt_title = path.basename(file_path).split('/')[-1]
        fig,ax = plt.subplots(num=plt_title[:-4])

        fig.suptitle(lev)
        ax.set_xlabel('Wavenumber /cm$^{-1}$')
        ax.set_ylabel('Signal-to-noise Ratio')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%9.3f'))
        ax.plot(x,spec)

        # Leave room at the bottom of the plot to print the line parameters
        fig.subplots_adjust(bottom=0.2)

        # Plot the fit, if there is one

        # If there is a wavenumber correction factor for this file, subtract it off in order
        # to match the wavenumbers in the .lin file
        if wavcorr != 0:
            fig.text(.5,.06,("wavcorr = %4.2e applied" %wavcorr), ha='center' )
            wnum = wnum*(1-wavcorr)
        try:
            (fitx,fity, linefit) = self.voigt_fit(specfile, wnum, x)
            # re-apply the wavenumber correction factor to the data
            fitx=fitx*(1+wavcorr)
            ax.plot(fitx,fity,ls='dotted')

            params = "sig = %9.3f cm$^{-1}$; Int. = %6.0f; FWHM = %6.3f cm$^{-1}$; Damp. = %6.3f "   \
                % (linefit['sig']*(1+wavcorr),linefit['xint'], linefit['width']/1000., (linefit['dmping']-1)/25)
            fig.text(.5, .02, params, ha='center')
        except:
            messagebox.showinfo("Info","No entry in linelist at this wavenumber")           

        fig.show()

    def MakeMenus(self):
        ulev_key =[]

        for wno in self.wavenumbers:      
            ulev_key.append(self.transition_ids[wno])
              
        self.reflev = StringVar()
        self.normlev = StringVar()
        self.dellev = StringVar()

        # Destroy any existing widgets

        if self.OldMenu:

            self.RefLevelButton.destroy()
            self.RefLevelMenu.destroy()
            self.RefFileLabel.destroy()
            self.Reffilemenu.destroy()
            self.NormalLevelButton.destroy()
            self.NormalLevelMenu.destroy()
            self.DelFileLabel.destroy()
            self.Dfilemenu.destroy()
            self.DelLevelButton.destroy()
            self.DelLevelMenu.destroy()
            self.ShowLevelLabel.destroy()
            self.ShowLevelMenu.destroy()
                
        # Create the rest of the buttons and dropdown menu widgets

        self.RefLevelButton = Button(self.label_frame, text="Ref. Level:", width=13, command=self.get_reference_level)
        self.RefLevelButton.grid(row=2,column=0,sticky=SW)
        self.RefLevelMenu = OptionMenu(self.label_frame, self.reflev, *ulev_key )
        self.RefLevelMenu.config(width=13)
        self.RefLevelMenu.grid(row=2,column=1)

        self.RefFileLabel = Label(self.label_frame, text="Rescale File:")
        self.RefFileLabel.grid(row=3,column=0, sticky=SW,pady=1)
        ref_file_clicked = StringVar()
        self.RefFileMenu = OptionMenu(self.label_frame, ref_file_clicked, command=self.get_reference_file, *self.all_spectrum_files)
        self.RefFileMenu.config(width=13)
        self.RefFileMenu.grid(row=3,column=1)

        self.NormalLevelButton = Button(self.label_frame, text="Rescale Level:", width=13,command=self.get_normal_level)
        self.NormalLevelButton.grid(row=4,column=0,sticky=SW)
        self.NormalLevelMenu = OptionMenu(self.label_frame,self.normlev, *ulev_key )
        self.NormalLevelMenu.config(width=13)
        self.NormalLevelMenu.grid(row=4,column=1)

        self.DelFileLabel = Label(self.label_frame, text="Delete file:")
        self.DelFileLabel.grid(row=5,column=0, sticky=SW,pady=1)
        dfile_clicked = StringVar()
        self.DfileMenu = OptionMenu(self.label_frame, dfile_clicked, command=self.Delfil, *self.all_spectrum_files)
        self.DfileMenu.config(width=13)
        self.DfileMenu.grid(row=5,column=1)     

        self.DelLevelButton = Button(self.label_frame, text="Delete level", width=13,command = self.DelLine)
        self.DelLevelButton.grid(row=6,column=0,sticky=SW,pady=1)
        self.DelLevelMenu = OptionMenu(self.label_frame, self.dellev, *ulev_key )
        self.DelLevelMenu.config(width=13)
        self.DelLevelMenu.grid(row=6,column=1)

        self.ShowLevelLabel = Label(self.file_frame,text="Plot Line to")
        self.ShowLevelLabel.grid(row=0,column=0)
        self.showlev = StringVar()
        self.ShowLevelMenu = OptionMenu(self.file_frame, self.showlev, *ulev_key)
        self.ShowLevelMenu.config(width=13)
        self.ShowLevelMenu.grid(row=0,column=1)

        self.ShowFileLabel = Label(self.file_frame,text="in file")
        self.ShowFileLabel.grid(row=1,column=0)
        self.ShowFileButton = Button(self.file_frame,text="Select File",command=self.PlotLevel)
        self.ShowFileButton.grid(row=1,column=1)


        self.OldMenu = 1

    def get_snrs_and_ints(self):
        """
        Description:
            Creates wavenumbers list, and SNRs, Intensities, and transition IDs
            dictionaries for upper level.
        """
        self.snrs = {} # {(spectrum, lower wavenumber) : snr}
        self.intensities = {}  # {(spectrum, lower wavenumber) : intensity}
        self.unc={}
        self.transition_ids = {}  # {transition wavenumber : lower wavenumber}
        root_npts={}
        self.w_maxI={}            # {(spectrum): Wavenumber of transition with max intensity in a spectrum}
        calunc_per_1000={}

        # Get all lines in files containing this level and cross reference against known linelist for wavenumbers.

        energy_levels_file = set_parameters()[6]

        # Check level file for given lower level key.
        levels = open(energy_levels_file, 'r').readlines()

        for spectrum in self.all_spectrum_files:

            with open(spectrum) as f:
                params = f.readline().split()
            resoln=float(params[0])  # Used for no. points/fwhm
            lowE=float(params[1])    # Used to get range for 7% unc.
            hiE=float(params[2])     
                                     # 7% from +- 5% discrepancy in D2 measurements
                                     # converts to unc of 2x5/sqrt(3) = 5.8%
                                     # added onto lamp unc. of 3.5% = 6.8% - round to 7 %
            calunc_per_1000[spectrum] = 70/(hiE-lowE) # Calibration uncertainty per 1000 cm-1.
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

                    if (level[6] == lower_level):
                        lower_level_key = level[6].strip('*')
                for line in self.known_lines:
                    line = line.split()
                    if (line[3].strip('*') == lower_level_key):
                        self.transition_ids[eval(line[0])] = lower_level  
                     
        self.wavenumbers = sorted(self.transition_ids.keys(), reverse=True)
       
        self.MakeMenus()

        for spectrum,lower_level in self.unc :
            wnum =  (list(self.transition_ids.keys())[list(self.transition_ids.values()).index(lower_level)])
            calunc = calunc_per_1000[spectrum]*(wnum-self.w_maxI[spectrum])/1000
            self.unc[spectrum,lower_level] = math.sqrt(calunc*calunc+self.unc[spectrum,lower_level])

    def show_values(self):
        """
        Description:
            Outputs SNRs and Intensities for upper level in each .I file
        """
        file_num = len(self.all_spectrum_files)

        width = 46+18*file_num

        self.log(f'SNRs and intensities of observed transitions from '
                 f'upper level {self.upper_level}:')
        self.log(''.join(['-']*width))

        file_name_line = f'| File Name:              | '
        for spectrum in self.all_spectrum_files: 
            # Strip off the extension for printing
            spectrum = spectrum.replace(self.spectrum_extension[1:],"")
            file_name_line += f'{spectrum:>14s}  | '
        file_name_line += f'Mean    | Stdev  |'

        self.log(file_name_line)
        self.log(''.join(['-']*width))
        self.log(f'| Wavenumber | L Level    | '+f'SNR | Intensity | '*file_num + f'        |        |')
        self.log(''.join(['-']*width))

        for wavenumber in self.wavenumbers:

            data_line = f"| {wavenumber:>10.3f} | {self.transition_ids[wavenumber]:>10s} | "
            # Calculate the weighted mean intensity and standard deviation
            try:   
                (Imean,Istdev) = self.StdDev(self.transition_ids[wavenumber])
                Istdev = Imean*Istdev
            except:
                # If it fails, it means that the line has been deleted in all the spectra
                Istdev = 0
                Imean = 0
                # Remove this value from self.wavenumbers if its not seen in any spectra
                self.wavenumbers.remove(wavenumber)

            for spectrum in self.all_spectrum_files:
                
                try:                    
                    tmp = Imean/self.snrs[spectrum, self.transition_ids[wavenumber]]
                    devn = math.sqrt(tmp*tmp + Istdev*Istdev)
                    if -3.*devn < Imean -self.intensities[spectrum, self.transition_ids[wavenumber]] < 3.*devn:
                        flag = " "
                    else:
                        flag = "*"
                except:
                    flag = " "

                try:
                    data_line += (f"{round(self.snrs[spectrum, self.transition_ids[wavenumber]]):>3d} | ")
                    data_line += flag
                    data_line += (f"{round(self.intensities[spectrum, self.transition_ids[wavenumber]]):8d} | ")

                except: data_line += f"--- | --------- | "

            data_line += (f"{round(Imean):7d} | ")
            data_line += (f"{round(Istdev):>6d} |")

            self.log(data_line)

        self.log(''.join(['-']*width)+'\n')

    def get_reference_level(self):
        """
        Description:
            Sets self.reference_level to entry contents after pressing Enter key (event).
            Normalizes the Intensities of each line with respect to the reference level,
            with its intensity set to 1000, then displays the wavenumbers of transitions,
            lower levels, SNRs, Intensities, LS configurations, and notes for observed lines.
        """
        self.reference_level = self.reflev.get()
        self.normalize_spectrum()
        self.log(f"Normalizing each spectrum with respect to level {self.reference_level} = 1000.\n")
        self.show_values()

    def normalize_spectrum(self):
        """
        Normalizes each emission spectrum file's intensities, with self.reference_level
        having an intensity of 1000.
        """
        for spectrum in self.all_spectrum_files:

            try:

                normal = self.intensities[spectrum, self.reference_level]

            except:

                self.log(f"Spectrum {spectrum} does not contain level {self.reference_level}\n")
                continue

            for lower_level in self.transition_ids.values():

                try: self.intensities[spectrum, lower_level] = 1000*self.intensities[spectrum, lower_level]/normal
                except KeyError: pass

    def get_reference_file(self, ref_file):
        """
        Description:
            Sets self.reference_file to entry contents after pressing Enter key (event).
        """
        self.reference_file = ref_file
        self.log(f"Reference spectrum file selected: {self.reference_file}\n")

    def get_normal_level(self):

        self.normal_level = self.normlev.get()
        self.log(f"Normalized levels using level {self.normal_level} from spectrum {self.reference_file}")
        self.add_comment("Why did you rescale this spectrum?");
        self.normalize_all_spectra()
        self.show_values()

    def normalize_all_spectra(self):
        """
        Description:
            Finds the weighted average intensity for the normal level for all files
            containing it, and uses this to rescale all other intensities.
        """
        weighted_int = 0
        sum_weight = 0
        sqsum_int = 0

        for spectrum in self.all_spectrum_files:

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
        # Rescale SNRs by 1/sqrt(1/SNR(meas)^2 + 1/SNR(calc)^2)
        # -- this is root sum of squares of uncertainties

        avg_int = weighted_int / sum_weight
        sqsum_int = sqsum_int**(0.5) / sum_weight

        unc_renorm = (1/(self.unc[self.reference_file,self.normal_level]**2) + sum_weight)**(-0.5)
        int_renorm = avg_int / self.intensities[self.reference_file, self.normal_level]
        
        for lower_level in self.transition_ids.values():

            try:

                self.intensities[self.reference_file, lower_level] *= int_renorm

                # If level is transfer level, then SNR = unc_renorm
                if lower_level == self.normal_level:

                    self.unc[self.reference_file, lower_level] = unc_renorm
                

                # Otherwise add uncertainties to unc. of ref. level
                else:

                    self.unc[self.reference_file, lower_level] = (self.unc[spectrum, lower_level]**(2) + unc_renorm**(2))**(0.5)

            except Exception as e: pass
        
    def Delfil(self,dfile):        # Put back ability to delete lines (GN, July 21)
        self.delfile=dfile
        self.log(f"Deleting value from file {self.delfile}")
        
    def DelLine(self):       # Put back ability to delete lines (GN, July 21)
        self.delval =self.dellev.get()
        del self.snrs[self.delfile,self.delval]
        del self.intensities[self.delfile,self.delval]
        self.log(f'{self.delfile}, {self.delval} deleted from list' )   
        self.add_comment("Why did you delete it? ");     
        self.show_values()

    def StdDev(self,key):
        """
        Calculate the mean and uncertainty of intensity using self.intensities and self.unc
        """
        avg_int = 0
        sum_weight = 0
        sqsum = 0
 
        for spectrum in self.all_spectrum_files:

            if (spectrum, key) in self.intensities.keys():
                weight = 1/(self.unc[spectrum,key])**2
                if weight > 555:           # 555 = 1/(0.06**2/2)
                                               # Cap on unc. from calibration, as line ratio.
                                               # So all lines with SNR>24 get same weight.
                                               # Omit lamp unc. here as it's the same for all spectra.
                    weight = 555
                        
                avg_int += self.intensities[spectrum, key] * weight
                sum_weight  += weight
                sqsum   += self.intensities[spectrum, key]**(2)

        if sum_weight >0:
            avg_int /= sum_weight

            # This line is fractional weighted standard deviation of intensities (u(I)/I in Sikstrom eqn 6)
            # Calibration uncertainty has already been taken into account in calculation of the individual uncs.
            Istdev = 1/math.sqrt(sum_weight)
            return(avg_int,Istdev)
        else:
            return(NONE)

    def display(self):
        """
        Description:
            Displays the wavenumbers of transitions, lower levels, SNRs,
            Intensities, LS configurations, and notes for observed lines.
        """
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
        delval=[]

        BFsq = 0

        self.log('Residuals:')

        for level_key in level_keys:
            avg_int[level_key] = 0
            sum_weight[level_key] = 0
            sqsum[level_key] = 0
            try:
                (avg_int[level_key],stdev[level_key]) =  self.StdDev(level_key)
            except:
#                delval.append(level_key,)
                avg_int[level_key] = 0
                stdev[level_key] = 0
            total_int += avg_int[level_key]
             
#
# Calculate the residual
#  
        self.get_calculations(E1)
        sumA=0

        theoval={}
        for (x,y) in E1.keys():
            if x==self.upper_level:
                for wno in self.wavenumbers:
                    if -self.discrim  < wno - float(y) < self.discrim:
                        
                        delval.append((x,y))  
                        theoval[x,wno]=E1[x,y]
        for (x,y) in E1.keys():            
            if x==self.upper_level and (x,y) not in delval:
                sumA += E1[x,y]
                self.log(f'{y:>10s} | {E1[x,y]:>7.4f}')
        frac_resid = sumA*self.lifetimes[self.upper_level]/1000.
        total_int *= (1+frac_resid)
       
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
            unc_aval = 100 * ((variance[level_key]) + self.life_unc**(2))**(0.5)
            try: 
                theo_val = theoval[self.upper_level,wavenumber]
            except:
                theo_val = 0

            self.log(f'| {wavenumber:>10.3f} | {level_key:>10s} | {round(sum_weight[level_key]):>6d} | '
                     f'{100*avg_int[level_key]:>10.3f}% | {100*math.sqrt(variance[level_key]):>10.1f}% | '
                     f'{aval:>12.3f} | {unc_aval:>10.0f}% | {theo_val:>12.3f} |')


        self.log(''.join(['-']*108)+'\n')  
        self.add_comment("Comment on results")

    def add_comment(self,title_text):
        """
        Add a comment to the output file
        """
        comment = simpledialog.askstring(title="Comment",prompt=title_text)
        try:
            self.log(comment) 
        except:
            comment = " "
            self.log(comment)

    def comment_box(self):
        """
        For use with comment button to include a prompt
        """
        self.add_comment("Comment")
            
    def make_widgets(self):
        """
        Description:
            Creates the widgets the user interacts with.
        """

# First create the frame and windows to put them all in

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
        Button(self.toolbar_frame, width=13,text='Quit', command=self.quit).grid(row=0,column=0)
        Button(self.toolbar_frame, width=13,text='Results', command=self.display).grid(row=0,column=1)

        self.label_frame = Frame(self.widgets_frame)
        self.label_frame.pack(side=TOP, fill=BOTH)

        self.file_frame = Frame(self.widgets_frame )
        self.file_frame.pack(side=BOTTOM,fill=Y )

# Create initial set of widgets to set the log file and get the upper level 

        Label(self.label_frame, text="Log file name:").grid(row=0,column=0, sticky=SW)
        self.file_name_entry = Entry(self.label_frame, foreground="Blue")
        self.file_name_entry.grid(row=0,column=1)
        self.file_name_entry.insert(0, set_parameters()[0])
        self.file_name_entry.bind("<Key-Return>", self.set_output)
    
        self.Upper = StringVar()
        Button(self.label_frame,text='Upper level',width=13, command=self.get_upper_level).grid(row=1,column=0, sticky=SW,pady=1)
        self.Dfilemenu = OptionMenu( self.label_frame, self.Upper, *self.levels )
        self.Dfilemenu.config(width=13)
        self.Dfilemenu.grid(row=1,column=1)


        Button(self.file_frame, text="Add comment",command=self.comment_box).grid(row=5,column=0)
        
     
if __name__ == '__main__':

    BranchingFractionCalc().mainloop()
