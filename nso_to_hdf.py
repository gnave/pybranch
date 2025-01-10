import numpy as np
import sys
from os.path import exists
from struct import unpack

def read_linelist(specfile):
    """
    This reads a .lin file from Xgremlin and returns the whole
    thing as an array of dictionaries, with the dictionary keys
    being the parameters of each line
    """    
    sp={}
    linel=[]

    flin = open(specfile+".lin","rb")
    #
    # Read the header information in the .lin file
    #
    nlin=unpack("i",flin.read(4))[0]          # No. of lines in the file
    # print(nlin, "lines in file")
    linlen = unpack("i",flin.read(4))[0]      # Total number of valid data in file

    # Next line is to make up total of 320 bytes that is the prefix in the file
    tmp = flin.read(312)

    #
    # Now read in all the lines and add to the list of dictionaries
    #
    for k in range(nlin) :
        sp['sig'],sp['xint'],sp['width'],sp['dmping'],sp['itn'],sp['ihold'] = unpack("dfffhh",flin.read(24))
        sp['tags'] = flin.read(4)
        sp['epstot'],sp['epsevn'],sp['epsodd'],sp['epsran'],sp['spare']=unpack("fffff",flin.read(20))
        sp['ident']=flin.read(32)
        linel.append(sp.copy())

    return(linel)

def read_header(specfile):
    #
    # Create the metadata from the header file.
    #
    header = {}
    hdr = open(specfile + ".hdr")

    for line in hdr:
        if line[0] == "/" or line[0:3] == "END":
            continue
        if line[0:8] == "continue" :
            val = val + line[9:32]
        else:
            key = line[0:8].rstrip()
            val = line[9:32]
        if line[0:2] == "id":
            val = line[9:80]
        header[key]=val
        header[key+"_comment"] = line[34:80]

    return(header)

def read_data(specfile):
    with open(specfile + ".dat","rb") as f:
         spec = np.fromfile(f,np.float32)
         if len(spec) - int(header["npo"]) != 0 :
             print(f'No. of points does not match npo: npo = {header["npo"]}, length = {len(spec)}')
         else:
             print (f'{len(spec)} points read')

    return(spec)

def read_response(specfile):
    return()

