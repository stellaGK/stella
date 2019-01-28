from stella_dirs import *
from stella_vmec import *
from numpy import *
from pylab import *
from struct import *
from scipy import *
from matplotlib import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.io import netcdf
plt.rcParams.update({'font.size': 28})
plt.rcParams['lines.linewidth'] = 2

import tabCompleter

from plotbox import *
from tabCompleter import *
import struct
import physcon as pc
import time
from stella_read import *


def inputlist(case):
    i      = 1
    inlist = []
    
    for f in listdir(outdir(case)):
        if f.endswith(".in"):
            inputname=f
            i = i +1
            inlist.append(f)
            print(os.path.join(outdir(case), f))
            
    if i > 1: morethanone = True
    
    return inlist

def interpol(vector_x, vector_y, value_x, der=0):

    tck     = interpolate.splrep(vector_x, vector_y, s=0)
    value_y = interpolate.splev(value_x,tck,der=der)

    return value_y

def nfield_periods(equil, svalue=None, dtheta=3, wrt=False):
    # dtheta (in units of pi) is the range the field
    # flux tube will extent along poloidally. I.e. 
    # A flux tube centered at (theta, zeta) = ( iota*zeta_center , zeta_center) will
    # cover the range iota*zeta_center +/- 3pi poloidally. 
    nfp = read_vmec_var(equil, varname='nfp')
    return nfp * dtheta / iota(equil, svalue)[1]

def euprof(fileprof=None, svalue=None):
    # This function reads a profile in the format of EUTERPE
    # and converts it to another with the parameters
    # that stella needs.
    profdata=loadtxt(fileprof, dtype='float')
    if shape(profdata)[1] == 9:
        S,DTI,TI,DTE,TE,DNI,NI,DNE,NE = arange(0,9)

    dsdrho = 2*sqrt(profdata[:,S])

    s       =  profdata[:,S]
    nine    =  profdata[:,NI]/profdata[:,NE]
    tite    =  profdata[:,TI]/profdata[:,TE]
    tprim_i = -profdata[:,DTI]*dsdrho
    tprim_e = -profdata[:,DTE]*dsdrho
    fprim_i = -profdata[:,DNI]*dsdrho
    fprim_e = -profdata[:,DNE]*dsdrho
    t_i     =  profdata[:,TI]/1000.  # stella uses keV for temp
    t_e     =  profdata[:,TE]/1000.  # stella uses keV for temp
    dens_e  =  profdata[:,NE]/1.0E19 # stella uses keV for temp
    dens_i  =  profdata[:,NI]/1.0E19 # stella uses keV for temp

    if svalue != None:
        print("&vmec_parameters")
        print("torflux="+str(format9(svalue))+'\n')
        print("&parameters")
        print("nine="  +str(format9(interpol(profdata[:,S],nine,svalue))))
        print("tite="  +str(format9(interpol(profdata[:,S],tite,svalue)))+'\n')
        print("&species_parameters")
        print("dens="  +str(format9(interpol(profdata[:,S],dens_i,svalue))))
        print("temp="  +str(format9(interpol(profdata[:,S],t_i,   svalue))))
        print("tprim=" +str(format9(interpol(profdata[:,S],tprim_i[:],svalue))))
        print("fprim=" +str(format9(interpol(profdata[:,S],fprim_i[:],svalue)))+'\n')
        print("\n")
        print("&species_parameters_2")        
        print("dens="  +str(format9(interpol(profdata[:,S],dens_e,svalue))))
        print("temp="  +str(format9(interpol(profdata[:,S],t_e,   svalue))))
        print("tprim=" +str(format9(interpol(profdata[:,S],tprim_e,svalue))))
        print("fprim=" +str(format9(interpol(profdata[:,S],fprim_e,svalue))))
    
##    f       = open('./stella_input_prof.dat',"w")
##    f.write('# (1)psi\t(2)nine\t\t(3)tite\t(4)dens(i)\t(5)temp(i)\t(6)tprim(i)\t(7)fprim(i)'+\
##            '\t(8)dens(e)\t(9)temp(e)\t(10)tprim(e)\t(11)fprim(e)\n')
##    for i in range(0, shape(profdata)[0]):
##        f.write(str(format2(s[i])))       
##        f.write(str(format2(nine[i])))
##        f.write(str(format2(tite[i])))
##        f.write(str(format2(dens_i[i])))
##        f.write(str(format2(t_i[i])))
##        f.write(str(format2(tprim_i[i])))
##        f.write(str(format2(fprim_i[i])))
##        f.write(str(format2(dens_e[i])))
##        f.write(str(format2(t_e[i])))
##        f.write(str(format2(tprim_e[i])))
##        f.write(str(format2(fprim_e[i])+'\n'))

