# -*- coding: utf-8 -*-
import numpy as np
import f90nml
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
#import matplotlib.colors as colors
# setup some plot defaults
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.rc('font', size=30)
rcParams.update({'text.latex.preamble' : r'\usepackage{bm}'})
rcParams.update({'figure.autolayout': True})

def read_variable(stelladata,varstring):
    try:
        var=stelladata[varstring].data
        var_present=True
    except KeyError:
            print('INFO: '+varstring+' not found in netcdf file')
            var=None
            var_present = False
    return var, var_present
    
def dimension_sizes(stelladata):
    nky = stelladata["ky"].size
    nkx = stelladata["kx"].size
    nzed = stelladata["zed"].size
    nvpa = stelladata["vpa"].size
    nmu = stelladata["mu"].size
    nspec = stelladata["species"].size
    ntubes = stelladata["tube"].size
    nalpha = stelladata["alpha"].size
    return nky, nkx, nzed, nvpa, nmu, nspec, ntubes, nalpha
    
def flux_data(stelladata):
    pflx = stelladata["pflx"].data
    vflx = stelladata["vflx"].data
    qflx = stelladata["qflx"].data
    time = stelladata["t"].data
    return pflx, vflx, qflx, time
    
def field_data(stelladata):
    phi2 = stelladata["phi2"].data
    apar2, apar2_present = read_variable(stelladata,"apar2")
    bpar2, bpar2_present = read_variable(stelladata,"bpar2")
    return phi2, apar2, apar2_present, bpar2, bpar2_present
    
def species_data(stelladata,filename):
    nspec = stelladata["species"].size
    
    charge = stelladata["charge"].data
    mass = stelladata["mass"].data
    dens = stelladata["dens"].data
    temp = stelladata["temp"].data
    tprim = stelladata["tprim"].data
    fprim = stelladata["fprim"].data
    vnew = stelladata["vnew"].data
    typeint = stelladata["type_of_species"].data
    typestring = get_species_string(filename)
    return charge, mass, dens, temp, tprim, fprim, vnew, typeint, typestring

def get_species_string(filename):
    with open(filename+'.in') as nml_file:
        namelist = f90nml.read(nml_file)
    
    nspec = namelist["species_knobs"]["nspec"]
    
    species_string = []
    species_string_plot = []
    ion_counter = 0
    for ispec in range(0,nspec):
        type = namelist["species_parameters_"+str(ispec+1)]["type"]
        if type == "ion":
            ion_counter += 1
            species_string_plot.append(type[0]+str(ion_counter))
        else:
            species_string_plot.append(type[0])
        species_string.append(type)
    #print(species_string)
    #print(species_string_plot)
    
    return species_string_plot
        
        
def plot_1d_list_pdf (xlist,ylist,marker_list,xlab, pdf,
  title='',ylab='',xlims=None,ylims=None,aspx=12,aspy=8, xticks = None, yticks = None,
  markersize=5, legend_title="", use_legend=False,loc_opt='upper right', ylab_list = None,
  bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=10, ncol_opt=1,
  legend_shadow=False,legend_frame=False, vlines = None,marker_fill_style = None,
  cartoon=False, linewidth=None, texts = None):

    fig=plt.figure(figsize=(aspx,aspy))
    nlist = len(ylist)
    if(ylab_list is None):
        ylab_list = [None for i in range(0,nlist)]
    for iy in range(0,nlist):
        plt.plot(xlist[iy],ylist[iy],marker_list[iy],markersize=markersize,label=ylab_list[iy],
        fillstyle = marker_fill_style, linewidth = linewidth)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    if(not xlims is None):
        plt.xlim(xlims[0],xlims[1])
    if(not ylims is None):
        plt.ylim(ylims[0],ylims[1])
    if(not vlines is None):
        for xin,xlabel,xcolor,xlinestyle in vlines:
            plt.axvline(x=xin, label=xlabel, color=xcolor,linestyle=xlinestyle,linewidth=linewidth)   
    if(not texts is None):
        for xin, yin, textin in texts:    
            print(xin,yin,textin)
            plt.text(xin,yin,textin)
    if(use_legend):
        plt.legend(title=legend_title,loc=loc_opt, bbox_to_anchor=bbox_to_anchor_opt,
        fontsize=legend_fontsize, frameon=legend_frame, handlelength=1, labelspacing=0.5,
        ncol=ncol_opt, columnspacing = 0.5 , handletextpad = 0.5, shadow=legend_shadow)
    if(not xticks is None):
        plt.xticks(xticks)
    if(not yticks is None):
        plt.yticks(yticks)    
    if (cartoon):
        plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='off')
        plt.box(False)
    pdf.savefig(fig)# pdf is the object of the current open PDF file to which the figures are appended
    plt.close (fig)
    return

def plot_fluxes_nc(filename, pflx, vflx, qflx, time, typestring):
    nspec = len(typestring)
    file = filename+".fluxes.pdf"
    pdf = PdfPages(file)
    marker_list = ['k','b','r','g','c','y']
    xlist = [time for ispec in range(0,nspec)]
    xlab = "$ t / (a /v_{\\rm{th,ref}}) $"

    # first the heat flux

    ylist = [qflx[:,ispec] for ispec in range(0,nspec)]
    ylab_list = [ "$Q_{" + typestring[ispec]+"}$" for ispec in range(0,nspec)]
    plot_1d_list_pdf (xlist,ylist,marker_list,xlab, pdf,
      title='',ylab='',xlims=None,ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
      markersize=5, legend_title="", use_legend=True,loc_opt='upper right', ylab_list = ylab_list,
      bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=20, ncol_opt=1,
      legend_shadow=False,legend_frame=False, vlines = None,marker_fill_style = None,
      cartoon=False, linewidth=None, texts = None)

    # now the particle flux  
    ylist = [pflx[:,ispec] for ispec in range(0,nspec)]
    ylab_list = [ "$\\Gamma_{" + typestring[ispec]+"}$" for ispec in range(0,nspec)]
    plot_1d_list_pdf (xlist,ylist,marker_list,xlab, pdf,
      title='',ylab='',xlims=None,ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
      markersize=5, legend_title="", use_legend=True,loc_opt='upper right', ylab_list = ylab_list,
      bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=20, ncol_opt=1,
      legend_shadow=False,legend_frame=False, vlines = None,marker_fill_style = None,
      cartoon=False, linewidth=None, texts = None)
      
    # now the momentum flux  
    ylist = [vflx[:,ispec] for ispec in range(0,nspec)]
    ylab_list = [ "$\\Pi_{" + typestring[ispec]+"}$" for ispec in range(0,nspec)]
    plot_1d_list_pdf (xlist,ylist,marker_list,xlab, pdf,
      title='',ylab='',xlims=None,ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
      markersize=5, legend_title="", use_legend=True,loc_opt='upper right', ylab_list = ylab_list,
      bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=20, ncol_opt=1,
      legend_shadow=False,legend_frame=False, vlines = None,marker_fill_style = None,
      cartoon=False, linewidth=None, texts = None)

    pdf.close()
    print(file)
    return None
    
def plot_fields_nc(filename,  phi2, apar2, apar2_present, bpar2, bpar2_present, time):
    file = filename+".fields.pdf"
    pdf = PdfPages(file)
    marker_list = ['k','b','r','g','c','y']
    xlist = [time]
    xlab = "$ t / (a /v_{\\rm{th,ref}}) $"

    
    ylist = [phi2]
    ylab_list = [ "$ \\langle(e\\Phi/ \\rho_\\ast T_{\\rm ref})^2\\rangle $"]
    if apar2_present:
        xlist.append(time)
        ylist.append(apar2)
        ylab_list.append("$ \\langle(A_{\\|} / \\rho_{\\ast} \\rho_{\\rm ref} B_{\\rm ref})^2\\rangle $")
    if bpar2_present:
        xlist.append(time)
        ylist.append(bpar2)
        ylab_list.append("$ \\langle(B_{\\|} / \\rho_{\\ast} B_{\\rm ref})^2\\rangle $")
    
    plot_1d_list_pdf (xlist,ylist,marker_list,xlab, pdf,
      title='',ylab='',xlims=None,ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
      markersize=5, legend_title="", use_legend=True,loc_opt='upper right', ylab_list = ylab_list,
      bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=20, ncol_opt=1,
      legend_shadow=False,legend_frame=False, vlines = None,marker_fill_style = None,
      cartoon=False, linewidth=None, texts = None)

    pdf.close()
    print(file)
    return None