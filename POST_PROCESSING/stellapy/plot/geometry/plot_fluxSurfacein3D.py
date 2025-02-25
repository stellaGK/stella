#!/usr/bin/python3   
from PIL import Image
import os, sys 
import pathlib
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D #@UnusedImport 

# Stellapy package
sys.path.append(os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent)+os.path.sep)   
from stellapy.data.geometry.get_RZBFromCosineAndSineTransformVMEC import get_RZBFromCosineAndSineTransformVMEC 
from stellapy.data.geometry.get_RZBFromCosineAndSineTransformVMEC import get_CosineAndSineTransformOfRZQ_forFieldLine
from stellapy.utils.files.ensure_dir import ensure_dir 
    
def plot_fluxSurfacein3D(wout_path, svalue=None, rho=0.7, ntheta=50, nzeta=50, quantity="B", \
        set_title=False, colorbar=False, show=True, force_limits=False, vmin=0, vmax=0, plot_field_line=False): 
     
    # Check the quantity: B, jacobian or lambda
    if quantity not in ["B", "J", "L"]: return 
     
    # Radial position
    if svalue==None: svalue = rho*rho
     
    # Read R(zeta,theta); Z(zeta, theta) and Q(zeta, theta) from the VMEC netcdf file
    R, Z, Q, zeta, _, nfp = get_RZBFromCosineAndSineTransformVMEC(wout_path, svalue, quantity, ntheta, nzeta)
    
    # Create the figure
    fig = plt.figure(figsize=(18, 9))
    ax = fig.add_subplot(projection='3d')  
    ax.set_rasterization_zorder(-10)
    
    # Iterate over the field periods
    for i in range(nfp+1):
     
        # Transform (R,Z) to (X,Y) 
        X = R*np.cos(zeta+i*2*np.pi/nfp)
        Y = R*np.sin(zeta+i*2*np.pi/nfp)
        
        # Plot
        ax, m = plot3D(ax, X, Y, Z, Q) 
        
        # Remember limits
        if not force_limits:
            vmin = np.min([vmin,np.min(X)])
            vmax = np.max([vmax,np.max(X)])
      
    # Plot the field line   
    if plot_field_line:   
        R, Z, Q, zeta, theta_pest, nfp = get_CosineAndSineTransformOfRZQ_forFieldLine(wout_path, svalue, quantity, nzeta=50, toroidal_turns=1, alpha0=0)
        X = R*np.cos(zeta+i*2*np.pi/nfp)
        Y = R*np.sin(zeta+i*2*np.pi/nfp)
        ax.plot(X, Y, Z, color='red', lw=10) 
        
    # Set limits
    ax.set_xlim([vmin, vmax])
    ax.set_ylim([vmin, vmax])    
    ax.set_zlim3d([vmin, vmax])
       
    # Set title
    if set_title:
        if quantity=="B": ax.set_title('Magnetic field strength')
        if quantity=="L": ax.set_title('Lambda')
        if quantity=="J": ax.set_title('Jacobian')
       
    # Intial viewing angle
    view = "narrow"; angle=5
    if view=="top": ax.view_init(90,angle)
    if view=="bottom": ax.view_init(270,angle)
    if view=="side": ax.view_init(0,angle)
    if view=="both": ax.view_init(45,angle)
    if view=="narrow": ax.view_init(15,angle)
       
    # Show figure 
    ax.set_box_aspect((np.ptp([vmin, vmax]), np.ptp([vmin, vmax]), np.ptp([vmin, vmax])))  
    if colorbar: fig.colorbar(m)
    ax.set_axis_off()
    if show: plt.show()
    return fig, ax

#---------------------
def plot3D(ax,x,y,z,f):
      
    # fourth dimention - colormap
    # create colormap ccording to x-value (can use any 50x50 array)
    color_dimension = f # change to desired fourth dimension
    minn, maxx = color_dimension.min(), color_dimension.max()
    norm = matplotlib.colors.Normalize(minn, maxx)
    m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
    m.set_array([])
    fcolors = m.to_rgba(color_dimension)
    
    # plot 
    ax.plot_surface(x,y,z, rstride=1, cstride=1, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')  
    return ax, m


#-------------------------------------------------------------------------------
#                        Run as the main script
#-------------------------------------------------------------------------------
if __name__ == "__main__":

    # VMEC file
    stella_folder = os.path.abspath(pathlib.Path(os.environ.get('STELLAPY')).parent.parent)
    figure_folder = stella_folder / pathlib.Path("stellapy_figures"); ensure_dir(figure_folder)
    vmec = stella_folder / pathlib.Path("DOCUMENTATION/VMEC_equilibria/wout_w7x_standardConfig.nc"); device = 'W7X'
    vmec = stella_folder / pathlib.Path("DOCUMENTATION/VMEC_equilibria/wout_161s1.nc"); device = 'NCSX'
    
    # Radial location
    rho = 0.8
    
    # Resolution
    n = 50
    dpi = 500
            
    # Plot the last closed flux surface     
    fig, ax = plot_fluxSurfacein3D(vmec, rho=rho, ntheta=n, nzeta=n, quantity="B", show=False)
     
    # Save the figure as eps 
    name = figure_folder/("rho"+str(int(rho*100))+"_n"+str(n)+"_"+device+".png")
    fig.savefig(name, format='png', dpi=dpi)
    print("Saved 3D flux surface for "+vmec.name)  
                
