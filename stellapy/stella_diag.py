from stella_dirs import *
from stella_utils import *
from numpy import *
from pylab import *
from struct import *
from scipy import *
from scipy.special import expit
from matplotlib import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.io import netcdf
import os.path
import h5py

plt.rcParams.update({'font.size': 28})
plt.rcParams['lines.linewidth'] = 2

import tabCompleter

from plotbox import *
from tabCompleter import *
import struct
import physcon as pc
import time
from os import listdir
from stella_read import *



def kspectra_movie(case, ):
    
    from stella_plots import movie_2d

    phi2_vs_kx_ky, k_x, k_y, n_time =\
                  phi2_vs_kxky(case), kx(case)[0], ky(case)[0], time(case)[1]
    
    phi2max = np.arange(n_time,dtype=float)
    phi2min = np.arange(n_time,dtype=float)
    
    for i in range(n_time):    
        phi2max[i] = np.absolute(phi2_vs_kx_ky[i,:,:].max())*rescale

    phi2min[:] = 0.0
    ylabel = '$k_x$'
    xlabel = '$k_y$'
    title  = '$|\\varphi(k_x, k_y)|^2$'

    movie_file = outdir(case) + '/phi2_vs_kxky.mp4'

    movie_2d(phi2_vs_kx_ky, k_y, k_x, phi2min,\
             phi2max, n_time-1, movie_file, xlabel, ylabel, title,cmp='YlGnBu')

def phi_t(case, kx_idx=None, ky_idx=None, t0=None, tube_idx=0, last=False):
    
    k_x, nakx, nakx_mid  = kx(case)
    k_y, naky            = ky(case)
    zeta, nzed, nzed_mid = zed(case)
    time_trace           = time(case)[0]
    
    phi_data = phi_vs_t(case)
    phi_data = np.concatenate((phi_data[:,:,:,nakx_mid:,:,:],\
                               phi_data[:,:,:,:nakx_mid,:,:]),axis=3)
    
    if kx_idx == None: kx_idx = 0
    if ky_idx == None: ky_idx = 0

    # shape(phi_vs_t(case))=(N_time, N_tubes, N_zed, N_kx, N_ky, 2)
    # Warning: check why the last dimension is 2, if phi[:,:,:,:,:,0]=0.0 !!!
    phi_vs_z_time_real = phi_data[:,tube_idx, :, kx_idx, ky_idx, 0]
    phi_vs_z_time_imag = phi_data[:,tube_idx, :, kx_idx, ky_idx, 1]


    print("Getting phi(z,t) for (kx, ky) = (", k_x[kx_idx], ',', k_y[ky_idx], ')')
    
    # The index for the with the nearest kx and ky
    if last == True:
        phi_vs_z_real_to_plot = phi_vs_z_time_real[shape(phi_vs_z_time_real)[1]-1, :]
        phi_vs_z_imag_to_plot = phi_vs_z_time_imag[shape(phi_vs_z_time_imag)[1]-1, :]
        time0                 = time_trace[shape(time_trace)[0]-1]

    if t0 != None:
        phi_vs_z_real_to_cut  = phi_vs_z_time_real[(time_trace >= t0),:]
        phi_vs_z_imag_to_cut  = phi_vs_z_time_imag[(time_trace >= t0),:]
        phi_vs_z_real_to_plot = phi_vs_z_real_to_cut[0,:]
        phi_vs_z_imag_to_plot = phi_vs_z_imag_to_cut[0,:]
        time_trace            = time_trace[(time_trace >= t0)]
        time0                 = time_trace[0]


    maxre  = phi_vs_z_real_to_plot.max()
    maxim  = phi_vs_z_imag_to_plot.max()
    maxphi = sqrt(phi_vs_z_real_to_plot**2.0+phi_vs_z_imag_to_plot**2.0).max()
        

    # Plotting |ph|
    
    ax = pl2d(xrange=[-pi,pi], yrange=[-maxphi,maxphi], xlabel='$\\zeta$', ylabel='$\\varphi^2$',\
         fig_size=(8.5, 7.5), title='t = '+str(time0) + ', (kx, ky) = (' +\
              str(format3(k_x[kx_idx])) + ', ' + str(format3(k_y[ky_idx])) + ')', ax=None, log=False)
    ax.plot(zeta, phi_vs_z_real_to_plot**2+phi_vs_z_imag_to_plot**2, '-',\
            color='black', linewidth=2, label='$\\varphi^2$')
    ax.autoscale()

    # Plotting Re(phi), Im(phi)
    ax = pl2d(xrange=[-pi,pi], yrange=[-maxre,maxre], xlabel='$\\zeta$', ylabel='$\\Re(\\varphi), \\Im({\\varphi})$',\
              fig_size=(8.5, 7.5), title='t = '+str(time0) + ', (kx, ky) = (' +\
              str(format3(k_x[kx_idx])) + ', ' + str(format3(k_y[ky_idx])) + ')',\
              ax=None, log=False)
    ax.plot(zeta, phi_vs_z_real_to_plot, '-',\
            color='black', linewidth=2, label='$\\Re({\\varphi})$')
    ax.legend(loc=1,labelspacing=0.0, prop={'size':24})
    
    ax2 = ax.twinx()
    ax2.plot(zeta, phi_vs_z_imag_to_plot, '-', color='blue', linewidth=2, label='$\\Im(\\varphi)$')
    ax2.set_ylim([-maxim, maxim])
    ax2.legend(loc=2,labelspacing=0.0, prop={'size':24}) 
#    ax2.autoscale()    
    show()
   

        

def phi2_kx_ky(case, trange=None, crange=None, log=True,\
               delta=1E-12, t0=None, last=False, movie=False):

    from stella_plots import movie_2d

    k_x, nakx, nakx_mid = kx(case)
    
    phi2_vs_kx_ky, k_y, n_time =\
                  phi2_vs_kxky(case), ky(case)[0], time(case)[1]

    phi2_vs_kx_ky = np.concatenate((phi2_vs_kx_ky[:, nakx_mid:,:],\
                                    phi2_vs_kx_ky[:,:nakx_mid ,:]),axis=1)
    time_trace = time(case)[0]

    if trange != None:
        # phi2_vs_kx_ky has dimensions (ntime, nakx, naky)
        # Here, we remove all values of phi2_vs_kx_ky of times out of the selected
        # interval.
        phi2_vs_kx_ky_to_avrg = phi2_vs_kx_ky[(time_trace > trange[0]) & (time_trace < trange[1]),:,:]
    else:
        phi2_vs_kx_ky_to_avrg = phi2_vs_kx_ky

    if last == True:
        phi2_vs_kx_ky_avrg = phi2_vs_kx_ky_to_avrg[shape(phi2_vs_kx_ky_to_avrg)[0]-1,:,:]
    else:
        phi2_vs_kx_ky_avrg = mean(phi2_vs_kx_ky_to_avrg, axis=0)

    if t0 != None:
        phi2_vs_kx_ky_to_avrg = phi2_vs_kx_ky[(time_trace >= t0),:,:]
        time_trace            = time_trace[(time_trace >= t0)]
        phi2_vs_kx_ky_avrg    = phi2_vs_kx_ky_to_avrg[0,:,:]

    ylabel = '$k_x\\rho_r$'
    xlabel = '$k_y\\rho_r$'
    title  = '$\\left<\\varphi^{2}\\right>(k_x, k_y)$ ' + 't = ' + str(time_trace[0])

    if crange != None:
        zmin = crange[0]
        zmax = crange[1]
    else:
        zmin = max(phi2_vs_kx_ky.min(), delta)
        zmax = phi2_vs_kx_ky.max()

    if movie:
        outfile = outdir(case) + '/phi2_vs_kxky.mp4'
        cmp='jet'
        step = 1
        fig = plt.figure(figsize=(12,8))
        x,y = np.meshgrid(k_y,k_x)

        
        im = plt.imshow(phi2_vs_kx_ky[0,:,:], cmap=cmp, vmin=zmin, vmax=zmax,
                        extent=[x.min(),x.max(),y.min(),y.max()],
                        interpolation='nearest', origin='lower', aspect='auto',\
                        norm=LogNorm())
        plt.colorbar()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
#        plt.title(title)

        ims = []
        ims.append([im])

        for i in range(1,n_time-1,step):
            im = plt.imshow(phi2_vs_kx_ky[i,:,:], cmap=cmp, vmin=zmin, vmax=zmax,
                            extent=[x.min(),x.max(),y.min(),y.max()],
                            interpolation='nearest', origin='lower', aspect='auto',\
                            norm=LogNorm())
            #title  = '$\\left<\\varphi^{2}\\right>(k_x, k_y)$'+ 't = ' + str(time_trace[i])
#            im.title(title)
            ims.append([im])
            
        ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True)
        ani.save(outfile)
        
    else:
        yran = [phi2_vs_kx_ky_avrg.min(), phi2_vs_kx_ky_avrg.max()]
        ax = cmap(xdata=k_y, ydata=k_x, zdata=phi2_vs_kx_ky_avrg,\
                  xlabel=xlabel, ylabel=ylabel, zlabel=title,\
                  ctics=False, contour=1, color=True, cmap='jet',\
                  epsname=None, crange=crange, fig_size=(8.5, 7.5))
        show()


def get_kperp2(case):
    # This function builds kperp2 as done internally by the code
    # See dist_fn.f90
    #
    # do iky = 1, naky
    #   if (zonal_mode(iky)) then
    #      do ikx = 1, nakx
    #         kperp2(iky,ikx,:,:) = akx(ikx)*akx(ikx)*gds22/(geo_surf%shat**2)
    #      end do
    #   else
    #      do ikx = 1, nakx
    #         kperp2(iky,ikx,:,:) = aky(iky)*aky(iky) &
    #              *(gds2 + 2.0*theta0(iky,ikx)*gds21 &
    #              + theta0(iky,ikx)*theta0(iky,ikx)*gds22)
    #      end do
    #   end if
    # end do
    kx                = kx_stella(case)
    ky                = ky_stella(case)
    z, nzed, mid_nzed = zed(case)
    
    kperp2 = empty((size))



def omega_k(case, last=False, view=True, write=False):
    # case should be a list of input files.
    for i in case:
        n_ky, n_kx = ky(i)[1], n_kx = kx(i)[1]

        if n_ky == 1 and n_kx == 1:
            omega_file = i+'.omega'
            omega_data = loadtxt(file_omega, dtype='float')
            omega_data_finite = omega_data.dropna()
#            omega_data_last   = omega_data_finite[]
            
            

def omega(case, last=False, view=True, yrange=None, xrange=None):
    # Note that:
    # time   ky   kx   Re[om]   Im[om]   Re[omavg]  Im[omavg]
    
    for i in arange(0,size(case)):
        datafile = outfile(case[i],'omega')
        n_time   = time(case[i])[1] # Including t=0
        n_ky     = ky(case[i])[1]
        n_kx     = kx(case[i])[1]

        data     = loadtxt(datafile, dtype='float').reshape(n_time-1, n_kx, n_ky, 7)

        if view:
            pl2y(ky(case[i])[0], y1data=data[n_time - 2, 0,:,3], y2data=data[n_time - 2, 0,:,4],\
                 xlabel='$k_y\\rho_{i}$', ylabel='$\\omega a/v_{\mathrm{th},\mathrm{i}}$',\
                 key1='$a\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                 key2='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$', xrange=[0,20],\
                 yrange=yrange,fig_size=(8.5, 7.5), wp=1, ax=None, mkt="o", mkc='red',\
                 title="$k_x="+str(kx(case[i])[0][0])+"$",\
                 hline1=None, hline2=None, vshadow=None, ls1=1, ls2=2)
            show()
            
        if n_kx == 1 and n_ky > 1:
                # The evolution of omega(ky, t) can be represented
                t      = time(case[i])[0][1:]
                kyrhoi = ky(case[i])[0]
                ky_re  = data[:, 0,:,3]
                ky_im  = data[:, 0,:,4]
                [X,Y]  = meshgrid(t, kyrhoi)
                Z      = ky_im.T

                surf(X, Y, Z, xlabel='$t$', ylabel='$k_y\\rho_{i}$',\
                     zlabel='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                     title="$k_x="+str(kx(case[i])[0][0])+"$", zrange=yrange, ax=None)
                show()

        if n_kx == 1 and n_ky == 1:
                # The evolution of omega(t) can be represented
                t      = time(case[i])[0][1:]
                kyrhoi = ky(case[i])[0]
                om_re  = data[:, 0,0,3]
                notnan = logical_not(isnan(om_re))
                om_re  = om_re[notnan]
                om_im  = data[:, 0,0,4]
                om_im  = om_im[notnan]
                t      = t[notnan]
                
                xlims = [t.min(), t.max()]
                ylims = [min(om_im.min(), abs(om_re).min())*1.05, max(om_im.max(), abs(om_re).max())*1.05]

                if view: 
                    pl2y(t, abs(om_re), om_im, xlabel='$t$', ylabel='$\\omega a/v_{\mathrm{th},\mathrm{i}}$',\
                         title="$(k_x, k_y)=("+str(kx(case[i])[0][0])+", " + str(ky(case[i])[0][0])+")$",\
                         xrange=xlims, yrange=ylims,\
                         key1='$a\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                         key2='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                         ax=None)

                

                    show()
            
    if last:
        return data[n_time-2, :, :, :]
    else:
        return data

def omega_all(case, last=False, view=True, yrange=None, xrange=None, scan='ky'):
    # This function takes a folder where various single (kx,ky)
    # runs have been performed and gather all results in one file,
    # taking the last instant with finite values of omega.
#    hidden = casestr(case).startswith('.')

    input_all=[outdir(case) + '/' + i for i in inputlist(case)]

    for i in input_all:
        # We assume that if no omega file related to a certain input is present
        # then the case is not valid.
        #print(outfile(case + '/' + i, quant='omega'))
        if not os.path.isfile(outfile(i, quant='omega')):
            print("No *omega data found for case with input " + i +'. Case ignored.') 
            input_all.remove(i)
            
    nval            = size(input_all)
    ns              = nspec([i])
    omega_all_data  = empty((nval, 7))
    phi2_all_data   = empty(nval)

    
    for i in arange(0, nval-1):
        omega_file = outfile(input_all[i], quant='omega')
        omega_all_data[i,:]=omega_last(omega_file)
        phi2_all_data[i]=phi2_vs_kxky(input_all[i])[shape(phi2_vs_kxky(input_all[i]))[0]-1,0,0]
        print(phi2_all_data[i])

    
    

def omega_last(omegafile):
    # This function returns the last row with finite values of omega
    # of a certaing omega file. The full path must be pass onto this function.
    # It assumes that only pair (kx,ky) are considered
    omegadata   = loadtxt(omegafile, dtype='float')
    notnan      = isfinite(omegadata[:,3])
    omegadata_f = omegadata[notnan,:]
    return omegadata_f[shape(omegadata_f)[0]-1,:]
        
        


def multi_lin_data(case=None, view=True, yrange=None, xrange=None, sort='ky',\
                   off=[], save=True, ql='per', last=True):
    # Note that:
    # ql == 'mik' applies formula (1) of Mikkelen et al. PoP 21 082302 (2014)
    # ql == 'per' applies formula (13)-(14) of Helander and Zocco PPCF 60 084006 (2018)
    
    # time   ky   kx   Re[om]   Im[om]   Re[omavg]  Im[omavg]

    if size(case) == 1:
        # When a directory is pass onto "case", several runs of single kx,ky
        # are assumed. The input files on that directory are read, neglecting
        # those with no *out.nc output.
        cases=[case[0] + '/' + i for i in inputlist(case[0])]
    else:
        # When a list of input files of interest is pass onto, we take
        # it as a set of cases to run the diagnostics ql_fluxes over.
        cases = case

    # We count first the number of input files with existing *out.nc
    count=0
    not_valid=[]

    print("\nChecking existence of *out.nc file corresponding inputs.")

    for i in cases:
        # Analysis is performed only if *out.nc file is present
        if os.path.isfile(outfile(i, 'out.nc')):           
            count = count + 1
            ns    = nspec(i)
        else:
            print("Case ", i, " removed from input list.")
            not_valid.append(i)
            
        for j in off:
            if j in i: not_valid.append(i)
            
    print("Number of found input files with corresponding output *.nc =", count)
            
    for j in not_valid:
        cases.remove(j)

    # Columns will contain kxrhoi, kyrhoi, re(omega) im(omega) nspec*3 fluxes
    ncol = 7 + ns*3
    lin_data = empty((size(cases), ncol), dtype='float')

    for i in cases:
        # Analysis is performed only if *out.nc file is present
        print("ql_fluxes diagnostics running over case", i)
        lin_data[cases.index(i),:] = ql_fluxes(case=i, view=False, formula=ql)[1]
        if lin_data[cases.index(i),6] < 100:
            print("Warning! phi2 found too low: setting omega, gamma and fluxes to NaN.")
            lin_data[cases.index(i),3]  = NaN
            lin_data[cases.index(i),2]  = NaN
            lin_data[cases.index(i),7:] = NaN            

    # We sort the data by increasing ky order
    olab   = '$\\omega a/v_{\\mathrm{th},i}$'
    glab   = '$\\gamma a/v_{\\mathrm{th},i}$'
    qlGlab = '$\\Gamma^{\\mathrm{ql}}$ [a.u.]'
    qlUlab = '$U^{\\mathrm{ql}}$ [a.u.]'
    qlQlab = '$Q^{\\mathrm{ql}}$ [a.u.]'
   
    
    if sort == 'ky':
        col = 1 ; xlab = '$k_{y}\\rho_i$'
    if sort == 'kx':
        col = 0 ; xlab = '$k_{x}\\rho_i$'
        
    lin_data = lin_data[lin_data[:,col].argsort()]

    if save == True:
        fn = outdir(cases[0]) + '/multi_lin_data.h5'
        hf = h5py.File(fn, 'w')
        hf.create_dataset('kx', data=lin_data[:,0])
        hf.create_dataset('ky', data=lin_data[:,1])
        hf.create_dataset('omega_ql_last', data=lin_data[:,2])
        hf.create_dataset('gamma_ql_last', data=lin_data[:,3])
        hf.create_dataset('omega_full_last', data=lin_data[:,4])
        hf.create_dataset('gamma_full_last', data=lin_data[:,5])        
        hf.create_dataset('phi2_ql_last', data=lin_data[:,6])
        hf.create_dataset('qlpflx', data=lin_data[:,7:7+ns])
        hf.create_dataset('qlvflx', data=lin_data[:,7+ns:7+2*ns])
        hf.create_dataset('qlqflx', data=lin_data[:,7+2*ns:7+3*ns])
        hf.create_dataset('vth_r', data=lin_data[:,7+2*ns:7+3*ns])
        hf.close()
        print("File saved: ", fn)

    if view == True:

        # Plotting
        xran = [min(lin_data[:,col]), max(lin_data[:,col])]

        # omega
        ax1 = pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=olab,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax1.autoscale()
        ax1.plot(lin_data[:,col], lin_data[:,2], 'o-', color='navy',\
                 linewidth=4, markerfacecolor='white')
        if last == True:
            ax1.plot(lin_data[:,col], lin_data[:,4], 'o-', color='gray',\
                     linewidth=2, markerfacecolor='white')
        # gamma
        ax2 = pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=glab,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax2.autoscale()
        ax2.plot(lin_data[:,col], lin_data[:,3], 'o-', color='crimson',\
                 linewidth=4, markerfacecolor='white')
        if last == True:
            ax2.plot(lin_data[:,col], lin_data[:,5], 'o-', color='gray',\
                     linewidth=2, markerfacecolor='white')            

    
        for k in arange(0, ns):
            fig = plt.figure(figsize=(18, 9))
            fig.subplots_adjust(left=0.1, wspace=0.35)

            title = 'species '+str(k+1)
            ax1 = fig.add_subplot(131)
            pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=qlGlab,\
                 fig_size=(8.5, 7.5), title=title, ax=ax1)
            ax1.autoscale()
            ax1.plot(lin_data[:,col], lin_data[:,7+k], 'o-', color='darkviolet',\
                     linewidth=4, markerfacecolor='white')
            
            ax2 = fig.add_subplot(132)
            pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=qlUlab,\
                 fig_size=(8.5, 7.5), title=title, ax=ax2)
            ax2.autoscale()
            ax2.plot(lin_data[:,col], lin_data[:,7+ns+k], 'o-', color='seagreen',\
                     linewidth=4, markerfacecolor='white')

            ax3 = fig.add_subplot(133)
            pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=qlQlab,\
                 fig_size=(8.5, 7.5), title=title, ax=ax3)
            ax3.autoscale()
            ax3.plot(lin_data[:,col], lin_data[:,7+2*ns+k], 'o-', color='darkorange',\
                     linewidth=4, markerfacecolor='white')
            
        show()


#            print(type(ks), type(fluxes))
            

                                 

#            multi_lin_data  = concatenate([array(fluxes[0],fluxes[1]), fluxes[2]])
#            print(multi_lin_data)
#            omega  = 

def fluxes(case=None, plot=False, si=False, time_avrg=None, tref=None, sign=1):
    i        = case
    fluxfile = outfile(i,'fluxes')
    fluxdata = loadtxt(fluxfile, dtype='float')
    nspecies = species(i)[1]
    
    if si:
        time_ref      = 1/ref_values(case, tref=tref)[5]
        fluxes_ref    = array(ref_values(case, tref=tref)[9:12])

        # We undo the normalization to get SI units
        fluxdata[:,0]=fluxdata[:,0]*time_ref
        for j in arange(0,3):
            for k in arange(nspecies):
                fluxdata[:,j+1+k*3]=sign*fluxdata[:,j+1+k*3]*fluxes_ref[j]
            
    if time_avrg :
        fluxes_avrg = time_average(fluxdata, interval=time_avrg)
    
    if plot:
        if si:
            xlab, ylab = '$t$ [$s$]',\
                         ['$\\Gamma_s$ [m$^{-2}$s$^{-1}$] ',\
                          '$\\Pi_s$ [kg s$^{-1}$]', '$Q_s$ [W m$^{-2}$]']
        else:
            xlab, ylab = '$t$ [$\\Omega_{r}^{-1}$]',\
                         ['$\\Gamma_s/\\Gamma_{gBs}$', '$\\Pi_s/\\Pi_{gBs}$', '$Q_s/Q_{gBs}$']
            
        for j in arange(0,3):
            ax = pl2d(xrange=None, yrange=None, xlabel=xlab, ylabel=ylab[j],\
                      fig_size=(8.5, 7.5))
            
            for k in arange(nspecies):
                xran = [min(fluxdata[:,0]), max(fluxdata[:,0])]
                yran = [min(fluxdata[:,j+1+k*3])*1.05, max(fluxdata[:,j+1+k*3])*1.05]
                
                ax.plot(fluxdata[:,0], fluxdata[:,j+1+k*3],\
                        '-', color='black', linewidth=3, label='species = ' + str(k))
                ax.set_xlim(xran)
                ax.set_ylim(yran)
                ax.ticklabel_format(style='sci', scilimits=(0,0))
                if time_avrg:
                    ax.axhline(y=fluxes_avrg[j+1+k*3],\
                               linestyle='--', linewidth=2, color='blue')

                    ax.barh(y=yran[0], width=(time_avrg[1]-time_avrg[0]),\
                            height=yran[1]-yran[0], \
                            left=time_avrg[0], align='edge', facecolor='yellow', alpha=0.3)
                    
            ax.legend(loc='best',labelspacing=0.0, prop={'size':26})
            
        show()

    if time_avrg : return(fluxes_avrg)
    else: return(fluxdata)


def time_average(data, time_column=0, interval=None):
    # Assuming the data is stored in an array such that
    # time is stored at data[:,time_average], this function averages
    # in in interval [t0,tf] the data store in data[:,j]
    # with j different from time_average.
    #
    # interval := [t0, tf]
    #
    data_avrg=empty(size(data[0,:]), dtype='float')
    data_avrg[:]=0
    num_val=0

    
    for i in arange(0, size(data[:, time_column])):
        if data[i, time_column] > interval[0] and data[i, time_column] < interval[1]:
            data_avrg[:]=data_avrg[:]+data[i,:]
            num_val = num_val+1

    return data_avrg/num_val
        



def ql_fluxes(case=None, view=False, last=True, formula='per'):
    i        = case
    tfull    = time(i)[0][1:]
    omegafile= outfile(i,'omega')
    fluxfile = outfile(i,'fluxes')
    n_time   = time(i)[1] # Including t=0
    n_ky     = ky(i)[1]
    n_kx     = kx(i)[1]
    ns       = nspec(i)
    density  = dens(i)[0]/dens(i)[0][0]
    # Full vectores to calculate the flux 
    flux    = empty((ns*3,size(tfull)), dtype='float')
    flux_ql = empty((ns*3,size(tfull)), dtype='float')
    flux_ql_last = empty(ns*3, dtype='float')
        
    fluxcols = ns*3+1

    if n_ky != 1 or n_kx != 1:
        print("Error: to use this function n_ky or n_kx cannot be > 1."); return

    omegadata = loadtxt(omegafile, dtype='float').reshape(n_time-1, n_kx, n_ky, 7)
    omega     = omegadata[:,0,0,3] 
    gamma     = omegadata[:,0,0,4]
    # First line of fluxdata is at time=0.0, omegadata begins after first time step
    print(fluxfile)
    fluxdata  = loadtxt(fluxfile, dtype='float')[1:,:]
    phi2data  = phi2_vs_kxky(i)[1:,0,0]
    kxrhoi    = kx(i)[0][0]
    kyrhoi    = ky(i)[0][0]
    
    # For plotting the particle fluxes
    if view:
        ax1 = pl2d(xrange=None, yrange=None, xlabel="$t v_{th,i}/a$",\
                   ylabel="$\\Gamma/\\Gamma_{gB}$", fig_size=(8.5, 7.5),\
                   title="$(k_x, k_y)\\rho_i=("+str(kx(i)[0][0])+", " + str(ky(i)[0][0])+")$",\
                   log=False)
        line = ['-', '--', ':']
        color = ['blue', 'crimson', 'green']

    for i in arange(0, 3*ns):
        flux[i,:]    = fluxdata[:,i+1]
        dens_i       = i-int(i/ns)*ns
        if formula == 'mik':
            flux_ql[i,:] = flux[i,:]*gamma/phi2data/density[dens_i]/kyrhoi**2.0
        elif formula == 'per':
            flux_ql[i,:] = flux[i,:]*(gamma**2.0+omega**2.0)/phi2data/density[dens_i]/kyrhoi**2.0/gamma
        # We remove not finite values
        
        notnan_flux      = isfinite(flux_ql[i,:])
        notnan_omega     = isfinite(omega[:])

        # We prepare all quantities to build the QL fluxes
        # up to the last value where these are finite
        # It can occur that this happens before than the instant
        # where omega finds the same problem.
        t_f        = tfull[notnan_flux]
        flux_ql_f  = flux_ql[i, notnan_flux]
        flux_f     = flux[i, notnan_flux]
        omega_ql   = omega[notnan_flux]
        gamma_ql   = gamma[notnan_flux]
        phi2data_f = phi2data[notnan_flux]

        omega_full = omega[notnan_omega]
        gamma_full = gamma[notnan_omega]
        omega_full_last = omega_full[size(omega_full)-1]
        gamma_full_last = gamma_full[size(gamma_full)-1]
        
        t_last_pos = size(t_f)-2
        
        t_last          = t_f[t_last_pos]
        flux_ql_last[i] = flux_ql_f[t_last_pos]
        omega_ql_last   = omega_ql[t_last_pos]
        gamma_ql_last   = gamma_ql[t_last_pos]
        phi2data_last   = phi2data[t_last_pos]
        
        if view and i < ns:
            ax1.plot(t_f, flux_ql_f, line[i], color=color[i], linewidth=4, label='species ='+str(i))

    if view:
        ax1.autoscale()
        ax1.legend(loc='best',labelspacing=0.0, prop={'size':24})
        show()

    data = list([kxrhoi, kyrhoi]) + list([omega_ql_last, gamma_ql_last, omega_full_last, gamma_full_last, phi2data_last]) +\
           list(flux_ql_last)
    # We return the number of species to know how many columns are expected
    return ns, array(data)
    
def phi(case, yrange=None, xrange=None):
    datafile    = outfile(case, quant='final_fields')
    
    vzed, nzed, iz0  = zed(case)[0], zed(case)[1], zed(case)[2]
    n_ky        = ky(case)[1]
    n_kx        = kx(case)[1]
    data        = loadtxt(datafile, dtype='float').reshape(nzed, n_kx, n_ky, 9)

    if n_kx == 1 and n_ky == 1:
        phi_re  = data[:, 0,0,4]
        phi_im  = data[:, 0,0,5]

        xlims = [vzed.min()*pi, vzed.max()*pi]
        ylims = [min(phi_im.min(), phi_re.min()), max(phi_im.max(), abs(phi_re).max())]
        
        pl2y(vzed*pi, phi_re, phi_im, xlabel='$\\zeta$',\
             ylabel='$\\tilde{\\varphi}(\\zeta)$',\
             title="$(k_x, k_y)=("+str(kx(case)[0][0])+", " + str(ky(case)[0][0])+")$",\
             xrange=xlims, yrange=ylims,\
             key1='$\\Re(\\tilde{\\varphi})(\\zeta)$',\
             key2='$\\Im(\\tilde{\\varphi})(\\zeta)$',\
             ax=None)        
        show()

def phi_last(case, ikx=0, iky=0):
    phi                      = phi_vs_t(case)[0]
    nt, nzed, nkx, nky, nphi = shape(phi)
    vzed, nzed, iz0          = zed(case)[0], zed(case)[1], zed(case)[2]
    
    phi_last = phi[nt-1,:,ikx,iky,:]
    phi_re   = phi_last[:,0]
    phi_im   = phi_last[:,1]
    
    xlims = [vzed.min(), vzed.max()]
    ylims = [min(phi_im.min(), phi_re.min()), max(phi_im.max(), abs(phi_re).max())]
        
    pl2y(vzed, phi_last[:,0], phi_last[:,1], xlabel='$\\zeta$',\
         ylabel='$\\tilde{\\varphi}(\\zeta)$',\
         title="$(k_x, k_y)=("+str(kx(case)[ikx][iky])+", " + str(ky(case)[ikx][iky])+")$",\
         xrange=xlims, yrange=ylims,\
         key1='$\\Re(\\tilde{\\varphi})(\\zeta)$',\
         key2='$\\Im(\\tilde{\\varphi})(\\zeta)$',\
         ax=None)
    show()
    

def vmecgeo(case):
    # Function to read and represent geometric quantities

    zeta  = zed(case)[0]
    quant = geo(case)


    # Labels
    bmag     = '$B/B_{\mathrm{ref}}$'
    gradpar  = '$L_{\\mathrm{ref}}\\nabla_{\|} z$'
    gbdrift  = '$2B_{\\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla y/B^3$'
    gbdrift0 = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla x/B^3$'
    cdrift   = '$2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\mathbf{\kappa}\\cdot\\nabla y/B^2$'
    cdrift0  = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\kappa\\cdot\\nabla x/B^2$'

    # Figures and axes
    
    fig = plt.figure(figsize=(18, 12))
    fig.subplots_adjust(left=0.1, wspace=0.35)
    
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    plxy(zeta, quant[0], xlabel='$\\zeta$', ylabel=bmag,      title='\\texttt{bmag}'   , ax=ax1)
    plxy(zeta, quant[1], xlabel='$\\zeta$', ylabel=gradpar,   title='\\texttt{gradpar}', ax=ax2)
    plxy(zeta, quant[2], xlabel='$\\zeta$', ylabel=gbdrift,   title='\\texttt{gbdrift}', ax=ax3)
    plxy(zeta, quant[3], xlabel='$\\zeta$', ylabel=gbdrift0,  title='\\texttt{gbdrift0}',ax=ax4)
    plxy(zeta, quant[4], xlabel='$\\zeta$', ylabel=cdrift,    title='\\texttt{cdrift}',  ax=ax5)
    plxy(zeta, quant[5], xlabel='$\\zeta$', ylabel=cdrift0,   title='\\texttt{cdrift0}', ax=ax6)    

    fig2 = plt.figure(figsize=(18, 6))
    fig2.subplots_adjust(left=0.1, wspace=0.35)
    
    ax7 = fig2.add_subplot(131)
    ax8 = fig2.add_subplot(132)
    ax9 = fig2.add_subplot(133)

    gs2   = '$|\\nabla y|^2$'
    gds21 = '$\hat{s}^2\\nabla x\\cdot\\nabla y$'
    gds22 = '$\hat{s}^2|\\nabla x|^2$' 
    
    plxy(zeta, quant[6], xlabel='$\\zeta$', ylabel=gs2,      title='\\texttt{gds2}'   , ax=ax7)
    plxy(zeta, quant[7], xlabel='$\\zeta$', ylabel=gds21,    title='\\texttt{gds21}',   ax=ax8)
    plxy(zeta, quant[8], xlabel='$\\zeta$', ylabel=gds22,    title='\\texttt{gds22}',   ax=ax9)    


    show()

    
#    bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22
    

def omega_max(case):
    return omega(case, last=True)[:,:,3].max(), omega(case, last=True)[:,:,4].max()


def ktkn(hpc=None, equil=None, process=False):

    hpc ='marconis'
    equil = 'w7xr003s'
    nLn = 21
    nLT = 21
    nky = 20
    run0=101
    irun=run0
    outfile=datadir() + 'marconis_w7xr003s_0101_0441.dat'

    omega_r     = empty((nLn, nLT, nky), dtype='float')
    omega_i     = empty((nLn, nLT, nky), dtype='float')
    omega_r_max = empty((nLn, nLT),      dtype='float')
    omega_i_max = empty((nLn, nLT),      dtype='float')    
    mLn         = empty((nLn, nLT),      dtype='float')
    mLT         = empty((nLn, nLT),      dtype='float')

    if process == True:
        # The runs are analyzed and the file compiling all the data
        # is written
        f = open(outfile, 'w')
        f.write('# (0) run   (1) i_Ln  (2) i_LT   (3) nprim   (4) tprim   '+\
                '(5) Max(Re(omega)) (6) Max(Im(omega))' + '\n')
        
        print("Writing file :", outfile + '\n')
        
        for i in arange (0, nLn):
            print('# (0) run   (1) i_Ln  (2) i_LT   (3) nprim   (4) tprim   '+\
                  '(5) Max(Re(omega)) (6) Max(Im(omega))')
            
            for j in arange (0, nLT):
                
                if hpc != None:
                    case  = hpc + '/' + equil + '_' + str(format8(irun))
                else:
                    case  = equil + '_' + str(format8(irun)) 
                    
                nprim = read_stella_float(case, 'fprim')[0]
                tprim = read_stella_float(case, 'tprim')[0]

                mLn[i, j]         = nprim
                mLT[i, j]         = tprim
                
                omega_r[i, j, :]  = omega(case, last=True)[0, :, 3]
                omega_i[i, j, :]  = omega(case, last=True)[0, :, 4]
                omega_r_max[i, j] = omega(case, last=True)[0, :, 3].max()
                omega_i_max[i, j] = omega(case, last=True)[0, :, 4].max()
                
                print(case, format8(i), format8(j), format2(nprim), format2(tprim),\
                      format2(omega_r[i, j, :].max()), format2(omega_i[i, j, :].max()))
                
                f.write(case+'\t'+str(format8(i))+'\t'+str(format8(j))+'\t'+\
                        str(format2(nprim))+'\t'+str(format2(tprim))+'\t'+\
                        str(format2(omega_r[i, j, :].max()))+'\t'+\
                        str(format2(omega_i[i, j, :].max()))+'\n')
                          
                irun = irun + 1

        f.close()

    if process == False:
        # Just read the existing file with all the data from the set of runs.
        print("Reading existing file : ", outfile)

        data = loadtxt(outfile, dtype='float', usecols=[3,4,5,6])
        mLn  = data[:,0].reshape((nLn, nLT))
        mLT  = data[:,1].reshape((nLn, nLT))
        omega_r_m = data[:,2].reshape((nLn, nLT))
        omega_i_m = data[:,3].reshape((nLn, nLT))      

        fig, ax = plt.subplots(figsize=(11.0, 8.4))
        col = ax.pcolormesh(mLn[:,0], mLT[0,:], omega_r_m.T, vmin=0, vmax=1.8)
#        CS1 = ax.contour((mLn[:,0], mLT[0,:], omega_r_m.T, levels=[0.1, 0.3, 0.5, 0.7, 0.9, 1.0], colors='white')
#        ax.clabel(CS1, inline=1, fontsize=14,  fmt='%1.1f', linewidths = 1.5)
        ax.set_xlabel('$a/L_{n_i}$')
        ax.set_ylabel('$a/L_{T_i}$')
        ax.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator(10))
        cbar=plt.colorbar(col)
        cbar.set_label('$\\Im({\\omega})^{\\textrm{max}}$')
        cbar.set_label('$a\\Im({\\omega})^{max}/v_{\\mathrm{th,i}}$')
        ax.set_aspect(1.0)
        title = '$k_{x}=0$; $0\\le k_{y}\le 10$'
        ax.set_title(title, y=1.02)
        ax.grid(color='grey', linestyle='-', linewidth=0.3)
                
                #        print mLn[:,0]
                #        print mLT[0,:]
                #        print omega_i_m
                #        mLT  = data[0,:].reshape(nLn, nLT)
                
                #    cmap(xdata=mLn[:,0], ydata=mLT[0,:], zdata=omega_r_m,\
                #         xlabel='$-a/L_{n_i}$', ylabel='$-a/L_{T_i}$', zlabel='$\\Im(\\omega)/v_{t}$')
                
        show()
