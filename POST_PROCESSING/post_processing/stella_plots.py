# -*- coding: utf-8 -*-

## some coding for commonly used plots

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import rc

from matplotlib.animation import FuncAnimation, FFMpegWriter

import imageio_ffmpeg as iioff

mpl.rcParams['animation.ffmpeg_path'] = iioff.get_ffmpeg_exe()

mpl.rcParams['text.usetex'] = False

font = {
    'family': 'serif',
    'serif' : ['DejaVu Serif', 'Times New Roman', 'serif'],
    'weight': 'bold',
    'size'  : 20,
}
rc('font', **font)
# setup some plot defaults
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rc('font', size=30)
# rcParams.update({'figure.autolayout': True})

def plot_1d(x,y,xlab,title='',ylab=''):

    fig = plt.figure(figsize=(12,8))
    plt.plot(x,y)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    return fig

def logyplot_1d(x,y,xlab,title='',ylab=''):

    fig = plt.figure(figsize=(12,8))
    plt.semilogy(x,y)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    return fig

def plot_2d(z,xin,yin,zmin,zmax,xlab='',ylab='',title='',cmp='RdBu'):

    fig = plt.figure(figsize=(12,8))
    x,y = np.meshgrid(xin,yin)
    plt.imshow(z, cmap=cmp, vmin=zmin, vmax=zmax,
               extent=[x.min(),x.max(),y.min(),y.max()],
               interpolation='nearest', origin='lower', aspect='auto')
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    return fig

def movie_2d(z, xin, yin, zmin, zmax, nframes, outfile,
             xlab='', ylab='', title='', step=1, cmp='RdBu',
             fps=10, interval=50, abs_flag=True):
    """
    Create/save an animation of even/odd parts of z with a time-varying colorbar.

    Arguments largely follow your original function. zmin and zmax should be
    indexable by frame (or be broadcastable scalars).
    """
    
    z_even_full = 0.5*(z + np.flip(z, axis=2))
    z_odd_full = 0.5*(z - np.flip(z, axis=2))
    
    if abs_flag == True :
        print('abs flag is True')
        z_even = np.sign(np.real(z_even_full)) * np.abs(z_even_full)
        z_odd  = np.sign(np.real(z_odd_full)) * np.abs(z_odd_full)
    else:
        z_even = np.real(z_even_full)#0.5*(z + np.flip(z, axis=1)))
        z_odd  = np.real(z_odd_full) #0.5*(z - np.flip(z, axis=1)))

    # Ensure frames list respects step
    frames = list(range(0, nframes, step))

    print ('shape even=', z_even.shape)
    # Meshgrid for extent
    xg, yg = np.meshgrid(xin, yin)
    extent = [xg.min(), xg.max(), yg.min(), yg.max()]

    # Figure / axes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    ax1.set_title('Even')
    ax2.set_title('Odd')
    ax1.set_xlabel(xlab); ax1.set_ylabel(ylab)
    ax2.set_xlabel(xlab); ax2.set_ylabel(ylab)

    # Create initial images once (frame 0)
    i0 = frames[0]
    im1 = ax1.imshow(z_even[i0], origin='lower', aspect='auto',
                     extent=extent, cmap='RdBu', animated=False)
    im2 = ax2.imshow(z_odd[i0], origin='lower', aspect='auto',
                     extent=extent, cmap='RdBu', animated=False)

    # Set initial clim from provided arrays/scalars:
    try:
        vmin0 = zmin[i0]
        vmax0 = zmax[i0]
    except Exception:
        # fallback if scalars
        vmin0 = np.asarray(zmin).item() if np.ndim(zmin) == 0 else zmin[0]
        vmax0 = np.asarray(zmax).item() if np.ndim(zmax) == 0 else zmax[0]

    im1.set_clim(vmin0, vmax0)
    im2.set_clim(vmin0, vmax0)

    # Create a single colorbar for both axes (use im1 as the mappable)
    cbar = fig.colorbar(im1, ax=[ax1, ax2])
    cbar.set_label('value')

    # Optional title for the whole figure
    if title:
        fig.suptitle(title)

    # Update function for FuncAnimation
    def update(frame_index):
        i = frames[frame_index]
        # update image arrays
        im1.set_data(z_even[i])
        im2.set_data(z_odd[i])

        # update color limits (choose shared scale, using zmin/zmax arrays)
        try:
            vmin = zmin[i]
            vmax = zmax[i]
        except Exception:
            # handle scalars or other shapes gracefully
            vmin = np.asarray(zmin).item() if np.ndim(zmin) == 0 else zmin[0]
            vmax = np.asarray(zmax).item() if np.ndim(zmax) == 0 else zmax[0]

        im1.set_clim(vmin, vmax)
        im2.set_clim(vmin, vmax)

        # Tell colorbar to update its normalization/ticks based on the mappable
        cbar.mappable.set_clim(vmin, vmax)
        cbar.update_normal(im1)

        # Return artists that have changed. Note: cbar is not blit-friendly,
        # so we run with blit=False below.
        return im1, im2, cbar

    # Build animation: frames=len(frames), blit=False so colorbar updates properly
    anim = FuncAnimation(fig, update, frames=len(frames),
                         interval=interval, blit=False)

    # Save the animation
    writer = FFMpegWriter(fps=fps, bitrate=1800)
    anim.save(outfile, writer=writer)

    plt.close(fig)  # close the figure to avoid display when running in scripts
    print(f"Saved movie to: {outfile}")

# def movie_2d(z,xin,yin,zmin,zmax,nframes,outfile,xlab='',ylab='',title='',step=1,cmp='RdBu'):

#     from matplotlib import animation

#     z_even = (np.real(0.5*(z+np.flip(z, axis=1)))**2 + np.imag(0.5*(z+np.flip(z, axis=1)))**2)**0.5
#     z_odd  = (np.real(0.5*(z-np.flip(z, axis=1)))**2 + np.imag(0.5*(z-np.flip(z, axis=1)))**2)**0.5
    
#     print('shape_even=', z_even.shape)
#     print('shape_odd=', z_odd.shape)
#     print('shape xin', xin.shape)
#     print('shape yin', yin.shape)
    
#     fig = plt.figure(figsize=(20,8))
#     x,y = np.meshgrid(xin,yin)
    
#     ims = []
#     ax_1 = plt.subplot(1,2,1)
#     ax_2 = plt.subplot(1,2,2)
#     ax_1.set_xlabel(xlab)
#     ax_1.set_ylabel(ylab)
#     ax_1.set_title('Even')
#     ax_2.set_xlabel(xlab)
#     ax_2.set_ylabel(ylab)
#     ax_2.set_title('Odd')
    
#     for i in range(0,nframes,step):
#         im = ax_1.imshow(z_even[i,:,:], cmap=cmp, vmin=zmin[i], vmax=zmax[i],
#                          extent=[x.min(),x.max(),y.min(),y.max()],
#                          interpolation='nearest', origin='lower', aspect='auto')
        
#         im2 = ax_2.imshow(z_odd[i,:,:], cmap=cmp, vmin=zmin[i], vmax=zmax[i],
#                         extent=[x.min(),x.max(),y.min(),y.max()],
#                         interpolation='nearest', origin='lower', aspect='auto')

#         ims.append([im, im2])

#         cbar = plt.colorbar(im, ax=[ax_1, ax_2])
#         # im = plt.imshow(z[i,:,:], cmap=cmp, vmin=zmin[i], vmax=zmax[i],
#         #        extent=[x.min(),x.max(),y.min(),y.max()],
#         #        interpolation='nearest', origin='lower', aspect='auto')
#         # ims.append([im])

#     writer = FFMpegWriter(fps=10, bitrate=1800)
#     ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True)
#     ani.save(outfile, writer=writer)

def movie_1d(x,y,xmin,xmax,ymin,ymax,nframes,outfile,xlab,ylab):

    from matplotlib import animation

    fig = plt.figure(figsize=(12,8))
    ax=plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
    line, = ax.plot([],[],lw=2)

    plt.xlabel(xlab)
    plt.ylabel(ylab)

    def init():
        line.set_data([],[])
        return line,

    def animate(i):
        line.set_data(x,y[i,:])
        return line,

    anim=animation.FuncAnimation(fig, animate, init_func=init, 
                                 frames=nframes, interval=20)

    anim.save(outfile)
