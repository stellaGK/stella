import numpy as np
import stella_data as sd
from matplotlib import pyplot as plt
from stella_data import time, ntime, nakx, naky, kx, ky, nzed
from stella_data import outdir, file_prefix
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi_vs_t as phi_kxky_ri
from stella_plots import plot_2d

ncopies = 5

ny = 2*naky - 1
nx = nakx
nz = nzed

Lx = 2*np.pi/kx[1]
Ly = 2*np.pi/ky[1]
Lz = 2*np.pi

dx = Lx/nx
dy = Ly/ny
dz = Lz/nz

# x is our x-grid
x = np.zeros((nx),dtype=np.float64)
for j in range(nx):
    x[j] = dx*j - Lx/2

# y is our y-grid
y = np.zeros((ny),dtype=np.float64)
for j in range(ny):
    y[j] = dy*j - Ly/2

# z is our z-grid
z = np.zeros((nz), dtype=np.float64)
for j in range(nz):
    z[j] = dz*j - Lz/2
    
# Initiate grid for phi - time, zed, kx, ky, complex phi
# From netcdf phi_kxky_ri has dim: time, ntubes, zed, kz, ky, (real, imaginary) parts of phi
# NB: ?? netcdf quantities are backwards so second dimension of phi_kxky_ri is ntubes and
#     ?? the final dimension is phi itself. The 0 is for real and 1 is for imaginary
phi_kxky = np.zeros((nzed,nakx,naky),dtype=np.complex128)
phi_kxky = nx*ny*(phi_kxky_ri[-1,0,:,:,:,0] + 1j*phi_kxky_ri[-1,0,:,:,:,1])
    
# Take the 1d inverse fft in x to get phi(x,ky,...)
phi_xky =  np.fft.ifft(phi_kxky,axis=1)
# phi_xky_all wil be the same as phi_xky but with the negative ky values included
# using the reality condition
phi_xky_all = np.zeros((nzed,nx,ny),dtype=np.complex128)
# Fill positive ky values the same
phi_xky_all[:,:,0:naky] = phi_xky
# Fill negative ky values with conjugates of phi (need the full ky spectrum)
for j in range(naky-1):
    phi_xky_all[:,:,ny-(j+1)] = np.conjugate(phi_xky[:,:,j+1])

# Take the 1d inverse fft in y to get phi(x,y)
phi_xy = np.real(np.fft.ifft(phi_xky_all,axis=2))

phi_avg = np.sum(phi_xy, axis=2)

# cmap = 'YlGnBu'
# xlab = '$x$'
# ylab = '$z$'
# title = 'avg $|phi|^2$'
# zmax = np.absolute(phi_avg).max()
# zmin = 0.0

# fig = plot_2d(phi_avg[:,:] ,z ,x ,zmin,zmax,ylab,xlab,title,cmap)

# file = 'phi2_vs_xz.pdf'
# pdf = PdfPages(file)
# pdf.savefig
# pdf.close()
# plt.close('all')

fout=open(outdir+file_prefix+'.phi_xvsz','w+')

fout.write('\n')
for ix in range(nx):
    for iz in range(nz):
        fout.write(str(x[ix])+' '+str(z[iz])+' '+str(phi_avg[iz,ix])+'\n' )
    fout.write('\n')
fout.write('\n')

fout.close()
