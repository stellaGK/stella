import numpy as np
from matplotlib import pyplot as plot
from stella_data import outdir, file_prefix
from stella_data import time, ntime, nakx, naky, kx_stella, ky, iz0, nakx_mid, zed, nzed
from stella_data import safety_factor, magnetic_shear, rho, diota_drho, rhostar, iota
from stella_data import phi_vs_t as phi_kxky_ri
from matplotlib.backends.backend_pdf import PdfPages

# !!! THIS DIAGNOSTIC HAS YET TO BE PROPERLY TESTED -- USE WITH CAUTION !!!

ncopies = 5

ny = 2*naky - 1
nx = nakx

Lx = 2*np.pi/kx_stella[1]
Ly = 2*np.pi/ky[1]
dx = Lx/nx
dy = Ly/ny

x = np.zeros((nx),dtype=np.float64)
for j in range(nx):
    x[j] = dx*j - Lx/2

y = np.zeros((ny),dtype=np.float64)
for j in range(ny):
    y[j] = dy*j - Ly/2

copyshift = np.zeros((ncopies),dtype=np.float64)
for i in range(ncopies):
    copyshift[i] = -Ly*(ncopies//2) + i*Ly
    
ymod = np.zeros((nx,ny,nzed,ncopies),dtype=np.float64)
for j in range(nx):
    for n in range(ncopies):
        for k in range(ny):
            for m in range(nzed):
                ymod[j,k,m,n] = y[k] + rho*(-np.pi + 2*m*np.pi/(nzed-1))*(iota/rhostar + x[j]*diota_drho) + copyshift[n]

ymod_max = np.max(ymod)
ymod_min = np.min(ymod)
            
phi_kxky = np.zeros((ntime,nzed,nakx,naky),dtype=np.complex128)
phi_kxky = nx*ny*(phi_kxky_ri[:,0,:,:,:,0] + 1j*phi_kxky_ri[:,0,:,:,:,1])
    
# take the 1d inverse fft in x to get phi(x,ky,...)
phi_xky =  np.fft.ifft(phi_kxky,axis=2)

# phi_xky_all wil be the same as phi_xky but with the negative ky values included
# using the reality condition
phi_xky_all = np.zeros((ntime,nzed,nx,ny),dtype=np.complex128)
phi_xky_all[:,:,:,0:naky] = phi_xky
for j in range(naky-1):
    phi_xky_all[:,:,:,ny-(j+1)] = np.conjugate(phi_xky[:,:,:,j+1])

# take the 1d inverse fft in y to get phi(x,y)
phi_xy = np.real(np.fft.ifft(phi_xky_all,axis=3))

# write to file phi(x,y,z=zmin,tfinal), phi(x,y,z=0,tfinal) and phi(x,y,zmax,tfinal)
fout=open(outdir+file_prefix+'.phi_xy','w+')
fout.write(str(x[0])+' '+str(x[-1])+' '+str(ymod_min)+' '+str(ymod_max)+' '+str(Ly)+'\n')
fout.write('\n')
for ix in range(nx):
    for icopy in range(ncopies):
        for iy in range(ny):
            fout.write(str(x[ix])+' '+str(y[iy])+' '+str(zed[0])+' '+str(phi_xy[-1,0,ix,iy])+' '+str(ymod[ix,iy,0,icopy])+'\n')
    fout.write('\n')
fout.write('\n')

for ix in range(nx):
    for icopy in range(ncopies):
        for iy in range(ny):
            fout.write(str(x[ix])+' '+str(y[iy])+' '+str(zed[iz0-1])+' '+str(phi_xy[-1,iz0-1,ix,iy])+' '+str(ymod[ix,iy,iz0-1,icopy])+'\n')
    fout.write('\n')
fout.write('\n')

for ix in range(nx):
    for icopy in range(ncopies):
        for iy in range(ny):
            fout.write(str(x[ix])+' '+str(y[iy])+' '+str(zed[-1])+' '+str(phi_xy[-1,-1,ix,iy])+' '+str(ymod[ix,iy,-1,icopy])+'\n')
    fout.write('\n')
fout.write('\n')

fout.close()