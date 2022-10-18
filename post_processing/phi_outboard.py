from scipy.io import netcdf
from array import array
import numpy as np
import numpy.matlib

basedir = '/marconi_work/FUA36_EMROT/stonge0/new_runs/rotation/smaller/d2Tdr2/g/'
infile = basedir + 'master.out.nc'
infile_nc = netcdf.netcdf_file(infile,'r')

keepZonal = False
upScale_factor = 2

def read_stella_float(infile, var):

  import numpy as np

  try:
    arr = infile.variables[var][:]
    flag = True
  except KeyError:
    print('INFO: '+var+' not found in netcdf file')
    arr =np.arange(1,dtype=float)
    flag = FLAG

  return arr, flag

def phi_vs_t(infile,var,ny,nx):
# t ntube z kx ky ri -> t z kx ky
  avt, present = read_stella_float(infile,var)
  arr = upScale_factor**2*ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  return arr

print('0')
naky = infile_nc.dimensions['ky']
nakx = infile_nc.dimensions['kx']
ny = 2*naky - 1
nakx2 = (nakx + 1)//2

nx_scale = nakx*upScale_factor
ny_scale = ny*upScale_factor

ky  = np.copy(infile_nc.variables['ky'][:])
kx  = np.copy(infile_nc.variables['kx'][:])

Lx = 2*np.pi/kx[1]
Ly = 2*np.pi/ky[1]
dx = Lx/nx_scale
dy = Ly/ny_scale

t  = np.copy(infile_nc.variables['t'][:])
nt = t.size

zed  = np.copy(infile_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)//2) - 1


print('2')

phi_kxky = phi_vs_t(infile_nc,'phi_vs_t',ny,nakx)
if not keepZonal:
    phi_kxky[:,:,:,0] = 0.0

phi_kxky_s = np.zeros((nt,nzed,nx_scale,naky),dtype=numpy.complex128)
phi_kxky_s[:,:,0:nakx2,:] = phi_kxky[:,:,0:nakx2,:]
phi_kxky_s[:,:,(nx_scale - nakx2+1):,:] = phi_kxky[:,:,nakx2:,:]

phi_xky = np.fft.ifft(phi_kxky_s,axis=2)

phi_xky2 = np.zeros((nt,nzed,nx_scale,ny_scale),dtype=numpy.complex128)

phi_xky2[:,:,:,0:naky] = phi_xky

for j in range (1,naky):
  phi_xky2[:,:,:,ny_scale - j] = np.conjugate(phi_xky[:,:,:,j])

phi_xy = np.real(np.fft.ifft(phi_xky2,axis=3))

bin_array = array('f',[0 for i in range((ny_scale+1)*(nx_scale+1))])
bin_array[0] = nx_scale
for i in range (nx_scale):
  bin_array[i+1] = dx*i

for n in range (nt) :
  for j in range (ny_scale):
    bin_array[(j+1)*(nx_scale+1)] = dy*j
    for i in range (nx_scale):
      bin_array[(j+1)*(nx_scale+1)+i+1] = phi_xy[n,omp,i,j]
  cout = open(basedir + 'phi_omp_nzf_' + str(n),'wb')
  bin_array.tofile(cout)
  cout.close()
