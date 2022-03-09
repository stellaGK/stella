from scipy.io import netcdf
from array import array
import numpy as np
import numpy.matlib

basedir = '/work/e607/e607/dstonge/rho_CBC/200/'
infile = basedir + 'master.out.nc'
infile_nc = netcdf.netcdf_file(infile,'r')

keepZonal = False

def read_stella_float(infile, var):

  import numpy as np

  try:
    #print('a')
    #arr = np.copy(infile.variables[var][:])
    arr = infile.variables[var][:]
    #print('b')
    flag = True
  except KeyError:
    print('INFO: '+var+' not found in netcdf file')
    arr =np.arange(1,dtype=float)
    flag = FLAG

  return arr, flag

def phi_vs_t(infile,var,ny,nx):
# t ntube z kx ky ri
  avt, present = read_stella_float(infile,var)
  arr = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  arr = np.fft.ifft(arr,axis=2)
  return arr

def mom_vs_t(infile,var,ny,nx):
# t nspec ntube z kx ky ri
  avt, present = read_stella_float(infile,var)
  avt_kxky = ny*nx*(avt[:,0,0,:,:,:,0] + 1j*avt[:,0,0,:,:,:,1])
  arr = np.fft.ifft(avt_kxky,axis=2)
  return arr

print('0')
naky = infile_nc.dimensions['ky']
nakx = infile_nc.dimensions['kx']
ny = 2*naky - 1

ky  = np.copy(infile_nc.variables['ky'][:])
kx  = np.copy(infile_nc.variables['kx'][:])

Lx = 2*np.pi/kx[1]
Ly = 2*np.pi/ky[1]
dx = Lx/nakx
dy = Ly/ny


t  = np.copy(infile_nc.variables['t'][:])
nt = t.size

zed  = np.copy(infile_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)//2) - 1


print('2')

phi_xky = phi_vs_t(infile_nc,'phi_vs_t',naky,nakx)
if not keepZonal:
    phi_xky[:,:,:,0] = 0.0
phi_xky2 = np.zeros((nt,nzed,nakx,ny),dtype=numpy.complex128)

phi_xky2[:,:,:,0:naky] = phi_xky


for j in range (1,naky):
  phi_xky2[:,:,:,ny-j] = np.conjugate(phi_xky[:,:,:,j])

phi_xy = np.real(np.fft.ifft(phi_xky2,axis=3))

bin_array = array('f',[0 for i in range((ny+1)*(nakx+1))])
bin_array[0] = nakx
for i in range (nakx):
  bin_array[i+1] = dx*i

for n in range (nt) :
  for j in range (ny):
    bin_array[(j+1)*(nakx+1)] = dy*j
    for i in range (nakx):
      bin_array[(j+1)*(nakx+1)+i+1] = phi_xy[n,omp,i,j]
  cout = open(basedir + 'phi_omp_' + str(n),'wb')
  bin_array.tofile(cout)
  cout.close()
