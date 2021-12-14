

from scipy.io import netcdf
import numpy as np
import numpy.matlib

basedir = '/Users/denisst-onge/stella/build/dirichlet/'
infile = basedir + 'master.out.nc'


infile_nc = netcdf.netcdf_file(infile,'r')


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
dy = Ly/naky


t  = np.copy(infile_nc.variables['t'][:])
nt = t.size

zed  = np.copy(infile_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]

jacob  = np.copy(infile_nc.variables['jacob'][:])
dl_over_b = np.squeeze(delzed*jacob)
dl_over_b[nzed-1] = 0.0
dl_over_b = dl_over_b/sum(dl_over_b)

dob = np.transpose(np.matlib.tile(dl_over_b,(naky,nakx,1)))

print('2')

phi_xky = phi_vs_t(infile_nc,'phi_vs_t',naky,nakx)
phi_xky2 = np.zeros((nt,nzed,nakx,ny),dtype=numpy.complex128)

phi_xky2[:,:,:,0:naky] = phi_xky


for j in range (1,naky):
  phi_xky2[:,:,:,ny-j] = np.conjugate(phi_xky[:,:,:,j])

phi_xy = np.real(np.fft.ifft(phi_xky2,axis=3))

for n in range (nt) :
  cout = open(basedir + 'phi_omp_' + str(n),'w')
  for i in range (nakx):
    for j in range (ny):
      cout.write('%f ' % (dx*i))
      cout.write('%f ' % (dy*j))
      cout.write('%f ' % phi_xy[n,omp,i,j])
      cout.write('\n')
    cout.write('\n')
  cout.close()
