

from scipy.io import netcdf
import numpy as np
import numpy.matlib

fpref='center'

basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/rho_scan2/r0.001ky/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001hr/'
center_file = basedir +fpref + '.out.nc'


center_nc = netcdf.netcdf_file(center_file,'r')


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
  #print('c')
  arr = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  #print('d')
  #arr = np.fft.ifft(avt_kxky,axis=2)
  #print('e')

  return arr

naky  = center_nc.dimensions['ky']
nakxc = center_nc.dimensions['kx']

ky  = np.copy(center_nc.variables['ky'][:])
kxc  = np.copy(center_nc.variables['kx'][:])

dx = 2*np.pi/kxc[1]/nakxc

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

zed  = np.copy(center_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]
jacobc  = np.copy(center_nc.variables['jacob'][:])

dl_over_bc = np.squeeze(delzed*jacobc/sum(delzed*jacobc))

dobc = np.transpose(np.matlib.tile(dl_over_bc,(naky,nakxc,1)))

phic_kxky = phi_vs_t(center_nc,'phi_vs_t',naky,nakxc)
phiZF_kxky= np.sum(dobc[:,:,0]*phic_kxky[:,:,:,0],1)

#print(phic_kxky[nt-1,omp,:,0])
#print(phiZF_kxky[nt-1,:])

phizfc_k= np.sum((1j*kxc)**0*dobc[:,:,0]*phic_kxky[:,:,:,0],1)
vxc_k   = np.sum((1j*kxc)**1*dobc[:,:,0]*phic_kxky[:,:,:,0],1)
exbc_k  = np.sum((1j*kxc)**2*dobc[:,:,0]*phic_kxky[:,:,:,0],1)


phizfc = np.real(np.fft.ifft(phizfc_k,axis=1))/naky
vxc    = np.real(np.fft.ifft(vxc_k,axis=1))/naky
exbc   = np.real(np.fft.ifft(exbc_k,axis=1))/naky
#vxc = np.fft.ifft(vxc_k,axis=2)

tave = int(0.7*float(nt))
p_ave  = np.mean(phizfc[tave:nt,:],0)
v_ave  = np.mean(vxc[tave:nt,:],0)
e_ave  = np.mean(exbc[tave:nt,:],0)

window_width=7
w2 = (window_width+1)/2
cumsum_vec = numpy.cumsum(numpy.insert(e_ave, 0, 0))
ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width



cout = open(basedir + fpref + '.exb_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] phi   ')
cout.write('[4] v     ')
cout.write('[5] exb   ')
cout.write('[6] phik  ')
cout.write('\n')
for i in range (nt):
  for j in range (nakxc):
    cout.write('%e ' % t[i])
    cout.write('%e ' % (dx*j))
    cout.write('%e ' % phizfc[i,j])
    cout.write('%e ' % vxc[i,j])
    cout.write('%e ' % exbc[i,j])
    cout.write('%e ' % abs(phizfc_k[i,j]))
    cout.write('\n')
  cout.write('\n')

cout.close()

cout = open(basedir + fpref + '.exb_ave','w')
cout.write('[1] x     ')
cout.write('[2] phi   ')
cout.write('[3] v     ')
cout.write('[4] exb   ')
cout.write('\n')
for j in range (w2,nakxc-w2):
  cout.write('%e ' % (dx*j))
  cout.write('%e ' % p_ave[j])
  cout.write('%e ' % v_ave[j])
  cout.write('%e ' % e_ave[j])
  cout.write('%e ' % ma_vec[j-w2])
  cout.write('\n')

cout.close()
