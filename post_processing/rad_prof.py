

from scipy.io import netcdf
import numpy as np
import numpy.matlib


basedir = '/work/e607/e607/dstonge/stella/rad_test/multi/bare_krook/'
right_file  = basedir + 'right.out.nc'
center_file = basedir + 'center.out.nc'
left_file   = basedir + 'left.out.nc'


right_nc  = netcdf.netcdf_file(right_file,'r')
center_nc = netcdf.netcdf_file(center_file,'r')
left_nc   = netcdf.netcdf_file(left_file,'r')


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

def phi_vs_t_to_x(infile,var,ny,nx):
# t ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  #print('c')
  avt_kxky = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  #print('d')
  arr = np.fft.ifft(avt_kxky,axis=2)
  #print('e')

  return arr

def mom_vs_t_to_x(infile,var,ny,nx):
# t nspec ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  avt_kxky = ny*nx*(avt[:,0,0,:,:,:,0] + 1j*avt[:,0,0,:,:,:,1])
  arr = np.fft.ifft(avt_kxky,axis=2)

  return arr

print('0')
naky  = center_nc.dimensions['ky']
nakxl =   left_nc.dimensions['kx']
nakxc = center_nc.dimensions['kx']
nakxr =  right_nc.dimensions['kx']
ky  = np.copy(center_nc.variables['ky'][:])

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

zed  = np.copy(center_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]

jacobl  = np.copy(  left_nc.variables['jacob'][:])
jacobc  = np.copy(center_nc.variables['jacob'][:])
jacobr  = np.copy( right_nc.variables['jacob'][:])

print('1')

dl_over_bl = np.squeeze(delzed*jacobl/sum(delzed*jacobl))
dl_over_bc = np.squeeze(delzed*jacobc/sum(delzed*jacobc))
dl_over_br = np.squeeze(delzed*jacobr/sum(delzed*jacobr))

dobl = np.transpose(np.matlib.tile(dl_over_bl,(naky,nakxl,1)))
dobc = np.transpose(np.matlib.tile(dl_over_bc,(naky,nakxc,1)))
dobr = np.transpose(np.matlib.tile(dl_over_br,(naky,nakxr,1)))

print('2')

#phil_xky = phi_vs_t_to_x(left_nc  ,'phi_vs_t',naky,nakxl)
phic_xky = phi_vs_t_to_x(center_nc,'phi_vs_t',naky,nakxc)
#phir_xky = phi_vs_t_to_x(right_nc ,'phi_vs_t',naky,nakxr)

print('3')

#densl_xky = mom_vs_t_to_x(left_nc  ,'density',naky,nakxl)
densc_xky = mom_vs_t_to_x(center_nc,'density',naky,nakxc)
#densr_xky = mom_vs_t_to_x(right_nc ,'density',naky,nakxr)


#uparl_xky = mom_vs_t_to_x(left_nc  ,'upar',naky,nakxl)
uparc_xky = mom_vs_t_to_x(center_nc,'upar',naky,nakxc)
#uparr_xky = mom_vs_t_to_x(right_nc ,'upar',naky,nakxr)

#templ_xky = mom_vs_t_to_x(left_nc  ,'temperature',naky,nakxl)
tempc_xky = mom_vs_t_to_x(center_nc,'temperature',naky,nakxc)
#tempr_xky = mom_vs_t_to_x(right_nc ,'temperature',naky,nakxr)

#vxl = 1j*ky*phil_xky
vxc = 1j*ky*phic_xky
#vxr = 1j*ky*phir_xky


print('4')
dens_zf = np.real(np.sum(dobc[:,:,0]*densc_xky[:,:,:,0],1))
upar_zf = np.real(np.sum(dobc[:,:,0]*uparc_xky[:,:,:,0],1))
temp_zf = np.real(np.sum(dobc[:,:,0]*tempc_xky[:,:,:,0],1))

temp_ave = np.mean(temp_zf[500:600,:],0)

cout = open(basedir + 'center.stuff','w')

print('5')
for i in range (0, nakxc - 1):
  cout.write('%d ' % i)
  cout.write('%f ' % dens_zf[nt-1,i])
  cout.write('%f ' % upar_zf[nt-1,i])
  cout.write('%f ' % temp_zf[nt-1,i])
  cout.write('%f ' % temp_ave[i])
  cout.write('\n')
#  print(("%d " % i), end='')
#  print( "%f " % dens_zf[nt-1,i], end='', file=cout )
#  print( "%f " % upar_zf[nt-1,i], end='', file=cout )
#  print( "%f " % temp_zf[nt-1,i], end='', file=cout )
#  print( "%f " % temp_ave[i]    , end='', file=cout )
#  print( "" ,file=cout )
cout.close()
exit()
