from scipy.io import netcdf
import numpy as np

fpref='master'
basedir = '/Users/denisst-onge/stella/build/RH/'
center_file = basedir +fpref + '.out.nc'


center_nc = netcdf.netcdf_file(center_file,'r')
zf_diag  = np.copy(center_nc.variables['zf_diag'][:])
zf_shape = np.shape(zf_diag)

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

#what measurements are in zf_diag, along with the time
header = ['time',
          'prl_str',
          'prl_str_rad',
          'mirror',
          'mirror_rad',
          'mgnx_drft',
          'mgny_drft',
          'mgn_drft_rad',
          'exb_nl',
          'exb_nl_rad',
          'zf_source'
          ]

#print the header
print('#'),
for i in range(len(header)):
    print('['+str(i+1) + ']'),
    print(header[i]),
print('')

#print the data
for i in range(nt):
  print(t[i]),
  for j in range(zf_shape[1]):
    print(zf_diag[i,j]),
  print('')
