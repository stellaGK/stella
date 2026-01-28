import numpy as np
from stella_plots import movie_2d
from stella_data import gzvs, vpa, zed, ntime, g2nozonal_vs_zvpas, nzed
from stella_data import gzvs_zonal
from stella_data import outdir

nvpa = vpa.size
print('nvpa', nvpa)
gmax = np.arange(ntime,dtype=float)
gmin = np.arange(ntime,dtype=float)

print('SHAPE=', gzvs_zonal.shape)

#g_zonal = np.zeros((ntime, nvpa,  nzed), dtype=complex)
g_zonal = np.array(gzvs_zonal[:, 0, :, :, 0] + 1j*gzvs_zonal[:, 0, :, :, 1], dtype=complex) #gzvs - g2nozonal_vs_zvpas

print('shape zonal', np.shape(g_zonal))

for i in range(ntime):
    gmax[i] = np.abs(g_zonal[i,:,:]).max()

gmin[:] = - gmax[:]
ylabel = '$v_{\parallel}$'
xlabel = '$z$'
title = '$\int d\mu \int d^2 \mathbf{R} g^2$'
movie_file = outdir + 'gzvs.mp4'

movie_2d(g_zonal,zed,vpa,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='RdBu', abs_flag=True)
