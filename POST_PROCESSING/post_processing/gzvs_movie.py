import numpy as np
from stella_plots import movie_2d, plot_2d
from stella_data import gzvs, vpa, zed, ntime, g2nozonal_vs_zvpas, nzed, time
from stella_data import gzvs_zonal, kx_stella
from stella_data import outdir
import matplotlib.pyplot as plt

nvpa = vpa.size
print('nvpa', nvpa)
gmax = np.arange(ntime,dtype=float)
gmin = np.arange(ntime,dtype=float)

kx_val = 0.05
ikx = next((i for i, value in enumerate(kx_stella) if value >=kx_val), None)
#(t, species, vpa, tube, zed, kx, ri)
g_zonal = np.array(gzvs_zonal[:, 0, :, 0, :, ikx, 0] + 1j*gzvs_zonal[:, 0, :, 0, :, ikx, 1], dtype=complex) 

for i in range(ntime):
    gmax[i] = np.abs(g_zonal[i,:,:]).max()
    
gmin[:] = -gmax[:]

ylabel = '$v_{\parallel}$'
xlabel = '$z$'
title = '$\int d\mu \int d^2 \mathbf{R} g$'
movie_file = outdir + 'gzvs.mp4'

#movie_2d(g_zonal,zed,vpa,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='RdBu', abs_flag=False)

time_avg = ntime//2
g_zonal_avg = np.sum(g_zonal[time_avg:,:,:],axis=0)
gmax = np.abs(g_zonal_avg).max()
gmin = - gmax

outfile = outdir + 'timeave.pdf'
plot_2d(g_zonal_avg,zed,vpa,gmin,gmax,'zed','vpa',title,cmp='RdBu', outfile=outfile)

