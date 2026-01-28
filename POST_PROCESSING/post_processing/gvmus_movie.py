import numpy as np
from stella_plots import movie_2d
from stella_data import gvmus, vpa, mu, ntime, g2nozonal_vs_vpamus

gmax = np.arange(ntime,dtype=float)
gmin = np.arange(ntime,dtype=float)

g_zonal = gvmus - g2nozonal_vs_vpamus 
for i in range(ntime):
    gmax[i] = np.absolute(g_zonal[i,0,:,:].max())
    
gmin[:] = - gmax[:]
xlabel = '$v_{\parallel}$'
ylabel = '$\mu$'
title = '$\int d^3 \mathbf{R} g^2$'
movie_file = 'gvmus.mp4'
movie_2d(g_zonal[:,0,:,:],vpa,mu,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu', abs_flag=True)
