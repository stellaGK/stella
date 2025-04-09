import numpy as np
from stella_plots import movie_2d
from stella_data import phi2_vs_kxky, kx, ky, ntime
from stella_data import es_heat_by_k
from stella_data import nakx, naky
from stella_data import outdir, file_prefix

phi2max = np.arange(ntime,dtype=float)
phi2min = np.arange(ntime,dtype=float)
for i in range(ntime):
    phi2max[i] = np.absolute(phi2_vs_kxky[i,:,:].max())
phi2min[:] = 0.0
ylabel = '$k_x$'
xlabel = '$k_y$'
title = '$\left|\\varphi(k_x,k_y)\\right|^2$'
movie_file = outdir+file_prefix+'.phi2_vs_kxky.mp4'
movie_2d(phi2_vs_kxky,ky,kx,phi2min,phi2max,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')

# make a movie with the zonal flow removed
phi2_non_zonal = np.arange(ntime*naky*nakx,dtype=float).reshape(ntime,nakx,naky)
phi2_non_zonal[:,:,0] = 0.0
phi2_non_zonal[:,:,1:] = phi2_vs_kxky[:,:,1:]
for i in range(ntime):
    phi2max[i] = np.absolute(phi2_non_zonal[i,:,:].max())
phi2min[:] = 0.0
ylabel = '$k_x$'
xlabel = '$k_y$'
title = '$\left|\\varphi(k_x,k_y)\\right|^2$'
movie_file = outdir+file_prefix+'.phi2_non-zonal_vs_kxky.mp4'
movie_2d(phi2_non_zonal,ky,kx,phi2min,phi2max,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')

qflx_max = np.arange(ntime,dtype=float)
qflx_min = np.arange(ntime,dtype=float)
for i in range(ntime):
    qflx_max[i] = np.absolute(es_heat_by_k[i,0,:,:].max())
qflx_min[:] = 0.0
ylabel = '$k_x$'
xlabel = '$k_y$'
title = '$Q_i$'
movie_file = outdir+file_prefix+'.qflx_vs_kxky.mp4'
movie_2d(es_heat_by_k[:,0,:,:],ky,kx,qflx_min,qflx_max,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')
