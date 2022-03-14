import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import ntime, time
from stella_data import xgrid
from stella_data import naky, nakx, nakx_mid
from stella_data import ky, kx_stella, kx
from stella_data import jacob, delzed
from stella_data import phi_vs_t
from stella_data import outdir, file_prefix
from stella_time import timeavg
from stella_plots import movie_1d, plot_1d

# get the integration weights for field-line average
dl_over_b = np.squeeze(delzed*jacob)
dl_over_b[-1] = 0.0
dl_over_b = dl_over_b/sum(dl_over_b)

# phi_vs_t = phi_vs_t[t, ntube, z, kx, ky, ri]
# choose first flux tube, ky=0 and complexify
phi_ky0 = phi_vs_t[:,0,:,:,0,0] + 1j*phi_vs_t[:,0,:,:,0,1]

# get the zonal component of phi by averaging over z
zonal_phi_kx = np.arange(ntime*nakx,dtype=complex).reshape(ntime,nakx)
for it in range (ntime):
  for ikx in range (nakx):
    zonal_phi_kx[it,ikx] = np.sum(dl_over_b*phi_ky0[it,:,ikx])

# compute the zonal flow and zonal flow shear from zonal_phi in k-space
zonal_flow_kx = np.arange(ntime*nakx,dtype=complex).reshape(ntime,nakx)
zonal_shear_kx = np.arange(ntime*nakx,dtype=complex).reshape(ntime,nakx)
for it in range (ntime):
  zonal_flow_kx[it,:] = 1j*kx_stella*zonal_phi_kx[it,:]
  zonal_shear_kx[it,:] = -kx_stella**2*zonal_phi_kx[it,:]

filter_fraction = 0.7
# filter high |kx| from zonal shear
zonal_flow_kx_filtered = np.copy(zonal_flow_kx)
zonal_shear_kx_filtered = np.copy(zonal_shear_kx)
# filter most negative kx values
zonal_flow_kx_filtered[:,nakx_mid:nakx-int((nakx-nakx_mid)*(1.0-filter_fraction))] = 0.0
zonal_shear_kx_filtered[:,nakx_mid:nakx-int((nakx-nakx_mid)*(1.0-filter_fraction))] = 0.0
# filter most positive kx values
zonal_flow_kx_filtered[:,nakx_mid-int(nakx_mid*filter_fraction):nakx_mid] = 0.0
zonal_shear_kx_filtered[:,nakx_mid-int(nakx_mid*filter_fraction):nakx_mid] = 0.0

# compute the spectral energy in the zonal phi
tmp = np.abs(zonal_phi_kx)**2
# re-order kx indices to run from most negative to most positive
zonal_phi_spectrum = np.concatenate((tmp[:,nakx_mid:],tmp[:,:nakx_mid]),axis=1)
# compute the spectral energy in the zonal flow
tmp = np.abs(zonal_flow_kx)**2
# re-order kx indices to run from most negative to most positive
zonal_flow_spectrum = np.concatenate((tmp[:,nakx_mid:],tmp[:,:nakx_mid]),axis=1)
# compute the spectral energy in the zonal shear
tmp = np.abs(zonal_shear_kx)**2
# re-order kx indices to run from most negative to most positive
zonal_shear_spectrum = np.concatenate((tmp[:,nakx_mid:],tmp[:,:nakx_mid]),axis=1)
# compute the spectral energy in the filtered zonal flow
tmp = np.abs(zonal_flow_kx_filtered)**2
# re-order kx indices to run from most negative to most positive
zonal_flow_spectrum_filtered = np.concatenate((tmp[:,nakx_mid:],tmp[:,:nakx_mid]),axis=1)
# compute the spectral energy in the filtered zonal shear
tmp = np.abs(zonal_shear_kx_filtered)**2
# re-order kx indices to run from most negative to most positive
zonal_shear_spectrum_filtered = np.concatenate((tmp[:,nakx_mid:],tmp[:,:nakx_mid]),axis=1)


# Fourier transform the zonal phi, flow and shear from kx to x space
zonal_phi_x = np.real(np.fft.ifft(zonal_phi_kx,axis=1))
zonal_flow_x = np.real(np.fft.ifft(zonal_flow_kx,axis=1))
zonal_shear_x = np.real(np.fft.ifft(zonal_shear_kx,axis=1))
zonal_flow_x_filtered = np.real(np.fft.ifft(zonal_flow_kx_filtered,axis=1))
zonal_shear_x_filtered = np.real(np.fft.ifft(zonal_shear_kx_filtered,axis=1))

# get the time average of the zonal components
zonal_phi_x_tavg = np.arange(nakx,dtype=float)
zonal_flow_x_tavg = np.arange(nakx,dtype=float)
zonal_shear_x_tavg = np.arange(nakx,dtype=float)
zonal_flow_x_filtered_tavg = np.arange(nakx,dtype=float)
zonal_shear_x_filtered_tavg = np.arange(nakx,dtype=float)
zonal_phi_spectrum_tavg = np.arange(nakx,dtype=float)
zonal_flow_spectrum_tavg = np.arange(nakx,dtype=float)
zonal_shear_spectrum_tavg = np.arange(nakx,dtype=float)
zonal_flow_spectrum_filtered_tavg = np.arange(nakx,dtype=float)
zonal_shear_spectrum_filtered_tavg = np.arange(nakx,dtype=float)
for ikx in range(nakx):
  zonal_phi_x_tavg[ikx] = timeavg(zonal_phi_x[:,ikx])
  zonal_flow_x_tavg[ikx] = timeavg(zonal_flow_x[:,ikx])
  zonal_shear_x_tavg[ikx] = timeavg(zonal_shear_x[:,ikx])
  zonal_flow_x_filtered_tavg[ikx] = timeavg(zonal_flow_x_filtered[:,ikx])
  zonal_shear_x_filtered_tavg[ikx] = timeavg(zonal_shear_x_filtered[:,ikx])
  zonal_phi_spectrum_tavg[ikx] = timeavg(zonal_phi_spectrum[:,ikx])
  zonal_flow_spectrum_tavg[ikx] = timeavg(zonal_flow_spectrum[:,ikx])
  zonal_shear_spectrum_tavg[ikx] = timeavg(zonal_shear_spectrum[:,ikx])
  zonal_flow_spectrum_filtered_tavg[ikx] = timeavg(zonal_flow_spectrum_filtered[:,ikx])
  zonal_shear_spectrum_filtered_tavg[ikx] = timeavg(zonal_shear_spectrum_filtered[:,ikx])

print('creating 1D plots of zonal quantities...')
# make 1D plots of the time-averaged zonal phi, flow and shear
xlabel = '$x/\\rho_{\\textnormal{ref}}$'
plot_1d(xgrid, zonal_phi_x_tavg, xlabel, ylab='zonal phi')
plot_1d(xgrid, zonal_flow_x_tavg, xlabel, ylab='zonal flow')
plot_1d(xgrid, zonal_shear_x_tavg, xlabel, ylab='zonal shear')
plot_1d(xgrid, zonal_flow_x_filtered_tavg, xlabel, ylab='filtered zonal flow')
plot_1d(xgrid, zonal_shear_x_filtered_tavg, xlabel, ylab='filtered zonal shear')
xlabel = '$k_x \\rho_{\\textnormal{ref}}$'
plot_1d(kx, zonal_phi_spectrum_tavg, xlabel, ylab='zonal phi spectrum')
plot_1d(kx, zonal_flow_spectrum_tavg, xlabel, ylab='zonal flow spectrum')
plot_1d(kx, zonal_shear_spectrum_tavg, xlabel, ylab='zonal shear spectrum')
plot_1d(kx, zonal_flow_spectrum_filtered_tavg, xlabel, ylab='filtered zonal flow spectrum')
plot_1d(kx, zonal_shear_spectrum_filtered_tavg, xlabel, ylab='filtered zonal shear spectrum')

outfile = outdir+file_prefix+'.zonal_vs_x_tavg.pdf'
pdf = PdfPages(outfile)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')

print('creating 1D movie of the zonal flow...')
xmin = xgrid[0]
xmax = xgrid[-1]
ymin = zonal_flow_x.min()
ymax = zonal_flow_x.max()
xlabel = '$x$'
ylabel = '$zonal flow$'
movie_1d(xgrid, zonal_flow_x, xmin, xmax, ymin, ymax, ntime-1, outdir+file_prefix+'.zonal_vs_x.mp4',xlabel,ylabel)

#p_ave   = np.mean(phizfc[tind:nti,:],0)
#v_ave   = np.mean(vxc[tind:nti,:],0)
#e_ave   = np.squeeze(np.mean(exbc[tind:nti,:],0))
#spec_ave=np.mean(np.abs(exbc_k[tind:nti]),0)

cout = open(outdir + file_prefix + '.zonal_vs_x_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] phi   ')
cout.write('[4] flow  ')
cout.write('[5] shear ')
cout.write('\n')
for it in range (ntime):
  for ikx in range (nakx):
    cout.write('%e ' % time[it])
    cout.write('%e ' % xgrid[ikx])
    cout.write('%e ' % zonal_phi_x[it,ikx])
    cout.write('%e ' % zonal_flow_x[it,ikx])
    cout.write('%e ' % zonal_shear_x[it,ikx])
    cout.write('\n')
  cout.write('\n')

cout.close()

# cout = open(basedir + fpref + '.exb_ave','w')
# cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nti-1]) + '\n')
# cout.write('[1] x     ')
# cout.write('[2] phi   ')
# cout.write('[3] v     ')
# cout.write('[4] exb   ')
# cout.write('[5] exb_k ')
# #cout.write('[5] exb smooth')
# cout.write('[6] spec ave')
# cout.write('\n')
# #for j in range (w2,nakxc-w2):
# for j in range (nakxc):
# #  cout.write('%e ' % radgrid[j,0])
#   cout.write('%e ' % (dx*j))
#   cout.write('%e ' % p_ave[j])
#   cout.write('%e ' % v_ave[j])
#   cout.write('%e ' % e_ave[j])
# #  cout.write('%e ' % ma_vec[j])
#   cout.write('%e ' % spec_ave[j])
#   cout.write('\n')
# cout.close()

# cout = open(basedir + fpref + '.exb_en','w')
# cout.write('[1] t     ')
# cout.write('[2] phi2  ')
# cout.write('[3] uzf2  ')
# cout.write('[4] exb2  ')
# cout.write('[5] exbk2  ')
# cout.write('\n')
# for i in range (nt):
#     cout.write('%e ' % t[i])
#     cout.write('%e ' % np.mean(phizfc[i,:]**2))
#     cout.write('%e ' % np.mean(vxc[i,:]**2))
#     cout.write('%e ' % np.mean(exbc[i,:]**2))
#     cout.write('%e ' % np.mean(np.abs(exbc_k[i,:])**2))
#     cout.write('\n')
# cout.close()

