import numpy as np
from matplotlib import pyplot as plt
from stella_data import time, ntime, nakx, naky, kx, ky, nspec
from stella_data import outdir, file_prefix
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi2_vs_kxky as phi2
from stella_data import qflux_vs_kxkys, entropy_vs_kx
from stella_time import timeavg

from stella_data import density, upar, temperature
from stella_data import jacob, delzed


fig, axs = plt.subplots(ncols=3, layout='constrained', sharex='col', figsize=(12,8))

dl_over_b = np.squeeze(delzed*jacob)
dl_over_b[-1] = 0.0
dl_over_b = dl_over_b/sum(dl_over_b)

den_full = density[:,0,0,:,:,0,0] + 1j*density[:,0,0,:,:,0,1]
upar_full = upar[:,0,0,:,:,0,0] + 1j*upar[:,0,0,:,:,0,1]
temp_full = temperature[:,0,0,:,:,0,0] + 1j*temperature[:,0,0,:,:,0,1]

den_avg = np.zeros((ntime, nakx), dtype=complex)
upar_avg = np.zeros((ntime, nakx), dtype=complex)
temp_avg = np.zeros((ntime, nakx), dtype=complex)

den_plot = np.zeros((ntime, nakx), dtype=float)
upar_plot = np.zeros((ntime, nakx), dtype=float)
temp_plot = np.zeros((ntime, nakx), dtype=float)

for it in range (ntime):
    for ikx in range (nakx):
        den_avg[it,ikx] = np.sum(den_full[it,:,ikx]*dl_over_b[:])
        upar_avg[it,ikx] = np.sum(upar_full[it,:,ikx]*dl_over_b[:])
        temp_avg[it,ikx] = np.sum(temp_full[it,:,ikx]*dl_over_b[:])

den_plot = np.abs(den_avg)
upar_plot = np.abs(upar_avg)
temp_plot = np.abs(temp_avg)

for ikx in range(nakx):
    axs[0].semilogy(time, den_plot[:,ikx])
    axs[1].semilogy(time, upar_plot[:,ikx])
    axs[2].semilogy(time, temp_plot[:,ikx])

axs[0].set_xlabel('$t (v_{t}/a)$')
axs[1].set_xlabel('$t (v_{t}/a)$')
axs[2].set_xlabel('$t (v_{t}/a)$')

file = outdir+file_prefix+'.stella_zona_den_upar_temp.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')

## PHI2
phi2avg = np.arange(nakx*naky,dtype=float).reshape(nakx,naky)
fout=open(outdir+file_prefix+'.kxky_spectra.txt','w+')
for ikx in range(phi2.shape[1]):
    for iky in range(phi2.shape[2]):
        phi2avg[ikx,iky] = timeavg(phi2[:,ikx,iky])
        fout.write(str(kx[ikx])+' '+str(ky[iky])+' '+str(phi2avg[ikx,iky])+'\n')
    fout.write('\n')
fout.close()

phi2ky = np.arange(naky,dtype=float)
f2out=open(outdir+file_prefix+'.ky_spectra.txt','w+')
f2out.write('1) ky, 2) |phi|^2 \n')
for iky in range(naky):
    phi2ky[iky] = np.sum(phi2avg[:,iky])
    f2out.write(str(ky[iky])+' '+str(phi2ky[iky])+'\n')
f2out.close()


# For plotting residual levels
fig, axs = plt.subplots(ncols=2, layout='constrained', sharex='col', figsize=(12,8))
for ikx in range(nakx):
    axs[0].semilogy(time,phi2[:,ikx,0])
    print('phi2(t=0), kx=', np.sqrt(phi2[0,ikx,0])/kx[ikx]**2 )
    axs[1].scatter(kx[ikx], np.sqrt(phi2[-1, ikx,0]/phi2[0,ikx,0]))
    print('ratio =', np.sqrt(phi2[-1, ikx,0]/phi2[0,ikx,0]), 'kx=',kx[ikx])
    
axs[1].set_xscale('log')
axs[0].set_xlabel('$t (v_{t}/a)$')
axs[1].set_label('$k_x \rho_i$') 
axs[0].set_xlim([time[0],time[ntime-1]])
axs[1].set_xlim(0,2.1 )
plt.title('$\Phi(k_x,k_y=0)$')

file = outdir+file_prefix+'.stella_kspectra_zonal.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')

# For plotting free energy 
fig=plt.figure(figsize=(12,8))

for ikx in range(1, nakx):
    plt.semilogy(time, entropy_vs_kx[:,ikx])

plt.semilogy(time, np.sum(entropy_vs_kx[:,:], axis=1), '--') 
plt.xlabel('$t (v_{t}/a)$')
plt.ylabel('Free energy')
file = outdir+file_prefix+'.stella_free_energy.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


# For plotting phi2 vs itime for all kx, ky
fig=plt.figure(figsize=(12,8))

for ikx in range(phi2.shape[1]-1):
    plt.semilogy(time,phi2[:,ikx+1,0],'--')
    
for iky in range(phi2.shape[2]-1):
    for ikx in range(phi2.shape[1]):
        plt.semilogy(time,phi2[:,ikx,iky+1])
        
        plt.xlabel('$t (v_{t}/a)$')
        plt.xlim([time[0],time[ntime-1]])
        plt.title('$\Phi^2(k_x,k_y)$')

file = outdir+file_prefix+'.stella_kspectra.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


# For plotting sum |phi|^2 over time 
fig=plt.figure(figsize=(12,8))

plt.semilogy(time,np.sum(phi2[:,:,0],axis=1),'--')
    
for iky in range(phi2.shape[2]-1):
    plt.semilogy(time,np.sum(phi2[:,:,iky+1],axis=1))
        
    plt.xlabel('$t (v_{t}/a)$')
    plt.xlim([time[0],time[ntime-1]])
    plt.title('$\sum_{k_x}\Phi^2(k_y)$')

file = outdir+file_prefix+'.phi2_ky_vs_t.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


