import numpy as np
from matplotlib import pyplot as plt
from stella_data import time, ntime, nakx, naky, kx, ky, nspec
from stella_data import outdir, file_prefix
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi2_vs_kxky as phi2
from stella_data import es_heat_by_k
from stella_time import timeavg

phi2avg = np.arange(nakx*naky,dtype=float).reshape(nakx,naky)
fout=open(outdir+file_prefix+'.kxky_spectra.txt','w+')
for ikx in range(phi2.shape[1]):
    for iky in range(phi2.shape[2]):
        phi2avg[ikx,iky] = timeavg(phi2[:,ikx,iky])
        fout.write(str(kx[ikx])+' '+str(ky[iky])+' '+str(phi2avg[ikx,iky])+'\n')
    fout.write('\n')
fout.close()

heat_flux_tavg = np.arange(nakx*naky*nspec,dtype=float).reshape(nspec,nakx,naky)
fout=open(outdir+file_prefix+'.heat_flux_kxky_spectra.txt','w+')
for ispec in range(es_heat_by_k.shape[1]):
    for ikx in range(es_heat_by_k.shape[2]):
        for iky in range(es_heat_by_k.shape[3]):
            heat_flux_tavg[ispec,ikx,iky] = timeavg(es_heat_by_k[:,ispec,ikx,iky])
            fout.write(str(ispec)+' '+str(kx[ikx])+' '+str(ky[iky])+' '+str(heat_flux_tavg[ispec,ikx,iky])+'\n')
        fout.write('\n')
    fout.write('\n')
fout.close()

heat_flux_tkxavg = np.arange(naky*nspec,dtype=float).reshape(nspec,naky)
fout=open(outdir+file_prefix+'.heat_flux_ky_spectra.txt','w+')
for ispec in range(es_heat_by_k.shape[1]):
    for iky in range(es_heat_by_k.shape[3]):
        heat_flux_tkxavg[ispec,iky] = np.sum(heat_flux_tavg[ispec,:,iky])
        fout.write(str(ispec)+' '+str(ky[iky])+' '+str(heat_flux_tkxavg[ispec,iky])+'\n')
    fout.write('\n')
fout.close()


phi2ky = np.arange(naky,dtype=float)
f2out=open(outdir+file_prefix+'.ky_spectra.txt','w+')
f2out.write('1) ky, 2) |phi|^2 \n')
for iky in range(naky):
    phi2ky[iky] = np.sum(phi2avg[:,iky])
    f2out.write(str(ky[iky])+' '+str(phi2ky[iky])+'\n')
ky0 = np.sum(phi2ky[1:]*ky[1:])/np.sum(phi2ky[1:])
ky1 = np.sum(phi2ky[1:])/np.sum(phi2ky[1:]/ky[1:])
f2out.write('#ky0: '+str(ky0)+' ky1: '+str(ky1)+'\n')
f2out.close()

fig=plt.figure(figsize=(12,8))

# distinguish between zonal/nonzonal modes
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


fig=plt.figure(figsize=(12,8))

# distinguish between zonal/nonzonal modes
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


fig=plt.figure(figsize=(12,8))
for iky in range(es_heat_by_k.shape[3]-1):
    plt.plot(time,np.sum(es_heat_by_k[:,0,:,iky+1],axis=1))
        
    plt.xlabel('$t (v_{t}/a)$')
    plt.xlim([time[0],time[ntime-1]])
    plt.title('$\sum_{k_x}Q(k_x,k_y)$')

file = outdir+file_prefix+'.qflx_ky_vs_t.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


