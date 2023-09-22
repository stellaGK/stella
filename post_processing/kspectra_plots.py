import numpy as np
from matplotlib import pyplot as plt
from stella_data import time, ntime, outdir, file_prefix, nakx, naky, kx, ky
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi2_vs_kxky as phi2
from stella_time import timeavg

phi2avg = np.arange(nakx*naky,dtype=float).reshape(nakx,naky)
fout=open(outdir+file_prefix+'.kxky_spectra','w+')
for ikx in range(phi2.shape[1]):
    for iky in range(phi2.shape[2]):
        phi2avg[ikx,iky] = timeavg(phi2[:,ikx,iky])
        fout.write(str(kx[ikx])+' '+str(ky[iky])+' '+str(phi2avg[ikx,iky])+'\n')
    fout.write('\n')
fout.close()

phi2ky = np.arange(naky,dtype=float)
f2out=open(outdir+file_prefix+'.ky_spectra','w+')
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

file = outdir+file_prefix+'.phi2_kxky_vs_t.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


phi2_ky = np.arange(naky*ntime,dtype=float).reshape(ntime,naky)
for it in range(ntime):
    for iky in range(naky):
        phi2_ky[it,iky] = np.sum(phi2[it,:,iky])

# distinguish between zonal/nonzonal modes
plt.semilogy(time, phi2_ky[:,0],'--')
for iky in range(1,naky):
    plt.semilogy(time, phi2_ky[:,iky])

plt.xlabel('$t (v_{th}/a)$')
plt.xlim(time[0],time[-1])
plt.ylim(max(1.0e-6,np.min(phi2_ky)),np.max(phi2_ky))
plt.title('$|\Phi_{k_y}|^2(t)$')
file = outdir+file_prefix+'.phi2_ky_vs_t.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


phi2_zonal = np.arange(ntime,dtype=float).reshape(ntime)
phi2_non_zonal = np.arange(ntime,dtype=float).reshape(ntime)
for it in range(ntime):
    phi2_zonal[it] = phi2_ky[it,0]
    phi2_non_zonal[it] = np.sum(phi2_ky[it,1:])
    
# distinguish between zonal/nonzonal modes
plt.semilogy(time, phi2_zonal,'--')
plt.semilogy(time, phi2_non_zonal)

plt.xlabel('$t (v_{th}/a)$')
plt.xlim(time[0],time[-1])
plt.ylim(max(1.0e-6,min(np.min(phi2_zonal),np.min(phi2_non_zonal))),max(np.max(phi2_zonal),np.max(phi2_non_zonal)))
plt.title('$|\Phi|^2(t)$')
file = outdir+file_prefix+'.phi2_vs_t.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


