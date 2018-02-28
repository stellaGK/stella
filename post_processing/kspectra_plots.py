from matplotlib import pyplot as plt
from stella_data import time, ntime, outdir
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi2_vs_kxky as phi2

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

file = outdir+'stella_kspectra.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')
