import numpy as np
from matplotlib import pyplot as plt
from stella_data import time, ntime, nakx, naky, kx, ky, nspec
from stella_data import outdir, file_prefix
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi2_vs_kxky as phi2
from stella_data import qflux_vs_kxkys, free_energy_vs_kx
from stella_time import timeavg

from stella_data import density, upar, temperature
from stella_data import jacob, delzed

def plot_zonal_upar_den_temp():
    print('Plotting zonal upar, den, temp')
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
    
    file = outdir+file_prefix+'.stella_zonal_den_upar_temp.pdf'
    pdf = PdfPages(file)
    for i in plt.get_fignums():
        pdf.savefig(i)
    pdf.close()
    plt.close('all')
    
    # ######
    # fig, axs = plt.subplots(ncols=3, layout='constrained', sharex='col', figsize=(12,8))
    # den_full_kx_ky = density[:,0,0,:,1:,:,0] + 1j*density[:,0,0,:,1:,:,1]
    # upar_full_kx_ky = upar[:,0,0,:,1:,:,0] + 1j*upar[:,0,0,:,1:,:,1]
    # temp_full_kx_ky = temperature[:,0,0,:,1:,:,0] + 1j*temperature[:,0,0,:,1:,:,1]
    # nztot = den_full_kx_ky.shape[1]

    # dens_zonal = np.sum(den_plot, axis=1)
    # upar_zonal = np.sum(upar_plot, axis=1)
    # temp_zonal = np.sum(temp_plot, axis=1)
    
    # den_nonzonal = np.zeros((ntime, nztot), dtype=float)
    # upar_nonzonal = np.zeros((ntime, nztot), dtype=float)
    # temp_nonzonal = np.zeros((ntime, nztot), dtype=float)
    
    # den_nonzonal = np.sum(np.sum(den_full_kx_ky, axis=3), axis=2)
    # upar_nonzonal = np.sum(np.sum(upar_full_kx_ky, axis=3), axis=2)
    # temp_nonzonal = np.sum(np.sum(temp_full_kx_ky, axis=3), axis=2)

    # den_avg_nonzonal = np.zeros((ntime), dtype=complex)
    # upar_avg_nonzonal = np.zeros((ntime), dtype=complex)
    # temp_avg_nonzonal = np.zeros((ntime), dtype=complex)
    
    # den_avg_nonzonal_plt = np.zeros((ntime), dtype=float)
    # upar_avg_nonzonal_plt = np.zeros((ntime), dtype=float)
    # temp_avg_nonzonal_plt = np.zeros((ntime), dtype=float)
    
    # for it in range (ntime):
    #     den_avg_nonzonal[it] = np.sum(den_nonzonal[it,:]*dl_over_b[:])
    #     upar_avg_nonzonal[it] = np.sum(den_nonzonal[it,:]*dl_over_b[:])
    #     temp_avg_nonzonal[it] = np.sum(den_nonzonal[it,:]*dl_over_b[:])
    
    # den_avg_nonzonal_plt = np.abs(den_avg_nonzonal)
    # upar_avg_nonzonal_plt = np.abs(upar_avg_nonzonal)
    # temp_avg_nonzonal_plt = np.abs(temp_avg_nonzonal)
    
    # axs[0].semilogy(time, dens_zonal, '--')
    # axs[0].semilogy(time, den_avg_nonzonal_plt)
    # axs[1].semilogy(time, upar_zonal, '--')
    # axs[1].semilogy(time, upar_avg_nonzonal_plt)
    # axs[2].semilogy(time, temp_zonal, '--')
    # axs[2].semilogy(time, temp_avg_nonzonal_plt)
    
    # axs[0].set_xlabel('$t (v_{t}/a)$')
    # axs[1].set_xlabel('$t (v_{t}/a)$')
    # axs[2].set_xlabel('$t (v_{t}/a)$')
    
    # file = outdir+file_prefix+'.stella_den_upar_temp.pdf'
    # pdf = PdfPages(file)
    # for i in plt.get_fignums():
    #     pdf.savefig(i)
    # pdf.close()
    # plt.close('all')
    
    # ######

def plot_zonal_den_upar_temp_vs_kx():
    fig, axs = plt.subplots(ncols=3, layout='constrained', sharex='col', figsize=(12,8))
    ntime_avg = 3*ntime//4
    avg_den = np.sum(den_plot[ntime_avg:,:], axis=0)/(ntime - ntime_avg)
    avg_upar = np.sum(upar_plot[ntime_avg:,:], axis=0)/(ntime - ntime_avg)
    avg_temp = np.sum(temp_plot[ntime_avg:,:], axis=0)/(ntime - ntime_avg)
    
    for ikx in range (nakx):
        axs[0].scatter(kx[ikx], avg_den[ikx])
        axs[1].scatter(kx[ikx], avg_upar[ikx])
        axs[2].scatter(kx[ikx], avg_temp[ikx])
    
    axs[0].set_xlabel('$k_x$')
    axs[1].set_xlabel('$k_x$')
    axs[2].set_xlabel('$k_x$')

    file = outdir+file_prefix+'.stella_zonal_den_upar_temp_vs_kx.pdf'
    pdf = PdfPages(file)
    for i in plt.get_fignums():
        pdf.savefig(i)
    pdf.close()
    plt.close('all')


# For plotting residual levels
def plot_rh_vs_kx():
    print('Plotting residual vs kx')
    mean = np.zeros((nakx), dtype=float)
    std = np.zeros((nakx), dtype=float)

    fig, axs = plt.subplots(ncols=2, layout='constrained', sharex='col', figsize=(12,8))
    
    t_start = 6.0E+01
    it_start = next((i for i, value in enumerate(time) if value >=t_start), None)

    ntime_avg = 200 #3*ntime//4  #7*ntime//8
    it_avg= next((i for i, value in enumerate(time) if value >=ntime_avg), None)
    time_end = ntime
    it_end = next((i for i, value in enumerate(time) if value >=time_end), None)
    
    #kx_val = 0.05
    #ikx = next((i for i, value in enumerate(kx[nakx//2 + 1:]) if value <=kx_val), None)
    #ikx2 = nakx//2 + 1 + ikx
    #print('kx=', kx)
    #print('kxvaal', ikx, kx[ikx], kx[ikx2], ikx2)
    for ikx in range(nakx):
        mean[ikx] = np.mean(np.sqrt(np.abs(phi2[it_avg:it_end:, ikx, 0])), axis=0)/np.sqrt(np.abs(phi2[it_start, ikx, 0]))
        #/np.sqrt(np.abs(phi2[it_start, ikx, 0])), axis=0)
        std[ikx] = np.std(np.sqrt(np.abs(phi2[it_avg:it_end, ikx, 0])/np.abs(phi2[it_start, ikx, 0])), axis=0)
        #std[ikx] = np.std(np.sqrt(phi2[ntime//2:, ikx, 0]/phi2[it_start, ikx, 0]), axis=0)
        
        axs[0].semilogy(time,phi2[:,ikx,0])
        #axs[1].scatter(kx[ikx], mean[ikx])
        
        axs[1].errorbar(kx[ikx], mean[ikx], yerr=std[ikx], fmt='o', capsize=5)

    print('mean =', mean) 
    print('kx')
    print(kx)
    print('std')
    print(std)
    #print(kx[ikx2])
    #print('mean')
    #print(mean[ikx2])
    #print('std')
    #print(std[ikx2])

    axs[1].set_xscale('log')
    axs[0].set_xlabel('$t (v_{t}/a)$')
    axs[1].set_label('$k_x \rho_i$') 
    axs[0].set_xlim([time[0],time[ntime-1]])
    plt.title('$\Phi(k_x,k_y=0)$')
    
    file = outdir+file_prefix+'.rh_kspectra_zonal.pdf'
    pdf = PdfPages(file)
    for i in plt.get_fignums():
        pdf.savefig(i)
    pdf.close()
    plt.close('all')
    data = np.column_stack((kx, mean, std))
    
    np.savetxt(outdir+"output.txt", data,
               header="kx mean std",
               fmt="%.6e")  
    
def plot_free_energy():
    print('Plotting free energy vs time')
    # For plotting free energy 
    fig=plt.figure(figsize=(12,8))
    
    for ikx in range(1, nakx):
        plt.semilogy(time, free_energy_vs_kx[:,ikx])
        
        plt.semilogy(time, np.sum(free_energy_vs_kx[:,:], axis=1), '--') 
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

def plot_phi2_ky_vs_t():
    print('Plotting phi2 ky vs t')
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


if __name__ == "__main__":
    #plot_zonal_upar_den_temp
    #plot_zonal_den_upar_temp_vs_kx
    plot_rh_vs_kx()

    #plot_free_energy()
    #plot_phi2_ky_vs_t()
