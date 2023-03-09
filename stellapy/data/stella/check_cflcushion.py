
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def check_cflcushion(delt=0.1, cfl_cushion_upper=0.5, cfl_cushion_lower=0.1, code_dt_max=0.1, nstep=100):
    """ 
    We always want to keep <code_dt> under <cfl_dt>, <cfl_cushion> is the factor that 
    we keep <code_dt> under <cfl_dt>. If <cfl_dt> increases so that it is less than <cfl_cushion>
    bigger than <code_dt>, then we will decrease <code_dt> to make sure it stays smaller.
    To not change <code_dt> at every time step, we use <delt_adjust> to reduce <code_dt> by 
    an extra factor <delt_adjust>.
    
    We also don't want <code_dt> to to be much smaller than <cfl_dt>, so if <code_dt>
    is even smaller than the <cfl_cushion> and <delt_adjust> factor under <cfl_dt>, 
    then increase it again.
    
    Print cfl_dt in stella to get an accurate idea of this quantity!.
    """
    
    # Define some characteristic delta t's as log10()
    vec_cfl_dt_discrete = [-1., -2., -3., -3., -3., -3., -2., -3., -1., -1] 
    vec_code_dt = [delt]; changes_in_delt = []
    print(0.1/0.22)
    print(0.1, 0.1/0.22*0.5)
    
    # Construct a continues vector of time steps
    vec_cfl_dt = []
    for i in range(len(vec_cfl_dt_discrete)-1):
        vec_cfl_dt += list(vec_cfl_dt_discrete[i] + np.array(range(nstep))/nstep*(vec_cfl_dt_discrete[i+1]-vec_cfl_dt_discrete[i]))
    vec_cfl_dt = 10**np.array(vec_cfl_dt) 
    vec_step = range(len(vec_cfl_dt))
    
    # Mimic the CFL decrease condition
    for i, cfl_dt in enumerate(vec_cfl_dt):
        if (vec_code_dt[-1] > cfl_dt*cfl_cushion_upper):
            print(10**((np.log10(cfl_cushion_upper)+np.log10(cfl_cushion_lower))/2))
            vec_code_dt.append(cfl_dt*10**((np.log10(cfl_cushion_upper)+np.log10(cfl_cushion_lower))/2))
            changes_in_delt.append(i)
            print()
            print(f"DECREASE! Because {vec_code_dt[-2]:6.2e} > {cfl_dt*cfl_cushion_upper:6.2e}")
            print(f"    {cfl_dt*cfl_cushion_upper:6.2e} = cfl_dt*cfl_cushion_upper")
            print(f"    {cfl_dt:6.2e} = cfl_dt")
            print(f"    {vec_code_dt[-2]:6.2e} = code_dt") 
            print(f"  ==> code_dt = {vec_code_dt[-1]}")
        elif (vec_code_dt[-1] < np.min([cfl_dt*cfl_cushion_lower, code_dt_max])):
            vec_code_dt.append(np.min([cfl_dt*10**((np.log10(cfl_cushion_upper)+np.log10(cfl_cushion_lower))/2), code_dt_max]))
            changes_in_delt.append(i)
            print()
            print(f"INCREASE! Because {vec_code_dt[-2]:6.2e} < {np.min([cfl_dt*cfl_cushion_lower, code_dt_max]):6.2e}")
            print(f"    {cfl_dt*cfl_cushion_lower:6.2e} = cfl_dt*cfl_cushion/delt_adjust")
            print(f"    {cfl_dt:6.2e} = cfl_dt")
            print(f"    {vec_code_dt[-2]:6.2e} = code_dt") 
            print(f"  ==> code_dt = {vec_code_dt[-1]}")
        else:
            vec_code_dt.append(vec_code_dt[-1])
            
    # Create a figure
    fig = plt.figure(figsize=(18, 9)); fig.set_tight_layout(False)
    grid_specifications = gridspec.GridSpec(1,1)
    grid_specifications.update(top=0.98, left=0.05, right=0.95, bottom=0.06, wspace=0.35, hspace=0.45)
    ax = plt.subplot(grid_specifications[0])
    
    # Plot dt(istep)
    ax.plot(vec_step, vec_cfl_dt, color='black', label='CFL dt')
    ax.plot(vec_step, vec_cfl_dt*cfl_cushion_upper, color='black', alpha=0.5, label='CFL dt*CFL cushion upper')
    ax.plot(vec_step, vec_cfl_dt*cfl_cushion_lower, color='black', alpha=0.2, label='CFL dt*CFL cushion lower')
    ax.plot(vec_step, vec_code_dt[1:], color='maroon', label='code dt')
    
    # Highlight the changes 
    if False:
        for change in changes_in_delt:
            ax.axvline(x=change, color='maroon', alpha=0.5, zorder=1)
    
    # Show figure
    ax.set_yscale('log')
    ax.autoscale()
    ax.legend(labelspacing=0.0, handlelength=1, shadow=True)
    plt.show()
    return
 
if __name__ == '__main__':
    from stellapy.plot.utils.style.load_styleFigures import load_styleFigures
    load_styleFigures()
    check_cflcushion()