  
from stellapy.data.moments.read_moments2D import get_moments2D
from stellapy.data.moments.read_moments3D import get_moments3D
from stellapy.data.moments.read_moments4D import get_moments4D 
from stellapy.data.moments.read_moments5D import get_moments5D 
from stellapy.data.moments.read_phaseshifts import get_phaseshifts
from stellapy.data.utils.calculate_attributeWhenReadFirstTime import calculate_attributeWhenReadFirstTime 

#===============================================================================
#                        CREATE THE MOMENTS OBJECT
#===============================================================================
  
class Moments:
    
    # Copy the data from <simulation> that is needed to construct <moments>
    def __init__(self, simulation): 
        self.path = simulation.path  
        self.vec = simulation.vec
        return
 
    # Get the 2D moments 
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_ts(self):           get_moments2D(self);     return self.upar_vs_ts
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_ts(self):           get_moments2D(self);     return self.dens_vs_ts
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_ts(self):           get_moments2D(self);     return self.temp_vs_ts
    @calculate_attributeWhenReadFirstTime 
    def upar2_vs_ts(self):          get_moments2D(self);     return self.upar2_vs_ts
    @calculate_attributeWhenReadFirstTime 
    def dens2_vs_ts(self):          get_moments2D(self);     return self.dens2_vs_ts
    @calculate_attributeWhenReadFirstTime 
    def temp2_vs_ts(self):          get_moments2D(self);     return self.temp2_vs_ts
 
    # Get the 3D moments 
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_tsz(self):          get_moments3D(self);     return self.upar_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_tsz(self):          get_moments3D(self);     return self.dens_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_tsz(self):          get_moments3D(self);     return self.temp_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def upar2_vs_tsz(self):         get_moments3D(self);     return self.upar2_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def dens2_vs_tsz(self):         get_moments3D(self);     return self.dens2_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def temp2_vs_tsz(self):         get_moments3D(self);     return self.temp2_vs_tsz
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_tskx(self):         get_moments3D(self);     return self.upar_vs_tskx
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_tskx(self):         get_moments3D(self);     return self.dens_vs_tskx
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_tskx(self):         get_moments3D(self);     return self.temp_vs_tskx
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_tsky(self):         get_moments3D(self);     return self.upar_vs_tsky
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_tsky(self):         get_moments3D(self);     return self.dens_vs_tsky
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_tsky(self):         get_moments3D(self);     return self.temp_vs_tsky
    
    # Get the 4D moments 
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_tskxky(self):       get_moments4D(self);     return self.upar_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_tskxky(self):       get_moments4D(self);     return self.dens_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_tskxky(self):       get_moments4D(self);     return self.temp_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def upar2_vs_tskxky(self):      get_moments4D(self);     return self.upar2_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def dens2_vs_tskxky(self):      get_moments4D(self);     return self.dens2_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def temp2_vs_tskxky(self):      get_moments4D(self);     return self.temp2_vs_tskxky
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_tskxky_zeta0(self): get_moments4D(self);     return self.upar_vs_tskxky_zeta0
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_tskxky_zeta0(self): get_moments4D(self);     return self.dens_vs_tskxky_zeta0
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_tskxky_zeta0(self): get_moments4D(self);     return self.temp_vs_tskxky_zeta0
    
    # Get the 5D moments 
    @calculate_attributeWhenReadFirstTime 
    def upar_vs_tszkxky(self):      get_moments5D(self);     return self.upar_vs_tszkxky
    @calculate_attributeWhenReadFirstTime 
    def dens_vs_tszkxky(self):      get_moments5D(self);     return self.dens_vs_tszkxky
    @calculate_attributeWhenReadFirstTime 
    def temp_vs_tszkxky(self):      get_moments5D(self);     return self.temp_vs_tszkxky
    
    # Get the phase shifts 
    @calculate_attributeWhenReadFirstTime 
    def tends(self):                               get_phaseshifts(self);     return self.tends
    @calculate_attributeWhenReadFirstTime 
    def tstarts(self):                             get_phaseshifts(self);     return self.tstarts
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_n0_vs_tkx(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_n0_vs_tkx
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_n1_vs_tkx(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_n1_vs_tkx
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_T0_vs_tkx(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_T0_vs_tkx
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_T1_vs_tkx(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_T1_vs_tkx 
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_n0_and_T0_vs_tkx(self):       get_phaseshifts(self);     return self.phase_shifts_n0_and_T0_vs_tkx
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_n1_and_T1_vs_tkx(self):       get_phaseshifts(self);     return self.phase_shifts_n1_and_T1_vs_tkx 
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_n0_vs_tky(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_n0_vs_tky
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_n1_vs_tky(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_n1_vs_tky
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_T0_vs_tky(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_T0_vs_tky
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_phi_and_T1_vs_tky(self):      get_phaseshifts(self);     return self.phase_shifts_phi_and_T1_vs_tky 
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_n0_and_T0_vs_tky(self):       get_phaseshifts(self);     return self.phase_shifts_n0_and_T0_vs_tky
    @calculate_attributeWhenReadFirstTime 
    def phase_shifts_n1_and_T1_vs_tky(self):       get_phaseshifts(self);     return self.phase_shifts_n1_and_T1_vs_tky 

#-------------------------
def load_momentsObject(self): 
    self.moments = Moments(self) 

            