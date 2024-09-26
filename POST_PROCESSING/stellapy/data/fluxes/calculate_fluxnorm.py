 

def calculate_fluxnorm(vmec_filename, poloidal_turns):   

    # Make sure poloidal_turns is rounded
    poloidal_turns = round(poloidal_turns, 1)

    # Detect device
    if vmec_filename=='wout_fake_CBC.nc':   device = "CBC"
    if "w7x" in  vmec_filename.lower():     device = "W7X"
    if "asd" in  vmec_filename.lower():     device = "ASDEX"
    if "ncsx" in vmec_filename.lower():     device = "NCSX"
    if "lhd" in  vmec_filename.lower():     device = "LHD"
    if "tj" in   vmec_filename.lower():     device = "TJII"
    
    # Saved data
    if device=="W7X": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [27.533845258457536, 52.917409056292016, 78.51436947923494, 105.74880657798518, 133.09408266298922, 158.6985935450743]
        sum_jacob_x_delzed_x_grho = [34.08906006949946, 62.583763154001, 92.13004726206152, 126.48023057300153, 162.58072449402445, 193.55249662286352]
    if device=="TJII": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [12.237612077152379, 21.95083847761314, 31.468522474720963, 42.99444528532781, 54.624007616596714, 64.14730837439862]
        sum_jacob_x_delzed_x_grho = [14.95255199457543, 27.721795298072283, 40.96346512691376, 54.03700622413496, 66.79193997973647, 79.88847006654052]
    if device=="CBC (VMEC)": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [24.068847821918464, 37.23541078900025, 49.23324978538552, 69.09228202954738, 96.41103797671367, 111.61802451320045]
        sum_jacob_x_delzed_x_grho = [27.75237822030718, 39.44474170337949, 49.847599668242474, 70.5507697314056, 104.12211745440091, 118.24267803774156]
    if device=="CBC (Miller)": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [1.0, 12.2124061436153 , 1.0, 36.6372184690110, 1.0, 61.0620307944062]
        sum_jacob_x_delzed_x_grho = [1.0, 12.2124061436153 , 1.0, 36.6372184690110, 1.0, 61.0620307944062]
    if device=="NCSX": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [20.638687914705017, 36.54374636951319, 52.97458078903025, 74.02497413538403, 94.6800179018151, 110.52263674328321]
        sum_jacob_x_delzed_x_grho = [26.711676220723486, 47.0653496900442, 68.97166108637191, 97.51072424024971, 124.06844035804295, 147.56203224229444]
    if device=="ASDEX": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [12.019818531206687, 19.256422201457678, 26.445239940579, 38.316786812963294, 50.483716395741915, 57.769234759840614]
        sum_jacob_x_delzed_x_grho = [14.050423938397024, 20.84564318372146, 27.66146757140219, 41.24466825155905, 55.7617983408809, 62.53684372933887]
    if device=="LHD": 
        turns=[0.5,1,1.5,2,2.5,3]
        sum_jacob_x_delzed = [24.10349137023479, 42.81321165733264, 61.71711772860515, 85.92476685446027, 109.91373980858012, 128.4407761841303]
        sum_jacob_x_delzed_x_grho = [25.650392202692647, 45.265877194512726, 64.92931721137728, 90.87316226446941, 116.22057346557887, 135.79675304353947]
            
    # Get the factor
    i = turns.index(poloidal_turns)
    print("Got flux_norm for", poloidal_turns, "turns in", device, ": ", sum_jacob_x_delzed[i]/sum_jacob_x_delzed_x_grho[i])
    return sum_jacob_x_delzed[i]/sum_jacob_x_delzed_x_grho[i]
    