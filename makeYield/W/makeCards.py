#need to cmsenv somewhere first

import os
import numpy as np

#lumi = [59.0,97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]
#lumi = [97.7,129.0,150,377.0,500.0,700.0,1000.0, 1200.0, 1500.0, 1700.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0]
lumi = [97.7,129.0,150,377.0,500.0,550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 1000.0, 1050.0, 1100.0, 1200.0, 1250.0, 1300.0, 1350.0, 1400.0, 1450.0,  1500.0, 1550.0, 1600.0, 1650.0,
1700.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0]

#lumi = [59.0,97.7, 129.0]

cmd1 = 'mkdir datacards_modified;'
os.system(cmd1)

for lu_no in range(len(lumi)):
    lu = lumi[lu_no]
    print(lu)
    
    command_a = "python card_modifiers/card_mod.py --luminosity " + str(lu)
    os.system(command_a)
    
    command_copy_dc = "cp Combined_Mod.txt datacards_modified/dc_"+ str(lu)+".txt"
    os.system(command_copy_dc)
