#need to cmsenv somewhere first

import os
import numpy as np

#lumi = [59.0,97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]
lumi = [59.0,97.7, 129.0]

cmd1 = 'mkdir datacards_modified;'
os.system(cmd1)

for lu_no in range(len(lumi)):
    lu = lumi[lu_no]
    print(lu)
    
    command_a = "python card_modifiers/Card_Mod.py --luminosity " + str(lu)
    os.system(command_a)
    
    command_copy_dc = "cp Combined_Mod.txt datacards_modified/dc_"+ str(lu)+".txt"
    os.system(command_copy_dc)
