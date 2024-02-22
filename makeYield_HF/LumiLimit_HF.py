#need to cmsenv somewhere first
# cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv

import os
import numpy as np

lumi = [97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]

HF_lim = []
W_lim = []
HFW_lim = []
ZTT_lim = []

ZTT_lim_mid = []
ZTT_lim_up = []
ZTT_lim_up_2 = []
ZTT_lim_low = []
ZTT_lim_low_2 = []
ZTT_lim = [ZTT_lim_low_2,ZTT_lim_low,ZTT_lim_mid,ZTT_lim_up,ZTT_lim_up_2]

lumi = [97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]



for lu_no in range(len(lumi)):
    lu = lumi[lu_no]
    print(lu)
    
    command_a = "python card_modifiers/HF_test.py --luminosity " + str(lu);
    os.system(command_a)
        
    Whether_Hybrid=True
    
    
    if(Whether_Hybrid):
            command_run = "combine -M HybridNew HF_Combined_Mod.txt --cl 0.9 -t -1 --rMin 0.0 --rMax 15.0 --expectedFromGrid=0.025  > out1.txt"
            os.system(command_run)
            command_run = "combine -M HybridNew HF_Combined_Mod.txt --cl 0.9 -t -1 --rMin 0.0 --rMax 15.0 --expectedFromGrid=0.16  > out2.txt"
            os.system(command_run)
            command_run = "combine -M HybridNew HF_Combined_Mod.txt --cl 0.9 -t -1 --rMin 0.0 --rMax 15.0 --expectedFromGrid=0.5  > out3.txt"
            os.system(command_run)
            command_run = "combine -M HybridNew HF_Combined_Mod.txt --cl 0.9 -t -1 --rMin 0.0 --rMax 15.0 --expectedFromGrid=0.84  > out4.txt"
            os.system(command_run)
            command_run = "combine -M HybridNew HF_Combined_Mod.txt --cl 0.9 -t -1 --rMin 0.0 --rMax 15.0 --expectedFromGrid=0.975  > out5.txt"
            os.system(command_run)
            
            for i_std in range(5):
                    with open("out"+str(int(i_std+1))+".txt", 'r') as f1:
                            if 'Limit: r <' in line:
                                    linsp = line.split()
                                    ZTT_lim[i_std].append(float(linsp[3]))
                                    
    else:
            command_run = "combine -M AsymptoticLimits HF_Combined_Mod.txt --cl 0.9 -t -1  > out1.txt"
            os.system(command_run)
            with open("out1.txt", 'r') as f1:
                for line in f1:
                    if '2.5%' in line:
                        linsp = line.split()
                        ZTT_lim[0].append(float(linsp[-1]))
                    if '16.0%' in line:
                        linsp = line.split()
                        ZTT_lim[1].append(float(linsp[-1]))
                    if '50.0%' in line:
                        linsp = line.split()
                        ZTT_lim[2].append(float(linsp[-1]))
                    if '84.0%' in line:
                        linsp = line.split()
                        ZTT_lim[3].append(float(linsp[-1]))
                    if '97.5%' in line:
                        linsp = line.split()
                        ZTT_lim[4].append(float(linsp[-1]))
            print(ZTT_lim_mid)
    
ZTT_lim = [ZTT_lim_mid,ZTT_lim_low_2,ZTT_lim_low,ZTT_lim_up,ZTT_lim_up_2]
print(ZTT_lim)
