#need to cmsenv somewhere first
# cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv

import os
import numpy as np

HF_lim = []
W_lim = []
HFW_lim = []
ZTT_lim = []

ZTT_lim_mid = []
ZTT_lim_up = []
ZTT_lim_up_2 = []
ZTT_lim_low = []
ZTT_lim_low_2 = []

lumi = [97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]

card_modifier_name = ["ZTT_tauh_test","ZTT_taumu_test","ZTT_taue_test"]
combined_card_name = ["ZTT_tauh_Combined_Mod","ZTT_taumu_Combined_Mod","ZTT_taue_Combined_Mod"]

# opening the file in read mode; output from makeYield
file_fam = open("ZTT_Lumi_Limit.txt", "r")
# using the for loop
lum_tau = []
lum_mu = []
lum_e = []
for linenum, line in enumerate(file_fam):
        line_mod=line
        x1 = line.split()
        if x1[5]=='tau_h':
            lum_tau.append([x1[1],x1[9],x1[11]])
        if x1[5]=='tau_mu':
            lum_mu.append([x1[1],x1[9],x1[11]])
        if x1[5]=='tau_e':
            lum_e.append([x1[1],x1[9],x1[11]])
file_fam.close()

lum_comb = [lum_tau,lum_mu,lum_e]

for lu_no in range(len(lum_tau)):
    lu = lum_comb[0][lu_no][0]
    print(lu)
    for k in range(3):
        
        
        command_a = "python card_modifiers/"+ card_modifier_name[k] +".py --luminosity " + str(59.0) + " --s " + lum_comb[k][lu_no][1] + " --b " + lum_comb[k][lu_no][2] + " --cuttype a";
        #command_b = "python card_modifiers/"+ card_modifier_name[k] +".py --luminosity " + str(59.0) + " --s " + lum_comb[k][lu_no][3] + " --b " + lum_comb[k][lu_no][4] + " --cuttype b";
        os.system(command_a)
        #os.system(command_b)
                                        
        #command_a_and_b = "combineCards.py "+combined_card_name[k]+"_a.txt "+combined_card_name[k]+"_b.txt > "+combined_card_name[k]+".txt"
        #os.system(command_a_and_b)
        #command_run = "combine -M AsymptoticLimits "+combined_card_name[k]+".txt --cl 0.9 -t -1  > out"+ str(k+1) +".txt";
        #print(command_run)
        
    command_h_mu_e = "combineCards.py "+combined_card_name[0]+"_a.txt "+combined_card_name[1]+"_a.txt "+combined_card_name[2]+"_a.txt > ZTT_Combined.txt"
    #command_h_mu_e = "cp " + combined_card_name[0]+".txt ZTT_Combined.txt"
    
    os.system(command_h_mu_e)
    command_run = "combine -M AsymptoticLimits ZTT_Combined.txt --cl 0.9 -t -1  > out1.txt"
    os.system(command_run)
    with open("out1.txt", 'r') as f1:
        for line in f1:
            if '2.5%' in line:
                linsp = line.split()
                ZTT_lim_low_2.append(float(linsp[-1]))
            if '16.0%' in line:
                linsp = line.split()
                ZTT_lim_low.append(float(linsp[-1]))
            if '50.0%' in line:
                linsp = line.split()
                ZTT_lim_mid.append(float(linsp[-1]))
            if '84.0%' in line:
                linsp = line.split()
                ZTT_lim_up.append(float(linsp[-1]))
            if '97.5%' in line:
                linsp = line.split()
                ZTT_lim_up_2.append(float(linsp[-1]))
    print(ZTT_lim_mid)
    
ZTT_lim = [ZTT_lim_mid,ZTT_lim_low_2,ZTT_lim_low,ZTT_lim_up,ZTT_lim_up_2]
print(ZTT_lim)
