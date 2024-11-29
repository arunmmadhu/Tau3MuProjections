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

ZTT_lim = [ZTT_lim_low_2,ZTT_lim_low,ZTT_lim_mid,ZTT_lim_up,ZTT_lim_up_2]

lumi = [59.0, 97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]

card_modifier_name = ["ZTT_tauha_test","ZTT_tauhb_test","ZTT_taumu_test","ZTT_taue_test"]
combined_card_name = ["Cat_1a_Mod","Cat_1b_Mod","Cat_2_Mod","Cat_3_Mod"]

# opening the file in read mode; output from makeYield
file_fam = open("ZTT_Lumi_Limit_2018_final.txt", "r")
# using the for loop
lum_taua = []
lum_taub = []
lum_mu = []
lum_e = []
for linenum, line in enumerate(file_fam):
        line_mod=line
        x1 = line.split()
        if x1[5]=='Cat_h_A':
            lum_taua.append([x1[1],x1[9],x1[11]])
        if x1[5]=='Cat_h_B':
            lum_taub.append([x1[1],x1[9],x1[11]])
        if x1[5]=='Cat_mu':
            lum_mu.append([x1[1],x1[9],x1[11]])
        if x1[5]=='Cat_e':
            lum_e.append([x1[1],x1[9],x1[11]])
file_fam.close()

lum_comb = [lum_taua,lum_taub,lum_mu,lum_e]

for lu_no in range(len(lum_taua)):
#for lu_no in [2]: #for extreme luminosities
    lu = lum_comb[0][lu_no][0]
    print(lu)
    
    Whether_Hybrid=True
    for k in range(4):
        
        
        command_a = "python card_modifiers/"+ card_modifier_name[k] +".py --luminosity " + str(59.0) + " --s " + lum_comb[k][lu_no][1] + " --b " + lum_comb[k][lu_no][2] + " --cuttype a";
        #command_b = "python card_modifiers/"+ card_modifier_name[k] +".py --luminosity " + str(59.0) + " --s " + lum_comb[k][lu_no][3] + " --b " + lum_comb[k][lu_no][4] + " --cuttype b";
        os.system(command_a)
        #os.system(command_b)
                                        
        #command_a_and_b = "combineCards.py "+combined_card_name[k]+"_a.txt "+combined_card_name[k]+"_b.txt > "+combined_card_name[k]+".txt"
        #os.system(command_a_and_b)
        #command_run = "combine -M AsymptoticLimits "+combined_card_name[k]+".txt --cl 0.9 -t -1  > out"+ str(k+1) +".txt";
        #print(command_run)
        
        #Whether_Hybrid = Whether_Hybrid and (float(lum_comb[k][lu_no][1])<35.0 and float(lum_comb[k][lu_no][2])<1.0)
        
    command_h_mu_e = "combineCards.py "+combined_card_name[0]+"_a.txt "+combined_card_name[1]+"_a.txt "+combined_card_name[2]+"_a.txt > ZTT_Combined.txt"
    #command_h_mu_e = "cp " + combined_card_name[0]+".txt ZTT_Combined.txt"
    os.system(command_h_mu_e)
    
    Whether_Hybrid=False
    Whether_Asymptotic=False
    Whether_Simple_Bayesian=True
    print("Whether_Hybrid: "+ str(Whether_Hybrid))
    print("Whether_Asymptotic: "+ str(Whether_Asymptotic))
    print("Whether_Simple_Bayesian: "+ str(Whether_Simple_Bayesian))
    
    
    if(Whether_Hybrid):
            
            rvals = []
            
            command_run = "combine -M AsymptoticLimits ZTT_Combined.txt --cl 0.9 -t -1  > out_asym.txt"
            os.system(command_run)
            with open("out_asym.txt", 'r') as f1:
                for line in f1:
                    if '2.5%' in line:
                        linsp = line.split()
                        rvals.append(float(linsp[-1]))
                    if '16.0%' in line:
                        linsp = line.split()
                        rvals.append(float(linsp[-1]))
                    if '50.0%' in line:
                        linsp = line.split()
                        rvals.append(float(linsp[-1]))
                    if '84.0%' in line:
                        linsp = line.split()
                        rvals.append(float(linsp[-1]))
                    if '97.5%' in line:
                        linsp = line.split()
                        rvals.append(float(linsp[-1]))
            rvals_min = [r_i * 0.55 for r_i in rvals]
            rvals_max = [r_i * 1.45 for r_i in rvals]
            print("rvals: ")
            print(rvals)
            
            
            print("Running Sigma -2")
            command_run = "combine -M HybridNew ZTT_Combined.txt --cl 0.9 -t -1 --rMin "+"{:.3f}".format(rvals_min[0])+" --rMax "+"{:.3f}".format(rvals_max[0])+" --rRelAcc=0.01 --expectedFromGrid=0.025  > out1.txt"
            os.system(command_run)
            print("Running Sigma -1")
            command_run = "combine -M HybridNew ZTT_Combined.txt --cl 0.9 -t -1 --rMin "+"{:.3f}".format(rvals_min[1])+" --rMax "+"{:.3f}".format(rvals_max[1])+" --rRelAcc=0.01 --expectedFromGrid=0.16  > out2.txt"
            os.system(command_run)
            print("Running Sigma Median")
            command_run = "combine -M HybridNew ZTT_Combined.txt --cl 0.9 -t -1 --rMin "+"{:.3f}".format(rvals_min[2])+" --rMax "+"{:.3f}".format(rvals_max[2])+" --rRelAcc=0.01 --expectedFromGrid=0.5  > out3.txt"
            os.system(command_run)
            print("Running Sigma +1")
            command_run = "combine -M HybridNew ZTT_Combined.txt --cl 0.9 -t -1 --rMin "+"{:.3f}".format(rvals_min[3])+" --rMax "+"{:.3f}".format(rvals_max[3])+" --rRelAcc=0.01 --expectedFromGrid=0.84  > out4.txt"
            os.system(command_run)
            print("Running Sigma +2")
            command_run = "combine -M HybridNew ZTT_Combined.txt --cl 0.9 -t -1 --rMin "+"{:.3f}".format(rvals_min[4])+" --rMax "+"{:.3f}".format(rvals_max[4])+" --rRelAcc=0.01 --expectedFromGrid=0.975  > out5.txt"
            os.system(command_run)
            
            for i_std in range(5):
                    with open("out"+str(int(i_std+1))+".txt", 'r') as f1:
                            for line in f1:
                                    if 'Limit: r <' in line:
                                            linsp = line.split()
                                            ZTT_lim[i_std].append(float(linsp[3]))
            print(ZTT_lim)
                                    
    if(Whether_Asymptotic):
            command_run = "combine -M AsymptoticLimits ZTT_Combined.txt --cl 0.9 -t -1  > out1.txt"
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
            print(ZTT_lim)
            
            
    if(Whether_Simple_Bayesian):
            command_run = "combine -M BayesianSimple ZTT_Combined.txt --cl 0.9 -t 100  > out1.txt"
            os.system(command_run)
            with open("out1.txt", 'r') as f1:
                for line in f1:
                    if 'median expected' in line:
                        linsp = line.split()
                        ZTT_lim[2].append(float(linsp[5]))
                    if '68% expected band' in line:
                        linsp = line.split()
                        ZTT_lim[1].append(float(linsp[4]))
                        ZTT_lim[3].append(float(linsp[8]))
                    if '95% expected band' in line:
                        linsp = line.split()
                        ZTT_lim[0].append(float(linsp[4]))
                        ZTT_lim[4].append(float(linsp[8]))
            print(ZTT_lim)
    
    
ZTT_lim = [ZTT_lim_mid,ZTT_lim_low_2,ZTT_lim_low,ZTT_lim_up,ZTT_lim_up_2]
print(ZTT_lim)
