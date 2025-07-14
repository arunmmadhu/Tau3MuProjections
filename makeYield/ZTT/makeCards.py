#!/usr/bin/env python3

import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
from ROOT import RooRealVar, RooFormulaVar, RooExponential, RooDataHist, RooArgList, RooAddPdf, RooFit, RooDataSet, RooGenericPdf, RooBifurGauss
import subprocess # to execute shell command
import argparse
import numpy as np
import CMS_lumi, tdrstyle
#from CMSStyle import CMS_lumi
import CMSStyle
import os
import re
from array import array


# CMS style
CMSStyle.cmsText = "CMS"
CMSStyle.extraText = "        Work in progress"
CMSStyle.relPosX = 0.070
CMSStyle.outOfFrame = False
CMSStyle.alignX_ = 1
CMSStyle.relPosX    = 0.04
#CMSStyle.relPosY    = 0.025
tdrstyle.setTDRStyle()

# references for T, B, L, R
H_ref = 800; 
W_ref = 800; 
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.18*H_ref 
L = 0.15*W_ref
R = 0.04*W_ref


class makeCards:
        
        def __init__(self):
                
                self.bdt_cv = None
                
                self.bgausmeanMC = None
                self.bgaussigmaMC_a = None
                self.bgaussigmaMC_b = None
                self.bgaus_distMC = None
                self.BDTNorm_MC = None
                self.BDT_distribution_MC = None
                self.MCSelector = None
                self.fullmc_unweighted = None
                self.fullmc = None
                
                self.a = None
                self.b = None
                self.c = None
                self.d = None
                self.quadratic = None
                self.expModel = None
                self.BDTNorm = None
                self.BDT_distribution_MC = None
                
                self.results_pdf = None
                self.results_mcpdf = None
                
                self.bkg_PFcut_norm = None
        
        
        
        def FitBDT(self,datafile_for_norm,datafile_for_shape,categ):
                
                fit_range_lo = 1.6
                fit_range_hi = 2.0
                
                signal_range_lo = 1.74
                signal_range_hi = 1.81
                
                treeName=''
                signalnorm = 1.0
                cat_label = ""
                if(categ == 'taue'):
                        treeName  = 'ztau3mutaue'
                        signalnorm = 0.00000856928
                        cat_label = r"$\tau_{e}$"
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                        cat_label = r"$\tau_{\mu}$"
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,1-prong}$"
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                        cat_label = r"$\tau_{h,3-prong}$"
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                        cat_label = "Inclusive"
                
                MiniTreeFile = ROOT.TFile.Open(datafile_for_shape)
                MiniTreeFile.cd()
                tree = MiniTreeFile.Get(treeName)
                
                MiniTreeFile_ForNorm = ROOT.TFile.Open(datafile_for_norm)
                MiniTreeFile_ForNorm.cd()
                tree_norm = MiniTreeFile_ForNorm.Get(treeName)
                
                tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                self.bdt_cv          = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 100)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 100)
                event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
                category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
                isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
                scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  signalnorm)  
                
                
                
                variables = ROOT.RooArgSet()
                variables.add(tripletMass)
                variables.add(self.bdt_cv)
                variables.add(dimu_OS1)
                variables.add(dimu_OS2)
                variables.add(event_weight)
                variables.add(category)
                variables.add(isMC)
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)&"
                omegavetoes="fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                
                # For fitting BDT Output in Data
                
                BDT_Score_Min=-0.3
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector', phivetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fulldata_shape = RooDataSet('fulldata_shape', 'fulldata_shape', tree,  variables, BlindDataSelector)
                fulldata_norm = RooDataSet('fulldata_norm', 'fulldata_norm', tree_norm,  variables, BlindDataSelector)
                #print("shape tree norm: ",fulldata_shape.numEntries()," and main tree norm: ",fulldata_norm.numEntries())
                data_scale = ROOT.RooRealVar('data_scale', 'data_scale', 1.0*fulldata_norm.numEntries()/fulldata_shape.numEntries())
                self.bkg_PFcut_norm = 1.0*fulldata_norm.numEntries()/fulldata_shape.numEntries()
                #dataset_vars_data = fulldata_shape.get()
                #dataset_vars_data.add(data_scale)
                
                fulldata = fulldata_shape
                
                #fulldata = RooDataSet('fulldata', 'fulldata', fulldata_shape,  dataset_vars_data, "",'data_scale')
                
                self.bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0);
                
                self.a = RooRealVar("a", "a", 1.0, 0.0, 10.0)
                self.b = RooRealVar("b", "b", 1.0, -10.0, 10.0)
                self.c = RooRealVar("c", "c", 1.0, -10.0, 10.0)
                self.d = RooRealVar("d", "d", 1.0, -10.0, 10.0)
                
                #quadratic = RooFormulaVar("quadratic", "a + b*self.bdt_cv + c*self.bdt_cv*self.bdt_cv", RooArgList(a, b, c, self.bdt_cv))
                #expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(quadratic)) #Exponential of the quadratic polynomial
                
                self.quadratic = RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv + d*bdt_cv*bdt_cv*bdt_cv", RooArgList(self.a, self.b, self.c, self.d, self.bdt_cv))
                self.expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(self.quadratic)) #Exponential of the cubic polynomial
                
                self.BDTNorm = RooRealVar("BDTNorm", "BDTNorm", 500.0, 0.1, 1000000000)
                self.BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(self.expModel), RooArgList(self.BDTNorm))
                #BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(quadratic), RooArgList(BDTNorm))
                
                self.results_pdf = self.BDT_distribution.fitTo(fulldata, RooFit.Range('BDT_Fit_Range'), RooFit.Save())
                
                
                
                
                # For fitting BDT Output in Signal

                self.MCSelector = RooFormulaVar('MCSelector', 'MCSelector', phivetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi) , RooArgList(variables))
                
                #signal always uses the normalized/final PF tree
                self.fullmc_unweighted = RooDataSet('mc', 'mc', tree_norm, variables, self.MCSelector)
                dataset_vars = self.fullmc_unweighted.get()
                dataset_vars.add(scale)
                
                self.fullmc = RooDataSet('mc', 'mc', self.fullmc_unweighted, dataset_vars, "",'scale')

                self.bdt_cv.setRange("BDT_MC_Fit_Range", -1.0, 1.0);

                self.bgausmeanMC = RooRealVar("bgausmeanMC", "bgausmeanMC", 0.5, 0.0, 0.9)
                self.bgaussigmaMC_a = RooRealVar("bgaussigmaMC_a", "bgaussigmaMC_a", 0.2, 0.000001, 1.0)
                self.bgaussigmaMC_b = RooRealVar("bgaussigmaMC_b", "bgaussigmaMC_b", 0.2, 0.000001, 1.0)

                self.bgaus_distMC = RooBifurGauss("bgaus_distMC", "bgaus dist MC", self.bdt_cv, self.bgausmeanMC, self.bgaussigmaMC_a, self.bgaussigmaMC_b)

                self.BDTNorm_MC = RooRealVar("BDTNorm_MC", "BDTNorm_MC", 500.0, 0.1, 50000)
                self.BDT_distribution_MC = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(self.bgaus_distMC), RooArgList(self.BDTNorm_MC))

                self.results_mcpdf = self.BDT_distribution_MC.fitTo(self.fullmc, RooFit.Range('BDT_MC_Fit_Range'), RooFit.Save())
                #results_mcpdf.Print()


                # Plot BDT for data and MC
                
                frame1 = self.bdt_cv.frame()
                frame1.SetTitle('')
                frame1.GetXaxis().SetTitle("BDT score, " + cat_label)
                
                frame2 = self.bdt_cv.frame()
                frame2.SetTitle('')
                frame2.GetXaxis().SetTitle("BDT score, " + cat_label)
                
                nbins = 100
                
                self.fullmc.plotOn(frame1, 
                              ROOT.RooFit.Binning(nbins), 
                              ROOT.RooFit.XErrorSize(0), 
                              ROOT.RooFit.LineWidth(2),
                              ROOT.RooFit.MarkerStyle(6),
                              ROOT.RooFit.MarkerColor(ROOT.kRed),
                              ROOT.RooFit.MarkerSize(0.75),
                              ROOT.RooFit.FillColor(ROOT.kCyan + 2)
                )
                fulldata.plotOn(frame2, 
                        ROOT.RooFit.Binning(nbins),
                        ROOT.RooFit.MarkerStyle(20),
                        ROOT.RooFit.MarkerColor(ROOT.kBlack), 
                        ROOT.RooFit.MarkerSize(0.75))
                
                self.BDT_distribution.plotOn(frame2, ROOT.RooFit.LineColor(ROOT.kBlue))
                
                # Create canvas
                canvas = ROOT.TCanvas("bdt_canvas", "BDT Canvas", 800, 800)
                canvas.cd()
                canvas.SetLeftMargin(L / W)
                canvas.SetRightMargin(R / W)
                canvas.SetTopMargin(T / H)
                canvas.SetBottomMargin(B / H)
                
                # Plot MC
                frame1.Draw()
                canvas.SetTicks(1, 1)
                CMSStyle.CMS_lumi(canvas, 5, 0)
                canvas.Update()
                frame1.Draw('sameaxis')
                canvas.SaveAs('bdt_fit_mc_' + categ + '.png')
                canvas.Clear()
                
                # Plot data + model (log scale)
                canvas.SetLogy()
                frame2.SetMinimum(1e-2)
                frame2.Draw()
                canvas.SetTicks(1, 1)
                CMSStyle.CMS_lumi(canvas, 5, 0)
                canvas.Update()
                frame2.Draw('sameaxis')
                canvas.SaveAs('bdt_fit_bkg_' + categ + '.png')
                canvas.SetLogy(0)



                #print("Certain BDT cut: ", self.fullmc.reduce('bdt_cv > 0.5').sumEntries())

                return frame1, frame2, self.fullmc, fulldata
                
                
        
        
        #To create datacards
        def MakeLumiScanCards(self,lumi,categ,analyzed_lumi):
                
                command_recreate_categ_dir = "rm -r lumi_limit_scans/{0}; mkdir lumi_limit_scans/{0}".format(categ)
                os.system(command_recreate_categ_dir)
                
                fit_range_lo = 1.6
                fit_range_hi = 2.0
                
                signal_range_lo = 1.74
                signal_range_hi = 1.81
                
                signalnorm = 1.0
                if(categ == 'taue'):
                        treeName  = 'ztau3mutaue'
                        signalnorm = 0.00000856928
                if(categ == 'taumu'):
                        treeName = 'ztau3mutaumu'
                        signalnorm = 0.00000822810
                if(categ == 'tauhA'):
                        treeName = 'ztau3mutauh_A'
                        signalnorm = 0.00000815958
                if(categ == 'tauhB'):
                        treeName = 'ztau3mutauh_B'
                        signalnorm = 0.00000815958
                if(categ == 'all'):
                        treeName   = 'ztautau'
                        signalnorm = 0.00000824176
                
                exp_fact = (signal_range_hi-signal_range_lo)/(fit_range_hi-fit_range_lo-(signal_range_hi-signal_range_lo))
                
                for point in bdt_points:  # For loop for bdt cuts in range [X_min;X_max]
                
                    self.bdt_cv.setRange("Integral_Range", point, 1.0)
                    
                    command_recreate_lumidir = "rm -r lumi_limit_scans/{0}/BDT_point_{1}; mkdir lumi_limit_scans/{0}/BDT_point_{1}".format(categ, str(point))
                    os.system(command_recreate_lumidir)
                    
                    BDT_cut = 'bdt_cv > %s' % point
                    
                    MC_dataset_with_BDT_cut = self.fullmc.reduce(BDT_cut)
                    
                    for lu_no in range(len(lumi)):
                        
                        lu = lumi[lu_no]
                        print(lu)
                        
                        sig_est = MC_dataset_with_BDT_cut.sumEntries() * lu/analyzed_lumi
                        #sig_est = self.BDTNorm_MC.getVal() * (self.BDT_distribution_MC.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range").getVal() ) * lu/analyzed_lumi
                        #bkg_est =bkg in signal region; sb_est = bkg in sideband
                        
                        bkg_est = self.bkg_PFcut_norm * exp_fact * self.BDTNorm.getVal() * (self.BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range").getVal() ) * lu/analyzed_lumi
                        
                        sb_est = self.bkg_PFcut_norm * self.BDTNorm.getVal() * (self.BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range").getVal() ) * lu/analyzed_lumi
                        
                        exp_fact_var = ROOT.RooRealVar("exp_fact_var", "exp_fact_var", exp_fact)
                        bdt_norm_var = ROOT.RooRealVar("bdt_norm_var", "bdt_norm_var", self.BDTNorm.getVal())
                        lumi_ratio = ROOT.RooRealVar("lumi_ratio", "lumi_ratio", lu / analyzed_lumi)
                        bkg_PF_normal_ratio = ROOT.RooRealVar("bkg_PF_normal_ratio", "bkg_PF_normal_ratio", self.bkg_PFcut_norm)
                        combined_norm_bkg = ROOT.RooProduct("combined_norm", "combined_norm", ROOT.RooArgList(exp_fact_var, bdt_norm_var, lumi_ratio, bkg_PF_normal_ratio))
                        bdt_integral_bkg = self.BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range")
                        bkg_est_var_withNorm = ROOT.RooProduct("bkg_est_var_withNorm", "bkg_est_var_withNorm", ROOT.RooArgList(combined_norm_bkg, bdt_integral_bkg))
                        
                        mc_sum = MC_dataset_with_BDT_cut.sumEntries()
                        mc_sum_var = ROOT.RooRealVar("mc_sum_var", "mc_sum_var", mc_sum)
                        sig_est_var_withNorm = ROOT.RooProduct("sig_est_var_withNorm", "sig_est_var_withNorm", ROOT.RooArgList(lumi_ratio, mc_sum_var))
                        
                        print( "bdt point: ", point," sig_est: ", sig_est, " bkg_est: ", bkg_est, " bkg_est_var_withNorm: ", bkg_est_var_withNorm.getVal(), " bkg_est_var_withNorm error: ", bkg_est_var_withNorm.getPropagatedError(self.results_pdf), " sig_est_var_withNorm ", sig_est_var_withNorm.getVal() )
                        
                        self.results_pdf.floatParsFinal().Print()
                        
                        exp_fact = (signal_range_hi-signal_range_lo)/(fit_range_hi-fit_range_lo-(signal_range_hi-signal_range_lo))
                        #print('   exp_fact   ', exp_fact)
                        exp_fact_different_pdf = getExtrapFactor('unfixed_exp', categ, point)
                        #print('   exp_fact_different_pdf   ', exp_fact_different_pdf )
                        
                        exp_uncert = 1.0 + abs(exp_fact_different_pdf-exp_fact)/exp_fact
                        
                        #print("bdt point: ", point," sig_est: ", sig_est, " bkg_est: ", bkg_est, " sb_est: ", sb_est, " exp_uncert: ", exp_uncert)
                        
                        command_mod_card = "python3 card_modifiers/Card_Mod.py --categ " + str(categ) + " --sig_exp " + str(sig_est) + " --bkg_exp " + str(bkg_est) + " --sb_exp " + str(sb_est) + " --ext_unc " + str(exp_uncert) + " --bkg_exp_err " + str(bkg_est_var_withNorm.getPropagatedError(self.results_pdf))
                        
                        
                        os.system(command_mod_card)
                        print(command_mod_card)
                        
                        command_copy_dc = "cp modified_dc_{0}.txt lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt".format(categ, str(point), str(lu))
                        os.system(command_copy_dc)
                        
        def CombineSubcategories(self,categ):
                pattern = re.compile(r'dc_(\d+(?:\.\d+)?)\.txt')
                
                
                files_dict = {}
                
                sub_cats = ['taue','taumu','tauhA','tauhB']
                
                dirs = []
                
                for subcatno in range(len(sub_cats)):
                        dirs.append(sub_cats[subcatno]+'/datacards_modified')
                
                os.system('mkdir combined')
                os.system('mkdir combined/datacards_modified')
                
                
                
                for directory in dirs:
                    for filename in os.listdir(directory):
                        match = pattern.match(filename)
                        #print("Matched: ",match)
                        if match:
                            number = match.group(1)
                            if number not in files_dict:
                                files_dict[number] = {}
                            files_dict[number][directory] = os.path.join(directory, filename)
                
                
                for number, files in files_dict.items():
                    if len(files) == len(dirs):  # Make sure all directories have this number
                        file1 = files['taue/datacards_modified']
                        file2 = files['taumu/datacards_modified']
                        file3 = files['tauhA/datacards_modified']
                        file4 = files['tauhB/datacards_modified']
                        output_file = "dc_%s.txt" % str(number)
                        
                
                #        command = "combineCards.py %s %s %s > %s"% (str(file1) ,str(file2) ,str(file3), str(output_file))  # rm ZTT for now
                        command = "combineCards.py %s %s %s %s > %s"% (str(file1),str(file2),str(file3) ,str(file4), str(output_file))
                        
                        print("Running command: %s" % (command))
                        os.system(command)
                        os.system('mv dc_*txt combined/datacards_modified')
                    else:
                        print("Skipping number %s, not all directories have the file." % str(number))
        




def executeDataCards_onCondor(lumi,categories,Whether_Hybrid,bdt_points):
        
        Cat_No = len(categories)
        
        for cat in range(Cat_No):
        
                for point in bdt_points:  # For loop for bdt cuts in range [X_min;X_max]
                        
                        for lu_no in range(len(lumi)):
                            lu = lumi[lu_no]
                            print('Luminosity: ', lu)
                                
                            if(Whether_Hybrid):
                                    
                                    #command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 5 --expectedFromGrid 0.5 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+'_'+categories[cat]+'_'+str(point),"lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt".format(categories[cat], str(point), str(lu)),str(lu)+'_'+categories[cat]+'_'+str(point))
                                    
                                    command_run = "combineTool.py -M HybridNew --generateNuisances=1 --generateExternalMeasurements=0 --fitNuisances=1 --testStat LHC -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -T 10000 --expectedFromGrid 0.5 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+'_'+categories[cat]+'_'+str(point),"lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt".format(categories[cat], str(point), str(lu)),str(lu)+'_'+categories[cat]+'_'+str(point))
                                    
                                    print("Run:   ", command_run)
                                    os.system(command_run)
                                    
                            else:
                                    
                                    command_run = "combineTool.py -M AsymptoticLimits --run blind  --cl 0.90 -n %s -d %s  --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name AsymTest%s " % (str(lu)+'_'+categories[cat]+'_'+str(point),"lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt".format(categories[cat], str(point), str(lu)),str(lu)+'_'+categories[cat]+'_'+str(point))
                                    
                                    #command_run = "combineTool.py -M AsymptoticLimits --run blind  --cl 0.90 -n %s -d %s  " % (str(lu)+'_'+categories[cat]+'_'+str(point),"lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt".format(categories[cat], str(point), str(lu)))
        
                                    print("Run:   ", command_run)
                                    os.system(command_run)
                                    
def CalculateUL_fromDatacard(lumi, categories, Whether_Hybrid, bdt_points):

    Cat_No = len(categories)

    for cat in range(Cat_No):
        category_dir = categories[cat]

        # Make sure category folder exists
        os.makedirs(category_dir, exist_ok=True)

        for lu in lumi:
            output_file_path = f"{category_dir}/limits_lumi_{lu}.txt"
            with open(output_file_path, 'w') as out_txt:

                for point in bdt_points:
                    print('Luminosity: ', lu, ' of category: ',categories[cat])

                    dc_path = f"lumi_limit_scans/{category_dir}/BDT_point_{point}/dc_{lu}.txt"
                    if not os.path.isfile(dc_path):
                        print(f"Missing datacard: {dc_path}")
                        continue

                    with open(dc_path, "r") as file:
                        sig = 0.0001
                        bkg = 0.0001
                        bkg_err = 0.0000

                        for linenum, line in enumerate(file):
                            if linenum == 14:
                                x1 = line.split()
                                sig = float(x1[1])
                                bkg = float(x1[2])
                                if sig < 0.0001:
                                        sig = 0.0001
                            
                            if linenum == 0:
                                x1 = line.split()
                                bkg_err = float(x1[2])
                            
                    print("sig:", sig, "bkg:", bkg)

                    cmd_ul = f"python3 ../../CLs_UL_Calculator_Integral.py {sig} {bkg} > out_UL_Calc.txt"
                    os.system(cmd_ul)

                    limit_val = 10000
                    with open("out_UL_Calc.txt") as f:
                        for line in f:
                            if "r_val" in line:
                                limit_val = line.split()[-1]
                                print(f"Limit UL_Calc: {limit_val}")

                    # Write result to the output file
                    out_txt.write(f"bdt {point}   sig {sig:.6f}   bkg {bkg:.6f}   limit {limit_val}  bkg_err {bkg_err}\n")

# GET limits from root file
def getLimits(file_name):
    file = TFile.Open(file_name)
    if not file or file.IsZombie():
        print(f"[ERROR] Could not open file: {file_name}")
        return []

    tree = file.Get("limit")
    if not tree or not isinstance(tree, TTree):
        print(f"[ERROR] 'limit' TTree not found or invalid in file: {file_name}")
        file.Close()
        return []

    limits = []
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        try:
            limits.append(tree.limit)
            print(">>>   %.2f" % limits[-1])
        except AttributeError:
            print(f"[ERROR] Entry {i} in tree is missing 'limit' branch")
            limits.append(999999.)

    file.Close()
    return limits[:6]
    
#check if combine output file is valid
def is_valid_root_file(filepath):
        if not os.path.exists(filepath):
            return False
        file_size = os.path.getsize(filepath)
        if file_size < 2 * 1024:
            return False
        f = ROOT.TFile.Open(filepath)
        if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
            return False
        if f.GetListOfKeys().GetSize() == 0:
            return False
        f.Close()
        return True
        
# PLOT upper limits
def ReadAndCopyMinimumBDTCard(lumi,categories,Whether_Hybrid,bdt_points):
 
    N = len(lumi)
    Cat_No = len(categories)
    
    label=[None] * Cat_No
    
    for cat in range(Cat_No):
    
            categ = categories[cat]
            
            cmd1 = "mkdir {0}".format(categ)
            os.system(cmd1)
            
            cmd2 = "rm -rf {0}/datacards_modified/; mkdir {0}/datacards_modified/;".format(categ)
            os.system(cmd2)
            
            if(categ == 'taue'):
                    analyzed_lumi = 59.83
            if(categ == 'taumu'):
                    analyzed_lumi = 59.83
            if(categ == 'tauhA'):
                    analyzed_lumi = 59.83
            if(categ == 'tauhB'):
                    analyzed_lumi = 59.83
            if(categ == 'all'):
                    analyzed_lumi = 59.83
            
            text_limits=open("TextLimits_"+categ+".txt","w")
            
            if(Whether_Hybrid):
                    limits_read = []
                    for point in bdt_points:  # For loop for bdt cuts in range [X_min;X_max]
                            limits_read_row = []
                            for i in range(N):
                                file_name = "higgsCombine"+str(lumi[i])+'_'+categories[cat]+'_'+str(point)+".HybridNew.mH120.quant0.500.root"
                                limit = [1000000.0,1000000.0,1000000.0,1000000.0,1000000.0]
                                
                                if is_valid_root_file(file_name):
                                        limit = getLimits(file_name)
                                        
                                        #  Check why some limits only have 1 value
                                        if len(limit)<1:
                                                limit = [1000000.0,1000000.0,1000000.0,1000000.0,1000000.0]
                                
                                print(" cat: ",categories[cat]," lumi: ",lumi[i]," bdt point: ",point," Limit: ",limit[0])
                                limits_read_row.append(limit[0])
                                
                            limits_read.append(limits_read_row)
                    
                    # Getting the BDT cut at minimum limit
                    transposed_matrix = list(zip(*limits_read))
                    min_values = [min(column) for column in transposed_matrix]
                    print("limits_read values: ",limits_read)
                    print("min values: ",min_values)
                    min_indices = []
                    for col_index, min_val in enumerate(min_values):
                            # Find the row in the original matrix where the minimum value is located in this column
                            for row_index, row in enumerate(limits_read):
                                if row[col_index] == min_val:
                                    min_indices.append((row_index, col_index))
                                    break  # Stop searching once the minimum is found
                    for bdt_index, lumi_index in min_indices:
                            #print(f"Minimum value found at row {row_index}, column {col_index}")
                            text_limits.write("bdt %.2f   lumi %.2f     median exp %.2f\n"%(bdt_points[bdt_index],lumi[lumi_index],limits_read[bdt_index][lumi_index]))
                            
                            command_copy_dc = "cp  lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt {0}/datacards_modified/dc_{2}.txt".format(categ, str(bdt_points[bdt_index]), str(lumi[lumi_index]))
                            print("Copy command: ",command_copy_dc)
                            os.system(command_copy_dc)
                    
                    output_dir = f"limit_scans_by_lumi/{categories[cat]}"
                    os.makedirs(output_dir, exist_ok=True)   
                    
                    for i in range(N):
                            y_vals = [limits_read_row[i] for limits_read_row in limits_read]
                            x_vals = bdt_points
                            
                            graph = TGraph(len(x_vals), array('d', x_vals), array('d', y_vals))
                            graph.SetLineColor(1)
                            graph.SetLineWidth(2)
                            graph.SetMarkerStyle(20)
                            
                            c = TCanvas(f"c_{i}", f"BDT Scan Lumi {lumi[i]}", W, H)
                            c.SetFillColor(0)
                            c.SetBorderMode(0)
                            c.SetFrameFillStyle(0)
                            c.SetFrameBorderMode(0)
                            c.SetLeftMargin(L / W)
                            c.SetRightMargin(R / W)
                            c.SetTopMargin(T / H)
                            c.SetBottomMargin(B / H)
                            c.SetGrid()
                            c.cd()
                            
                            frame = c.DrawFrame(min(x_vals), 0.0, max(x_vals)*1.1, max(y_vals)*1.2)
                            frame.GetYaxis().SetTitle("B(#tau #rightarrow #mu#mu#mu) UL (10^{-7})")
                            frame.GetXaxis().SetTitle("MVA cut value")
                            frame.GetXaxis().SetTitleSize(0.05)
                            frame.GetYaxis().SetTitleSize(0.05)
                            frame.GetXaxis().SetLabelSize(0.04)
                            frame.GetYaxis().SetLabelSize(0.04)
                            frame.GetYaxis().SetTitleOffset(0.9)
                            frame.GetXaxis().SetNdivisions(508)
                            
                            frame.SetMinimum(0.0)
                            frame.SetMaximum(30.0)
                            
                            frame.GetXaxis().SetLimits(min(x_vals)*0.95 ,max(x_vals)*1.1)
                            
                            graph.Draw("PL same")
                            
                            legend = TLegend(0.15, 0.75, 0.5, 0.85)
                            legend.SetBorderSize(0)
                            legend.SetFillStyle(0)
                            legend.SetTextSize(0.041)
                            legend.SetTextFont(42)
                            legend.AddEntry(graph, "Expected median limit", "l")
                            legend.Draw()
                            
                            latex = TLatex()
                            latex.SetNDC()
                            latex.SetTextSize(0.04)
                            latex.SetTextFont(42)
                            text = f"Category: Z#rightarrow#tau#tau_{{3#mu}}, L = {lumi[i]} fb^{{-1}}"
                            latex.DrawLatex(0.15, 0.88, text)
                            
                            c.Update()
                            c.SaveAs(f"{output_dir}/limit_scan_lumi_{lumi[i]}.png")
                            c.Close()
                    
                    #print(min_indices)
            
            else:
                    limits_read = []
                    for point in bdt_points:  # For loop for bdt cuts in range [X_min;X_max]
                            limits_read_row = []
                            for i in range(N):
                                file_name1 = "higgsCombine"+str(lumi[i])+'_'+categories[cat]+'_'+str(point)+".AsymptoticLimits.mH120.root"
                                
                                limit = [1000000.0,1000000.0,1000000.0,1000000.0,1000000.0]
                                
                                if is_valid_root_file(file_name1):
                                        limit = getLimits(file_name1)
                                        
                                        #  Check why some limits only have 1 value
                                        if len(limit)<5:
                                                limit = [1000000.0,1000000.0,1000000.0,1000000.0,1000000.0]
                                
                                print(" cat: ",categories[cat]," lumi: ",lumi[i]," bdt point: ",point," Limit: ",limit[2])
                                limits_read_row.append(limit[2])
                            limits_read.append(limits_read_row)
                    
                    # Getting the BDT cut at minimum limit
                    transposed_matrix = list(zip(*limits_read))
                    min_values = [min(column) for column in transposed_matrix]
                    min_indices = []
                    for col_index, min_val in enumerate(min_values):
                            # Find the row in the original matrix where the minimum value is located in this column
                            for row_index, row in enumerate(limits_read):
                                if row[col_index] == min_val:
                                    min_indices.append((row_index, col_index))
                                    break  # Stop searching once the minimum is found
                    for bdt_index, lumi_index in min_indices:
                            #print(f"Minimum value found at row {row_index}, column {col_index}")
                            text_limits.write("bdt %.2f   lumi %.2f     median exp %.2f\n"%(bdt_points[bdt_index],lumi[lumi_index],limits_read[bdt_index][lumi_index]))
                            
                            command_copy_dc = "cp  lumi_limit_scans/{0}/BDT_point_{1}/dc_{2}.txt {0}/datacards_modified/dc_{2}.txt".format(categ, str(bdt_points[bdt_index]), str(lumi[lumi_index]))
                            print("Copy command: ",command_copy_dc)
                            os.system(command_copy_dc)
                            
                    output_dir = f"limit_scans_by_lumi/{categories[cat]}"
                    os.makedirs(output_dir, exist_ok=True)   
                    
                    for i in range(N):
                            y_vals = [limits_read_row[i] for limits_read_row in limits_read]
                            x_vals = bdt_points
                            
                            graph = TGraph(len(x_vals), array('d', x_vals), array('d', y_vals))
                            graph.SetLineColor(1)
                            graph.SetLineWidth(2)
                            graph.SetMarkerStyle(20)
                            
                            c = TCanvas(f"c_{i}", f"BDT Scan Lumi {lumi[i]}", W, H)
                            c.SetFillColor(0)
                            c.SetBorderMode(0)
                            c.SetFrameFillStyle(0)
                            c.SetFrameBorderMode(0)
                            c.SetLeftMargin(L / W)
                            c.SetRightMargin(R / W)
                            c.SetTopMargin(T / H)
                            c.SetBottomMargin(B / H)
                            c.SetGrid()
                            c.cd()
                            
                            frame = c.DrawFrame(min(x_vals) * 0.95, 0.0, max(x_vals)*1.1, max(y_vals)*1.2)
                            frame.GetYaxis().SetTitle("B(#tau #rightarrow #mu#mu#mu) UL (10^{-7})")
                            frame.GetXaxis().SetTitle("MVA cut value")
                            frame.GetXaxis().SetTitleSize(0.05)
                            frame.GetYaxis().SetTitleSize(0.05)
                            frame.GetXaxis().SetLabelSize(0.04)
                            frame.GetYaxis().SetLabelSize(0.04)
                            frame.GetYaxis().SetTitleOffset(0.9)
                            frame.GetXaxis().SetNdivisions(508)
                            
                            #Vary the frame ymax from 30 to 5 as lumi goes from 59.83 to 3000
                            frame_max = 30.0 + (15.0 - 30.0) * ((lumi[i] - 59.83) / (3000.0 - 59.83))
                            
                            frame.SetMinimum(0.0)
                            frame.SetMaximum(frame_max)
                            
                            graph.Draw("PL same")
                            
                            legend = TLegend(0.15, 0.75, 0.5, 0.85)
                            legend.SetBorderSize(0)
                            legend.SetFillStyle(0)
                            legend.SetTextSize(0.041)
                            legend.SetTextFont(42)
                            legend.AddEntry(graph, "Expected median limit", "l")
                            legend.Draw()
                            
                            latex = TLatex()
                            latex.SetNDC()
                            latex.SetTextSize(0.04)
                            latex.SetTextFont(42)
                            text = f"Category: Z#rightarrow#tau#tau_{{3#mu}}, L = {lumi[i]} fb^{{-1}}"
                            latex.DrawLatex(0.15, 0.88, text)
                            
                            c.Update()
                            c.SaveAs(f"{output_dir}/limit_scan_lumi_{lumi[i]}.png")
                            c.Close()
                            
                    #print(min_indices)

# PLOT upper limits
def ReadAndCopyMinimumBDTCard_usingUL(lumi,categories,Whether_Hybrid,bdt_points):
 
    N = len(lumi)
    Cat_No = len(categories)
    
    label=[None] * Cat_No
    
    for cat in range(Cat_No):
    
            categ = categories[cat]
            
            # Make sure category folder exists
            os.makedirs(categ, exist_ok=True)
            
            cmd2 = "rm -rf {0}/datacards_modified/; mkdir {0}/datacards_modified/;".format(categ)
            os.system(cmd2)
            
            if(categ == 'taue'):
                    analyzed_lumi = 59.83
            if(categ == 'taumu'):
                    analyzed_lumi = 59.83
            if(categ == 'tauhA'):
                    analyzed_lumi = 59.83
            if(categ == 'tauhB'):
                    analyzed_lumi = 59.83
            if(categ == 'all'):
                    analyzed_lumi = 59.83
            
            os.makedirs(f"{categ}/datacards_modified", exist_ok=True)
            
            
            N = len(lumi)
            M = len(bdt_points)
        
            limits_matrix = []  # [bdt_index][lumi_index]
        
            for bdt_point in bdt_points:
                limits_row = []
                for i in range(N):
                    lumi_val = lumi[i]
                    filename = f"{categ}/limits_lumi_{lumi_val}.txt"
                    limit_val = 1e6
        
                    if os.path.isfile(filename):
                        with open(filename) as f:
                            for line in f:
                                    tokens = line.strip().split()
                                    if len(tokens) < 10:
                                        continue
                                    try:
                                        bdt_val = float(tokens[1])
                                        if abs(bdt_val - bdt_point) < 1e-5:
                                            limit_val = float(tokens[-3])
                                            break
                                    except:
                                        continue
                    else:
                        print(f"Missing file: {filename}")
                    
                    limits_row.append(limit_val)
                limits_matrix.append(limits_row)
        
            # Transpose: now limits[lumi_index][bdt_index]
            transposed = list(zip(*limits_matrix))
        
            # Open output text file
            text_limits = open(f"TextLimits_{categ}.txt", "w")
            
            output_dir = f"limit_scans_by_lumi_ZTT_UL/{categories[cat]}"
            os.makedirs(output_dir, exist_ok=True)  
        
            for lumi_index, row in enumerate(transposed):
                y_vals = list(row)
                x_vals = bdt_points
        
                min_limit = min(y_vals)
                bdt_index = y_vals.index(min_limit)
                best_bdt = bdt_points[bdt_index]
                lumi_val = lumi[lumi_index]
        
                print(f"Best BDT for lumi {lumi_val}: {best_bdt} with limit {min_limit}")
        
                # Copy best datacard
                src = f"lumi_limit_scans/{categ}/BDT_point_{best_bdt}/dc_{lumi_val}.txt"
                dst = f"{categ}/datacards_modified/dc_{lumi_val}.txt"
                os.makedirs(f"{categ}/datacards_modified", exist_ok=True)
                os.system(f"cp {src} {dst}")
                
                # Step 4: Read sig, bkg, bkg_err from line with best_bdt
                sig_val = bkg_val = bkg_err_val = -1.0
                filename = f"{categ}/limits_lumi_{lumi_val}.txt"
                if os.path.isfile(filename):
                        with open(filename) as f:
                            for line in f:
                                tokens = line.strip().split()
                                if len(tokens) < 10:
                                    continue
                                try:
                                    bdt_val = float(tokens[1])
                                    if abs(bdt_val - best_bdt) < 1e-5:
                                        sig_val = float(tokens[3])
                                        bkg_val = float(tokens[5])
                                        limit_val = float(tokens[7])
                                        bkg_err_val = float(tokens[9])
                                        break
                                except:
                                    continue
                
                
                
                # Write to text file
                text_limits.write("bdt %.4f   lumi %.2f   median_exp %.3f   sig %.4f   bkg %.4f   bkg_err %.4f\n" %(best_bdt, lumi_val, min_limit, sig_val, bkg_val, bkg_err_val))
                
                graph = TGraph(len(x_vals), array('d', x_vals), array('d', y_vals))
                graph.SetLineColor(1)
                graph.SetLineWidth(2)
                graph.SetMarkerStyle(20)
                
                c = TCanvas(f"c_{i}", f"BDT Scan Lumi {lumi[lumi_index]}", W, H)
                c.SetFillColor(0)
                c.SetBorderMode(0)
                c.SetFrameFillStyle(0)
                c.SetFrameBorderMode(0)
                c.SetLeftMargin(L / W)
                c.SetRightMargin(R / W)
                c.SetTopMargin(T / H)
                c.SetBottomMargin(B / H)
                c.SetGrid()
                c.cd()
                
                frame = c.DrawFrame(min(x_vals) * 0.95, 0.0, max(x_vals)*1.1, max(y_vals)*1.2)
                frame.GetYaxis().SetTitle("B(#tau #rightarrow #mu#mu#mu) UL (10^{-7})")
                frame.GetXaxis().SetTitle("MVA cut value")
                frame.GetXaxis().SetTitleSize(0.05)
                frame.GetYaxis().SetTitleSize(0.05)
                frame.GetXaxis().SetLabelSize(0.04)
                frame.GetYaxis().SetLabelSize(0.04)
                frame.GetYaxis().SetTitleOffset(1.2)
                frame.GetXaxis().SetNdivisions(508)
                
                #Vary the frame ymax from 30 to 5 as lumi goes from 59.83 to 3000
                frame_max = 30.0 + (15.0 - 30.0) * ((lumi[lumi_index] - 59.83) / (3000.0 - 59.83))
                
                frame.SetMinimum(0.0)
                frame.SetMaximum(frame_max)
                
                graph.Draw("PL same")
                
                legend = TLegend(0.19, 0.75, 0.5, 0.85)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.041)
                legend.SetTextFont(42)
                legend.AddEntry(graph, "Expected median limit", "l")
                legend.Draw()
                
                latex = TLatex()
                latex.SetNDC()
                latex.SetTextSize(0.04)
                latex.SetTextFont(42)
                text = f"Category: Z#rightarrow#tau#tau_{{3#mu}}, L = {lumi[lumi_index]} fb^{{-1}}"
                latex.DrawLatex(0.19, 0.86, text)
                
                CMSStyle.CMS_lumi(c, 5, 0)
                
                c.Update()
                c.SaveAs(f"{output_dir}/limit_scan_lumi_{lumi[lumi_index]}.png")
                c.Close()
                
            text_limits.close()
            


# Get Extrapolation Factor
def getExtrapFactor(pdftype, categ, bdtcut):
    file_path_1 = 'Slopes_' + categ + '_' + pdftype + '.txt'
    
    closest_cut = None       # To store the closest cut
    closest_ext_fact = None  # To store the extrapolation factor corresponding to the closest cut
    min_diff = float('inf')  # Initialize with a large number
    
    with open(file_path_1, 'r') as file1:
        # Read lines from the file
        for line1 in file1:
            # Split the line into components
            parts_1 = line1.split()
            
            cut = float(parts_1[1])
            ext_fact = float(parts_1[10])
            n_sideband = round(float(parts_1[5]))
            print(' round  ', float(parts_1[5]), parts_1)
            # Calculate the difference between the current cut and bdtcut
            diff = abs(cut - bdtcut)
            
            # Update the closest cut and its extrapolation factor if this cut is closer
            #print('diff ', diff, ' min_diff   ', min_diff, '  n_sideband   ', n_sideband)
            if diff < min_diff and n_sideband > 0.5:
                min_diff = diff
                closest_cut = cut
                closest_ext_fact = ext_fact
    
    # Return the closest cut and its extrapolation factor
    print('------>  ', closest_ext_fact)
    return closest_ext_fact







def MakeAndSaveExpFactors(datafile,categ,bdt_points):
        
        
        fit_range_lo = 1.6
        fit_range_hi = 2.0
        
        signal_range_lo = 1.74
        signal_range_hi = 1.81
        
        MiniTreeFile = ROOT.TFile.Open(datafile)
        MiniTreeFile.cd()
        
        treeName=''
        signalnorm = 1.0
        cat_label = ""
        if(categ == 'taue'):
                treeName  = 'ztau3mutaue'
                signalnorm = 0.00000856928
                cat_label = r"$\tau_{e}$"
        if(categ == 'taumu'):
                treeName = 'ztau3mutaumu'
                signalnorm = 0.00000822810
                cat_label = r"$\tau_{\mu}$"
        if(categ == 'tauhA'):
                treeName = 'ztau3mutauh_A'
                signalnorm = 0.00000815958
                cat_label = r"$\tau_{h,1-prong}$"
        if(categ == 'tauhB'):
                treeName = 'ztau3mutauh_B'
                signalnorm = 0.00000815958
                cat_label = r"$\tau_{h,3-prong}$"
        if(categ == 'all'):
                treeName   = 'ztautau'
                signalnorm = 0.00000824176
                cat_label = "Inclusive"
        
        tree = MiniTreeFile.Get(treeName)
        
        tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
        bdt_cv               = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
        dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 2)
        dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 2)
        event_weight         = ROOT.RooRealVar('weight'                     , 'event_weight'        ,  0,  5)  # this weight includes also the scale  mc signal scale
        category             = ROOT.RooRealVar('category'                   , 'category'            ,  0,  5)
        isMC                 = ROOT.RooRealVar('isMC'                       , 'isMC'                ,  0,  1000000)
        scale                = ROOT.RooRealVar('scale'                      , 'scale'               ,  signalnorm)  
        
        
        
        variables = ROOT.RooArgSet()
        variables.add(tripletMass)
        variables.add(bdt_cv)
        variables.add(dimu_OS1)
        variables.add(dimu_OS2)
        variables.add(event_weight)
        variables.add(category)
        variables.add(isMC)
        
        phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)&"
        omegavetoes="fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
        
        
        for bdt_cut in bdt_points:
        
                
                # For fitting BDT Output in Data
                
                BDT_Score_Min=-0.3
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector',' bdt_cv > ' + str(bdt_cut)+' & '+ phivetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fulldata = RooDataSet('fulldata', 'fulldata', tree,  variables, BlindDataSelector)
                
                bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0);
                
                
                tripletMass.setRange('left' , fit_range_lo    , signal_range_lo)
                tripletMass.setRange('right', signal_range_hi , fit_range_hi)
                tripletMass.setRange('full' , fit_range_lo    , fit_range_hi)
                
                tripletMass.setRange("SB1",fit_range_lo,1.75)
                tripletMass.setRange("SB2",1.80,fit_range_hi)
                tripletMass.setRange("fullRange",fit_range_lo,fit_range_hi)
                tripletMass.setRange("SIG",signal_range_lo,signal_range_hi)
                
                nbkg = ROOT.RooRealVar('nbkg', 'nbkg', 1000, 0, 500000)
                slope = ROOT.RooRealVar('slope', 'slope', 1.0, -100, 100)
                expo = ROOT.RooExponential('bkg_expo', 'bkg_expo', tripletMass, slope)
                pdfmodel = ROOT.RooAddPdf('bkg_extended_expo', 'bkg_extended_expo', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
                results_pdf = pdfmodel.fitTo(fulldata, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save())
                
                SG_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "SIG").getVal()
                SB_integral = pdfmodel.createIntegral(ROOT.RooArgSet(tripletMass), ROOT.RooArgSet(tripletMass), "left,right").getVal()
                print('>>>>>>>>>>>>>>>>>>>>> ',fulldata.numEntries())
                with open("Slopes_%s_%s"%(categ,"unfixed_exp")+".txt", "a") as f:
                        f.write("Cut: %s nbkg: %s n_sideband: %s expected_bkg: %s SG/SB ratio: %s \n"%(bdt_cut,nbkg.getVal(),fulldata.numEntries(),nbkg.getVal()*SG_integral,SG_integral/SB_integral))







if __name__ == "__main__":
        
        # Enable batch mode
        ROOT.gROOT.SetBatch(True)
        
        #categories = ['tauhA']
        categories = ['taue','taumu','tauhA','tauhB']
        #categories = ['tauhA','tauhB','all']
        #categories = ['combined'] # Can only be run after the other 4 categories are read and copied
        
        datafile_bdt_shape = "../../../../Combine_Tree_ztau3mutau_orig_PostBDT.root"
        
        datafile_main = "../../../../Combine_Tree_ztau3mutau_PF_PostBDT.root"#The final PF cuts I would prefer to use
        
        #lumi = np.round(np.arange(100,4500,500), 1)
        #lumi = np.insert(lumi, 0 , 59.83)
        #lumi = np.append(lumi, 137.0) # only 2016-2018
        #lumi = np.append(lumi, 400.0)
        #lumi = np.append(lumi, 2000.0)
        #lumi = np.append(lumi, 3000.0)
        #lumi = np.append(lumi, 4500.0)
        #lumi = np.sort(lumi)
        #lumi = np.round(lumi,1)
        
        #lumi = np.round([59.83],1)
        
        lumi = np.round([  59.83,  97.7,  137,  400,  600, 1100, 1600, 2000, 2600, 3000, 3600, 4100, 4500],1)
        lumi = np.sort(lumi)
        
        #lumi = np.round([  59.83,  3000],1)
        #lumi = np.round([  3000.0],1)
        #lumi = np.round([  59.83],1)
        #lumi = np.sort(lumi)
        
        cmd1 = 'mkdir lumi_limit_scans;'
        os.system(cmd1)
        
        #bdt_points = np.round(np.arange(0.2, 0.8 + 0.04, 0.04), 2)
        
        #bdt_points = np.round(np.arange(-0.1, 0.82, 0.04), 2)
        
        bdt_points = np.round(np.arange(-0.24, 0.82, 0.04), 2)
        
        Cat_No = len(categories)
        
        #To create datacards
        WhetherFitBDTandMakeCards = False
        
        for cat in range(Cat_No):
                categ = categories[cat]
                
                datafile_for_norm = datafile_main
                datafile_for_shape = datafile_main
                analyzed_lumi = 1.0
                if(categ == 'taue' and WhetherFitBDTandMakeCards):
                        analyzed_lumi = 59.83
                        datafile_for_norm = datafile_main
                        datafile_for_shape = datafile_bdt_shape
                if(categ == 'taumu'):
                        analyzed_lumi = 59.83
                        datafile_for_norm = datafile_main
                        datafile_for_shape = datafile_main
                if(categ == 'tauhA'):
                        analyzed_lumi = 59.83
                        datafile_for_norm = datafile_main
                        datafile_for_shape = datafile_bdt_shape
                if(categ == 'tauhB'):
                        analyzed_lumi = 59.83
                        datafile_for_norm = datafile_main
                        datafile_for_shape = datafile_bdt_shape
                if(categ == 'all'):
                        analyzed_lumi = 59.83
                        datafile_for_norm = datafile_main
                        datafile_for_shape = datafile_main
                if(categ == 'combined'):
                        analyzed_lumi = 59.83
                        datafile_for_norm = datafile_main
                        datafile_for_shape = datafile_main
                        
                if(WhetherFitBDTandMakeCards and (not categ == 'combined')):
                        open("Slopes_%s_%s"%(categ,"unfixed_exp")+".txt", 'w').close()
                        MakeAndSaveExpFactors(datafile_for_shape,categ,bdt_points)
                        BDTFit_Cat = makeCards()
                        BDTFit_Cat.FitBDT(datafile_for_norm,datafile_for_shape,categ)
                        BDTFit_Cat.MakeLumiScanCards(lumi,categ,analyzed_lumi)
                        
                if(WhetherFitBDTandMakeCards and categ == 'combined'):
                        BDTFit_Cat = makeCards()
                        BDTFit_Cat.CombineSubcategories(categ)
                
                
        #executeDataCards_onCondor(lumi,categories,False,bdt_points)
        #ReadAndCopyMinimumBDTCard(lumi,categories,False,bdt_points)
        
        #CalculateUL_fromDatacard(lumi,categories,False,bdt_points)
        ReadAndCopyMinimumBDTCard_usingUL(lumi,categories,True,bdt_points)
        
        #executeDataCards_onCondor(lumi,categories,True,bdt_points)
        #ReadAndCopyMinimumBDTCard(lumi,categories,True,bdt_points)
        
        
        