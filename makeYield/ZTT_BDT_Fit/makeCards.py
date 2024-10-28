import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
from ROOT import RooRealVar, RooFormulaVar, RooExponential, RooDataHist, RooArgList, RooAddPdf, RooFit, RooDataSet, RooGenericPdf, RooBifurGauss
#import CMS_lumi, tdrstyle
import subprocess # to execute shell command
import argparse
import numpy as np
#import tdrstyle
#from CMSStyle import CMS_lumi
import os

class makeCards:
        
        def __init__(self):
                
                self.bdt_cv = None
                
                self.bgausmeanMC = None
                self.bgaussigmaMC_a = None
                self.bgaussigmaMC_b = None
                self.bgaus_distMC = None
                self.BDTNorm_MC = None
                self.BDT_distribution_MC = None
                
                self.a = None
                self.b = None
                self.c = None
                self.quadratic = None
                self.expModel = None
                self.BDTNorm = None
                self.BDT_distribution_MC = None
        
        
        
        def FitBDT(self,datafile,categ):
                
                fit_range_lo = 1.6
                fit_range_hi = 2.0
                
                signal_range_lo = 1.74
                signal_range_hi = 1.81
                
                MiniTreeFile = ROOT.TFile.Open(datafile)
                MiniTreeFile.cd()
                
                treeName=''
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
                
                tree = MiniTreeFile.Get(treeName)
                
                tripletMass          = ROOT.RooRealVar('tripletMass'                , '3#mu mass'           , fit_range_lo, fit_range_hi, 'GeV')
                self.bdt_cv          = ROOT.RooRealVar('bdt_cv'                     , 'bdt_cv'              , -1 , 1)
                dimu_OS1             = ROOT.RooRealVar('dimu_OS1'                   , 'dimu_OS1'            ,  0 , 2)
                dimu_OS2             = ROOT.RooRealVar('dimu_OS2'                   , 'dimu_OS2'            ,  0 , 2)
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
                variables.add(scale)
                
                phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)"
                omegavetoes="&fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
                
                
                # For fitting BDT Output in Data
                
                BDT_Score_Min=-0.3
                
                BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector', phivetoes+omegavetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fulldata = RooDataSet('data', 'data', tree,  variables, BlindDataSelector)
                
                self.bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0);
                
                self.a = RooRealVar("a", "a", 1.0, -10.0, 10.0)
                self.b = RooRealVar("b", "b", 1.0, -10.0, 10.0)
                self.c = RooRealVar("c", "c", 1.0, -100.0, 10.0)
                
                #quadratic = RooFormulaVar("quadratic", "a + b*self.bdt_cv + c*self.bdt_cv*self.bdt_cv", RooArgList(a, b, c, self.bdt_cv))
                #expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(quadratic)) #Exponential of the quadratic polynomial
                
                self.quadratic = RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv", RooArgList(self.a, self.b, self.c, self.bdt_cv))
                self.expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(self.quadratic)) #Exponential of the quadratic polynomial
                
                self.BDTNorm = RooRealVar("BDTNorm", "BDTNorm", 500.0, 0.1, 50000)
                self.BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(self.expModel), RooArgList(self.BDTNorm))
                #BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(quadratic), RooArgList(BDTNorm))
                
                results_pdf = self.BDT_distribution.fitTo(fulldata, RooFit.Range('BDT_Fit_Range'), RooFit.Save())
                results_pdf.Print()
        
                
                
                
                # For fitting BDT Output in Signal
                
                BDT_Score_Min=-0.4
                
                MCSelector = RooFormulaVar('MCSelector', 'MCSelector', phivetoes+omegavetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
                
                fullmc = RooDataSet('mc', 'mc', tree, variables, MCSelector,'scale')
                
                self.bdt_cv.setRange("BDT_MC_Fit_Range", -1.0, 1.0);
                
                self.bgausmeanMC = RooRealVar("bgausmeanMC", "bgausmeanMC", 0.5, 0.0, 0.9)
                self.bgaussigmaMC_a = RooRealVar("bgaussigmaMC_a", "bgaussigmaMC_a", 0.2, 0.000001, 1.0)
                self.bgaussigmaMC_b = RooRealVar("bgaussigmaMC_b", "bgaussigmaMC_b", 0.2, 0.000001, 1.0)
                
                self.bgaus_distMC = RooBifurGauss("bgaus_distMC", "bgaus dist MC", self.bdt_cv, self.bgausmeanMC, self.bgaussigmaMC_a, self.bgaussigmaMC_b)
                
                self.BDTNorm_MC = RooRealVar("BDTNorm_MC", "BDTNorm_MC", 500.0, 0.1, 50000)
                self.BDT_distribution_MC = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(self.bgaus_distMC), RooArgList(self.BDTNorm_MC))
                
                results_mcpdf = self.BDT_distribution_MC.fitTo(fullmc, RooFit.Range('BDT_MC_Fit_Range'), RooFit.Save())
                results_mcpdf.Print()
                
                
                # Plot BDT for data and MC
                
                frame1 = self.bdt_cv.frame()
                frame1.SetTitle('')
                
                frame2 = self.bdt_cv.frame()
                frame2.SetTitle('')
                
                nbins = 100
                
                fullmc.plotOn(frame1, 
                              ROOT.RooFit.Binning(nbins), 
                              ROOT.RooFit.XErrorSize(0), 
                              ROOT.RooFit.LineWidth(2),
                              ROOT.RooFit.MarkerStyle(6),
                              ROOT.RooFit.MarkerColor(ROOT.kRed ),
                              ROOT.RooFit.MarkerSize(0.75),
                              ROOT.RooFit.FillColor(ROOT.kCyan  + 2)
                )
                self.BDT_distribution_MC.plotOn(frame1, ROOT.RooFit.LineColor(ROOT.kRed ))
                
                fulldata.plotOn(frame2, 
                        ROOT.RooFit.Binning(nbins),
                        ROOT.RooFit.MarkerStyle(20),
                        ROOT.RooFit.MarkerColor(ROOT.kBlack), 
                        ROOT.RooFit.MarkerSize(0.75))
                        
                #BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(fulldata.sumEntries("1", "BDT_Fit_Range"), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('BDT_Fit_Range') )
                #BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(BDTNorm.getVal()*BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "BDT_Fit_Range").getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('BDT_Fit_Range') )
                self.BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) )
                
                frame1.Draw()
                ROOT.gPad.SaveAs('bdt_fit_mc_'+categ+'.png')
                ROOT.gPad.Clear()
                
                frame2.Draw()
                ROOT.gPad.SaveAs('bdt_fit_bkg_'+categ+'.png')
                
                
                return frame1, frame2, fullmc, fulldata
        
        
        #To create datacards
        def MakeLumiScanCards(self,lumi,categ):
                
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
                
                num_points = 20
                bdt_points = np.round(np.linspace(0.0,1.0,num_points), 2)
                
                
                for lu_no in range(len(lumi)):
                    lu = lumi[lu_no]
                    print(lu)
                    
                    command_recreate_lumidir = "rm -r lumi_limit_scans/{0}/lumi_{1}; mkdir lumi_limit_scans/{0}/lumi_{1}".format(categ, int(lu))
                    os.system(command_recreate_lumidir)
                    
                    for point in bdt_points:  # For loop for bdt cuts in range [X_min;X_max]
                        
                        self.bdt_cv.setRange("Integral_Range", point, 100.0)
                        
                        print "test 1"
                        
                        sig_est = signalnorm * self.BDTNorm_MC.getVal() * (self.BDT_distribution_MC.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range").getVal() )
                        bkg_est = exp_fact * self.BDTNorm.getVal() * (self.BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range").getVal() )
                        sb_est = exp_fact * self.BDTNorm.getVal() * (self.BDT_distribution.createIntegral(ROOT.RooArgSet(self.bdt_cv), ROOT.RooArgSet(self.bdt_cv), "Integral_Range").getVal() )
                        
                        print "test 2"
                        
                        command_mod_card = "python card_modifiers/Card_Mod.py --categ " + str(categ) + " --sig_exp " + str(sig_est) + " --bkg_exp " + str(bkg_est) + " --sb_exp " + str(sb_est) 
                        os.system(command_mod_card)
                        
                        command_copy_dc = "cp modified_dc_{0}.txt lumi_limit_scans/{0}/lumi_{1}/dc_{2}.txt".format(categ, int(lu), str(lu))
                        os.system(command_copy_dc)
        

if __name__ == "__main__":
    
        # Enable batch mode
        ROOT.gROOT.SetBatch(True)
        
        datafile = "Combine_Tree_ztau3mutau.root"
        category = 'taumu'
        
        
        BDTFit_Cat = makeCards()
        BDTFit_Cat.FitBDT(datafile,category)
        #print "is valid BDT_distribution_MC: ",isinstance(BDTFit_Cat.BDT_distribution_MC, ROOT.RooAddPdf), BDTFit_Cat.BDT_distribution_MC
        
        lumi = np.round(np.arange(100,3050,100), 0)
        lumi = np.insert(lumi, 0 , 59.8)
        
        cmd1 = 'mkdir lumi_limit_scans;'
        os.system(cmd1)
        
        cmd2 = 'mkdir datacards_modified;'
        os.system(cmd2)
        
        BDTFit_Cat.MakeLumiScanCards(lumi,category)
        

        
        
        