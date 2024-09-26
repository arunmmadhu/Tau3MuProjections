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

        
def FitBDT(datafile,categ):
        
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
        variables.add(scale)
        
        phivetoes="(fabs(dimu_OS1 - 1.020)>0.020)&(fabs(dimu_OS2 - 1.020)>0.020)"
        omegavetoes="&fabs(dimu_OS1 - 0.782)>0.020&fabs(dimu_OS2 - 0.782)>0.020&"
        
        
        # For fitting BDT Output in Data
        
        BDT_Score_Min=-0.3
        
        BlindDataSelector = RooFormulaVar('DataSelector', 'DataSelector', phivetoes+omegavetoes+' isMC == 0 & (tripletMass<=%s || tripletMass>=%s) & (tripletMass>=%s & tripletMass<=%s) ' %(signal_range_lo,signal_range_hi,fit_range_lo,fit_range_hi) , RooArgList(variables))
        
        fulldata = RooDataSet('data', 'data', tree,  variables, BlindDataSelector)
        
        bdt_cv.setRange("BDT_Fit_Range", BDT_Score_Min, 1.0);
        
        a = RooRealVar("a", "a", 1.0, -10.0, 10.0)
        b = RooRealVar("b", "b", 1.0, -10.0, 10.0)
        c = RooRealVar("c", "c", 1.0, -100.0, 10.0)
        
        #quadratic = RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv", RooArgList(a, b, c, bdt_cv))
        #expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(quadratic)) #Exponential of the quadratic polynomial
        
        quadratic = RooFormulaVar("quadratic", "a + b*bdt_cv + c*bdt_cv*bdt_cv", RooArgList(a, b, c, bdt_cv))
        expModel = RooGenericPdf("expModel", "exp(quadratic)", RooArgList(quadratic)) #Exponential of the quadratic polynomial
        
        BDTNorm = RooRealVar("BDTNorm", "BDTNorm", 500.0, 0.1, 50000)
        BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(expModel), RooArgList(BDTNorm))
        #BDT_distribution = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(quadratic), RooArgList(BDTNorm))
        
        results_pdf = BDT_distribution.fitTo(fulldata, RooFit.Range('BDT_Fit_Range'), RooFit.Save())
        results_pdf.Print()

        
        
        
        # For fitting BDT Output in Signal
        
        BDT_Score_Min=-0.4
        
        MCSelector = RooFormulaVar('MCSelector', 'MCSelector', phivetoes+omegavetoes+' isMC !=0 & (isMC == 211 | isMC == 210231 | isMC == 210232 | isMC == 210233 ) & (tripletMass>=%s & tripletMass<=%s) ' %(fit_range_lo,fit_range_hi) , RooArgList(variables))
        
        fullmc = RooDataSet('mc', 'mc', tree, variables, MCSelector,'scale')
        
        #BDTOutput_x = RooRealVar("BDTOutput_x", "BDT Output", -1.0, 1.0)
        bdt_cv.setRange("BDT_MC_Fit_Range", -1.0, 1.0);
        
        bgausmeanMC = RooRealVar("bgausmeanMC", "bgausmeanMC", 0.5, 0.0, 0.9)
        bgaussigmaMC_a = RooRealVar("bgaussigmaMC_a", "bgaussigmaMC_a", 0.2, 0.000001, 1.0)
        bgaussigmaMC_b = RooRealVar("bgaussigmaMC_b", "bgaussigmaMC_b", 0.2, 0.000001, 1.0)
        
        bgaus_distMC = RooBifurGauss("bgaus_distMC", "bgaus dist MC", bdt_cv, bgausmeanMC, bgaussigmaMC_a, bgaussigmaMC_b)
        
        BDTNorm_MC = RooRealVar("BDTNorm_MC", "BDTNorm_MC", 500.0, 0.1, 50000)
        BDT_distribution_MC = RooAddPdf("BDT_distribution", "BDT_distribution",RooArgList(bgaus_distMC), RooArgList(BDTNorm_MC))
        
        results_mcpdf = BDT_distribution_MC.fitTo(fullmc, RooFit.Range('BDT_MC_Fit_Range'), RooFit.Save())
        results_mcpdf.Print()
        
        
        # Plot BDT for data and MC
        
        frame1 = bdt_cv.frame()
        frame1.SetTitle('')
        
        frame2 = bdt_cv.frame()
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
        BDT_distribution_MC.plotOn(frame1, ROOT.RooFit.LineColor(ROOT.kRed ))
        
        fulldata.plotOn(frame2, 
                ROOT.RooFit.Binning(nbins),
                ROOT.RooFit.MarkerStyle(20),
                ROOT.RooFit.MarkerColor(ROOT.kBlack), 
                ROOT.RooFit.MarkerSize(0.75))
                
        #BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(fulldata.sumEntries("1", "BDT_Fit_Range"), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('BDT_Fit_Range') )
        #BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) , ROOT.RooFit.Normalization(BDTNorm.getVal()*BDT_distribution.createIntegral(ROOT.RooArgSet(bdt_cv), ROOT.RooArgSet(bdt_cv), "BDT_Fit_Range").getVal(), ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange('BDT_Fit_Range') )
        BDT_distribution.plotOn(frame2,  ROOT.RooFit.LineColor(ROOT.kBlue) )
        
        frame1.Draw()
        ROOT.gPad.SaveAs('bdt_fit_mc_'+categ+'.png')
        ROOT.gPad.Clear()
        
        frame2.Draw()
        ROOT.gPad.SaveAs('bdt_fit_bkg_'+categ+'.png')
        
        return BDT_distribution_MC, BDT_distribution, BDTNorm_MC, BDTNorm, frame1, frame2, fullmc, fulldata


if __name__ == "__main__":
    
        # Enable batch mode
        ROOT.gROOT.SetBatch(True)
        
        # specify the luminosities here
        lumi_values = [59.0]
        lumi_size = len(lumi_values)
        
        datafile = "Combine_Tree_ztau3mutau.root"
        
        signalnorm = 0.00000824176
        
        BDT_distribution_MC, BDT_distribution, BDTNorm_MC, BDTNorm, frame1, frame2, fullmc, fulldata = FitBDT(datafile,'taumu')
        
        #FitBDT(datafile,'taue',signalnorm)
        
        
        
        
        
        
        
        