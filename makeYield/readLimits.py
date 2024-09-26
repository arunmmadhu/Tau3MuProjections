#!/usr/bin/env python

import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
import CMS_lumi, tdrstyle
import subprocess # to execute shell command
import argparse
import numpy as np
import tdrstyle
from CMSStyle import CMS_lumi
import os


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--category', required=True               , help='W, HF, ZTT',  dest='category'          , default='W') 
parser.add_argument('-r', '--runType' , default='plot'              , help='Option run or plot', dest='runType'   )
parser.add_argument('-m', '--Method' , default='A'              , help='Option A=Asymptotic, H=HybridNew', dest='Method' )
args = parser.parse_args()



ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
CMS_lumi.cmsText = "CMS, work in progress"
CMS_lumi.extraText = ""
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


#lumi = [97.7,129.0,150,377.0,500.0,700.0,1000.0, 1200.0, 1500.0, 1700.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0 ]
#$lumi = [59.0,97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]
#lumi = [59.0,97.7, 129.0]
lumi = np.round(np.arange(500,2000,20), 0)

def executeDataCards_onCondor(lumi,categories,Whether_Hybrid):
        
        Cat_No = len(categories)
        
        for cat in range(Cat_No):
        
                for lu_no in range(len(lumi)):
                        
                    lu = lumi[lu_no]
                    print(lu)
                        
                    if(Whether_Hybrid):
                            
                            print("Median")
                            command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.5 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+categories[cat],categories[cat]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                            os.system(command_run)
                            
                            print("Running Sigma -2")
                            command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.025 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+categories[cat],categories[cat]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                            os.system(command_run)
                            print("Running Sigma -1")
                            command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.16 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+categories[cat],categories[cat]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                            os.system(command_run)
                            print("Running Sigma +1")
                            command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.84 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+categories[cat],categories[cat]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                            os.system(command_run)
                            print("Running Sigma +2")
                            command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.975 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+categories[cat],categories[cat]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                            os.system(command_run)
                            
                    else:
                            
                            command_run = "combineTool.py -M AsymptoticLimits  -n %s -d %s --cl 0.90  --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+categories[cat],categories[cat]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                            os.system(command_run)
            

# GET limits from root file
def getLimits(file_name):
 
    file = TFile(file_name)
    tree = file.Get("limit")
 
    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]
 
    return limits[:6]

# PLOT upper limits
def plotUpperLimits(lumi,categories,Whether_Hybrid):
 
    N = len(lumi)
    Cat_No = len(categories)
    
    label=[None] * Cat_No
    
    for cat in range(Cat_No):
            if categories[cat]=='HF':
                    label[cat] = 'Heavy Flavor'
            if categories[cat]=='ZTT':
                    label[cat] = 'ZTT'
            if categories[cat]=='W':
                    label[cat] = 'W'
    
    yellow=[None] * Cat_No
    green=[None] * Cat_No
    median=[None] * Cat_No
    
    for cat in range(Cat_No):
    
            yellow[cat] = TGraph(2*N)   
            green[cat] = TGraph(2*N)  
            median[cat] = TGraph(N)  
            
            text_limits=open("TextLimits.txt","w")
            
            if(Whether_Hybrid):
                    for i in range(N):
                        file_name1 = "higgsCombine"+str(lumi[i])+categories[cat]+".HybridNew.mH120.123456.quant0.025.root"
                        limit1 = getLimits(file_name1)
                        file_name2 = "higgsCombine"+str(lumi[i])+categories[cat]+".HybridNew.mH120.123456.quant0.160.root"
                        limit2 = getLimits(file_name2)
                        file_name3 = "higgsCombine"+str(lumi[i])+categories[cat]+".HybridNew.mH120.123456.quant0.500.root"
                        limit3 = getLimits(file_name3)
                        file_name4 = "higgsCombine"+str(lumi[i])+categories[cat]+".HybridNew.mH120.123456.quant0.840.root"
                        limit4 = getLimits(file_name4)
                        file_name5 = "higgsCombine"+str(lumi[i])+categories[cat]+".HybridNew.mH120.123456.quant0.975.root"
                        limit5 = getLimits(file_name5)
                        
                        yellow.SetPoint( 2*N-1-i, lumi[i], limit1[2]) # - 2 sigma
                        green[cat].SetPoint(  2*N-1-i, lumi[i], limit2[2]) # - 1 sigma
                        median[cat].SetPoint(    i,    lumi[i], limit3[2]) #    median
                        green[cat].SetPoint(     i,    lumi[i], limit4[2]) # + 1 sigma
                        yellow[cat].SetPoint(    i,    lumi[i], limit5[2]) # + 2 sigma
                        
                        text_limits.write("bdt %.2f     median[cat] exp %.2f\n"%(lumi[i],limit3[2]))
            
            else:
                    for i in range(N):
                        file_name1 = "higgsCombine"+str(lumi[i])+categories[cat]+".AsymptoticLimits.mH120.root"
                        limit1 = getLimits(file_name1)
                        
                        yellow[cat].SetPoint( 2*N-1-i, lumi[i], limit1[0]) # - 2 sigma
                        green[cat].SetPoint(  2*N-1-i, lumi[i], limit1[1]) # - 1 sigma
                        median[cat].SetPoint(    i,    lumi[i], limit1[2]) #    median
                        green[cat].SetPoint(     i,    lumi[i], limit1[3]) # + 1 sigma
                        yellow[cat].SetPoint(    i,    lumi[i], limit1[4]) # + 2 sigma
                        
                        text_limits.write("bdt %.2f     median[cat] exp %.2f\n"%(lumi[i],limit1[2]))

    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetGrid()
    c.SetFrameLineWidth(2);
    c.SetTickx();
    c.SetTicky();
    c.SetLogy();
    c.cd()

    frame = c.DrawFrame(1.4,0.001, 4.1, 1.2)
    frame.GetYaxis().CenterTitle()

    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(False)
    frame.GetYaxis().SetTitle("Expected Limit (10^{-7})")
    frame.GetXaxis().SetTitle("Integrated Luminosity fb^{-1}")




    #frame.SetMinimum(median[cat].GetHistogram().GetMinimum()*0.8)
    #frame.SetMaximum(median[cat].GetHistogram().GetMinimum()+(5.0/3.0)*(median[cat].GetHistogram().GetMaximum()-median[cat].GetHistogram().GetMinimum()))
    
    frame.SetMinimum(0.01)
    frame.SetMaximum(5.0)
    
    frame.GetXaxis().SetLimits(min(lumi) ,max(lumi)*1.2)


    for cat in range(Cat_No):
            yellow[cat].SetFillColor(ROOT.kOrange)
            yellow[cat].SetLineColor(ROOT.kOrange)
            yellow[cat].SetFillStyle(1001)
            yellow[cat].Draw('F')
         
            green[cat].SetFillColor(ROOT.kGreen+1)
            green[cat].SetLineColor(ROOT.kGreen+1)
            green[cat].SetFillStyle(1001)
            green[cat].Draw('Fsame')
         
            median[cat].SetLineColor(1)
            median[cat].SetLineWidth(2)
            median[cat].SetLineStyle(2)
            median[cat].Draw('Lsame')
            median[cat].Draw()
 
#    CMS_lumi.CMS_lumi(c,13,11)
    ROOT.gPad.SetTicks(1,1)
    CMS_lumi(ROOT.gPad, 5, 0)
    ROOT.gPad.Update()
    frame.Draw('sameaxis')
 
    x1 = 0.15
    x2 = x1 + 0.24
    y2 = 0.86
    y1 = 0.70
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    for cat in range(Cat_No):
            legend.AddEntry(median[cat], label[cat],'L')
    
#    legend.AddEntry(median[cat], "HybridNew CL_{s} expected upper limit",'L')
#    legend.AddEntry(green[cat], "#pm 1 std. deviation",'f')
#    legend.AddEntry(yellow[cat],"#pm 2 std. deviation",'f')


    legend.Draw()
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    Text = 'Projected Expected UL'
    

    latex.SetTextAlign(1)
    latex.DrawLatex(0.15, 0.85, Text)
    latex.Draw('same') 
    print " "
    if len(categories) < 2:
            c.SaveAs("Limit_scan_"+categories[0]+".png")
    else:
            c.SaveAs("Limit_scan.png")
    c.Close()
    

        


# MAIN
def main():


        
#    categories = ['ZTT']
#    categories = ['W']
#    categories = ['HF']
    
    Whether_Hybrid=False
    if(args.Method == 'H'):
            Whether_Hybrid=True
    
    categories=[args.category]
    
    if(args.runType == 'run'):
            executeDataCards_onCondor(lumi,categories,Whether_Hybrid)
    if(args.runType == 'plot'):
            plotUpperLimits(lumi,categories,Whether_Hybrid)

if __name__ == '__main__':
    main()
