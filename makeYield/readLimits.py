#!/usr/bin/env python3

import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TPaveLabel, TPaveText, TLatex
import subprocess # to execute shell command
import argparse
import numpy as np
import tdrstyle
#from CMSStyle import CMS_lumi
import CMSStyle
import os
import ctypes
import math


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--category', required=True               , help='W, HF, ZTT',  dest='category'          , default='W') 
parser.add_argument('-r', '--runType' , default='plot'              , help='Option run or plot', dest='runType'   )
parser.add_argument('-m', '--Method' , default='A'                  , help='Option A=Asymptotic, H=HybridNew', dest='Method' )
args = parser.parse_args()



ROOT.gROOT.SetBatch(ROOT.kTRUE)

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


#lumi = [97.7,129.0,150,377.0,500.0,700.0,1000.0, 1200.0, 1500.0, 1700.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0 ]
#lumi = [59.0,97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]
#lumi = [59.0,97.7, 129.0]

analyzed_lumi = 1.0
#lumi_W = np.round(np.arange(100,3050,100), 0)
#lumi_W = np.insert(lumi_W, 0 , 59.8)
#lumi_W = np.append(lumi_W, 2000)
#lumi_W = np.append(lumi_W, 3000)
#lumi_W = np.append(lumi_W, 4500)
#lumi_W = np.sort(lumi_W)
lumi_W = lumi = np.round([  59.83,  97.7,  137,  400,  600, 1100, 1600, 2000, 2600, 3000, 3600, 4100, 4500],1)
lumi_W = np.sort(lumi_W)
analyzed_lumi_W = 59.8

#lumi_HF = np.round(np.arange(100,3050,100), 0)
#lumi_HF = np.insert(lumi_HF, 0 , 97.7)
#lumi_HF = np.append(lumi_HF, 4500)
#lumi_HF = np.sort(lumi_HF)
lumi_HF = np.round([ 97.7,  137,  400,  600, 1100, 1600, 2000, 2600, 3000, 3600, 4100, 4500],1)
lumi_HF = np.sort(lumi_HF)
analyzed_lumi_HF = 97.7

#lumi_ZTT = np.round(np.arange(100,4500,500), 0)
#lumi_ZTT = np.insert(lumi_ZTT, 0 , 59.8)
#lumi_ZTT = np.append(lumi_ZTT, 2000)
#lumi_ZTT = np.append(lumi_ZTT, 3000)
#lumi_ZTT = np.append(lumi_ZTT, 4500)
#lumi_ZTT = np.sort(lumi_ZTT)
lumi_ZTT = lumi = np.round([  59.83,  97.7,  137,  400,  600, 1100, 1600, 2000, 2600, 3000, 3600, 4100, 4500],1)
lumi_ZTT = np.sort(lumi_ZTT)
analyzed_lumi_ZTT = 59.8

lumi_Combo = np.round(np.arange(100,3050,100), 0)
analyzed_lumi_Combo = 59.8




def executeDataCards_onCondor(lumi_W,lumi_HF,lumi_ZTT,categories,Whether_Hybrid):
        
        lumi = lumi_W
        Cat_No = len(categories)
        
        for cat in range(Cat_No):
        
                if categories[cat]=='HF':
                    lumi = lumi_HF
                if categories[cat]=='ZTT':
                    lumi = lumi_ZTT
                if categories[cat]=='W':
                    lumi = lumi_W
                
                for lu_no in range(len(lumi)):
                        
                    lu = lumi[lu_no]
                    print('Luminosity: ', lu)
                    
                    subcat = []
                    
                    if(categories[cat]=="W"):
                            #subcat = ["W"]
                            subcat = ['CatA','CatB','CatC','combined']
                    if(categories[cat]=="HF"):
                            subcat = ["HF"]
                    if(categories[cat]=="ZTT"):
                            #subcat = ['taue','taumu','tauhA','tauhB','all','combined']
                            subcat = ['taue','taumu','tauhA','tauhB','combined']
                            #subcat = ['combined']
                    if(categories[cat]=="Combo"):
                            subcat = ["Combo"]
                    
                    Cat_No_sub = len(subcat)
                    
                    for cat_sub in range(Cat_No_sub):
                    
                            if(Whether_Hybrid):
                                    
                                    print("Median")
                                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.5 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub],categories[cat]+"/"+subcat[cat_sub]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub])
                                    print("Run:   ", command_run)
                                    os.system(command_run)
                                    
                                    print("Running Sigma -2")
                                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.025 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub],categories[cat]+"/"+subcat[cat_sub]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                                    os.system(command_run)
                                    print("Running Sigma -1")
                                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.16 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub],categories[cat]+"/"+subcat[cat_sub]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                                    os.system(command_run)
                                    print("Running Sigma +1")
                                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.84 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub],categories[cat]+"/"+subcat[cat_sub]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                                    os.system(command_run)
                                    print("Running Sigma +2")
                                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.975 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub],categories[cat]+"/"+subcat[cat_sub]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+categories[cat])
                                    os.system(command_run)
                                    
                            else:
                                    
                                    command_run = "combineTool.py -M AsymptoticLimits --run blind  --cl 0.90 -n %s -d %s  --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name Asymptotic%s " % (str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub],categories[cat]+"/"+subcat[cat_sub]+"/datacards_modified/dc_"+str(lu)+".txt",str(lu)+"_"+categories[cat]+"_"+subcat[cat_sub])
        
                                    print("Run:   ", command_run)
                                    os.system(command_run)
            

# GET limits from root file
def getLimits(file_name):
 
    file = TFile(file_name)
    tree = file.Get("limit")
 
    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print(">>>   %.2f" % limits[-1])
 
    return limits[:6]

# PLOT upper limits
def plotUpperLimits(lumi_W,lumi_HF,lumi_ZTT,analyzed_lumi_W,analyzed_lumi_HF,analyzed_lumi_ZTT,categories,Whether_Hybrid,WhetherIndividualPlots):
 
    #placeholder
    N = len(lumi_HF) #dummy value
    lumi = lumi_W #dummy value
    analyzed_lumi = analyzed_lumi_W #dummy value
    
    N1 = len(lumi_HF)
    N2 = len(lumi_ZTT)
    N3 = len(lumi_W)
    Cat_No = len(categories)
    
    # Whether to make individual plots
    WhetherIndividualPlots = True
    
    WhetherMultipleBroadCategories = False
    if (Cat_No>1):
            WhetherMultipleBroadCategories = True
    
    label=[None] * Cat_No
    
    subcat = [None] * Cat_No
    subcat_label = [None] * Cat_No
    Cat_No_sub = [None] * Cat_No
    
    for cat in range(Cat_No):
            if categories[cat]=='HF':
                    label[cat] = 'HF'
                    subcat[cat] = ["HF"]
                    subcat_label[cat] = ["HF"]
            if categories[cat]=='ZTT':
                    label[cat] = 'ZTT'
                    if(WhetherMultipleBroadCategories):
                            subcat[cat] = ['combined']
                            subcat_label[cat] = ['ZTT']
                    else:
                            #subcat[cat] = ['combined','taue','taumu','tauhA','tauhB','all']
                            #subcat_label[cat] = ["ZTT Combined", r"$\tau_{e}$", r"$\tau_{\mu}$", r"$\tau_{h,1-prong}$", r"$\tau_{h,3-prong}$", "Inclusive"]
                            subcat[cat] = ['combined','taue','taumu','tauhA','tauhB']
                            subcat_label[cat] = ["ZTT Combined", r"$\tau_{e}$", r"$\tau_{\mu}$", r"$\tau_{h,1-prong}$", r"$\tau_{h,3-prong}$"]
                            
                    #subcat = ['taue','taumu','tauhA','tauhB','all','combined']
                    #subcat_label = [r"$\tau_{e}$", r"$\tau_{\mu}$", r"$\tau_{h,1-prong}$", r"$\tau_{h,3-prong}$", "Inclusive", "Combined"]
                    #subcat = ['combined']
                    #subcat_label = ['combined']
            if categories[cat]=='W':
                    label[cat] = 'W'
                    #subcat[cat] = ["W"]
                    #subcat_label[cat] = ["W"]
                    if(WhetherMultipleBroadCategories):
                            subcat[cat] = ['combined']
                            subcat_label[cat] = ['W']
                    else:
                            subcat[cat] = ['combined','CatA','CatB','CatC']
                            subcat_label[cat] = ["W Combined","Category A", "Category B", "Category C"]
                    
            Cat_No_sub[cat] = len(subcat[cat])
    
    WhetherMultipleSmallCategories = False
    if (Cat_No_sub[0]>1):
            WhetherMultipleSmallCategories = True
    
    yellow = [[None] * Cat_No_sub[cat] for cat in range(Cat_No)]
    green = [[None] * Cat_No_sub[cat] for cat in range(Cat_No)]
    median = [[None] * Cat_No_sub[cat] for cat in range(Cat_No)]

    
    for cat in range(Cat_No):
    
            for cat_sub in range(Cat_No_sub[cat]):
                    
                    if categories[cat]=='HF':
                            N=N1
                            lumi = lumi_HF
                    if categories[cat]=='ZTT':
                            N=N2
                            lumi = lumi_ZTT
                    if categories[cat]=='W':
                            N=N3
                            lumi = lumi_W
                    
                    
                    yellow[cat][cat_sub] = TGraph(2*N)   
                    green[cat][cat_sub] = TGraph(2*N)  
                    median[cat][cat_sub] = TGraph(N)  
                    
                    text_limits=open("TextLimits.txt","w")
                    
                    
                    if(Whether_Hybrid):
                            for i in range(N):
                                file_name1 = "higgsCombine"+str(lumi[i])+"_"+categories[cat]+"_"+subcat[cat][cat_sub]+".HybridNew.mH120.123456.quant0.025.root"
                                limit1 = getLimits(file_name1)
                                file_name2 = "higgsCombine"+str(lumi[i])+"_"+categories[cat]+"_"+subcat[cat][cat_sub]+".HybridNew.mH120.123456.quant0.160.root"
                                limit2 = getLimits(file_name2)
                                file_name3 = "higgsCombine"+str(lumi[i])+"_"+categories[cat]+"_"+subcat[cat][cat_sub]+".HybridNew.mH120.123456.quant0.500.root"
                                limit3 = getLimits(file_name3)
                                file_name4 = "higgsCombine"+str(lumi[i])+"_"+categories[cat]+"_"+subcat[cat][cat_sub]+".HybridNew.mH120.123456.quant0.840.root"
                                limit4 = getLimits(file_name4)
                                file_name5 = "higgsCombine"+str(lumi[i])+"_"+categories[cat]+"_"+subcat[cat][cat_sub]+".HybridNew.mH120.123456.quant0.975.root"
                                limit5 = getLimits(file_name5)
                                
                                yellow[cat][cat_sub].SetPoint( 2*N-1-i, lumi[i], limit1[2]) # - 2 sigma
                                green[cat][cat_sub].SetPoint(  2*N-1-i, lumi[i], limit2[2]) # - 1 sigma
                                median[cat][cat_sub].SetPoint(    i,    lumi[i], limit3[2]) #    median
                                green[cat][cat_sub].SetPoint(     i,    lumi[i], limit4[2]) # + 1 sigma
                                yellow[cat][cat_sub].SetPoint(    i,    lumi[i], limit5[2]) # + 2 sigma
                                
                                text_limits.write("bdt %.2f     median[cat] exp %.2f\n"%(lumi[i],limit3[2]))
                    
                    else:
                            for i in range(N):
                                #file_name1 = "higgsCombine" + str(lumi[i]) + "_" +categories[cat]+ "_" + "_".join(subcat[cat][cat_sub]) + ".AsymptoticLimits.mH120.root"
                                #This subcat[cat][0] needs to be corrected. Works for WhetherMultipleBroadCategories.
                                file_name1 = "higgsCombine"+str(lumi[i])+"_"+categories[cat]+"_"+subcat[cat][cat_sub]+".AsymptoticLimits.mH120.root"
                                
                                print("filename: ",file_name1)
                                
                                limit1 = getLimits(file_name1)
                                
                                yellow[cat][cat_sub].SetPoint( 2*N-1-i, lumi[i], limit1[0]) # - 2 sigma
                                green[cat][cat_sub].SetPoint(  2*N-1-i, lumi[i], limit1[1]) # - 1 sigma
                                median[cat][cat_sub].SetPoint(    i,    lumi[i], limit1[2]) #    median
                                green[cat][cat_sub].SetPoint(     i,    lumi[i], limit1[3]) # + 1 sigma
                                yellow[cat][cat_sub].SetPoint(    i,    lumi[i], limit1[4]) # + 2 sigma
                                
                                text_limits.write("bdt %.2f     median[cat] exp %.2f\n"%(lumi[i],limit1[2]))

    
    
    #if((not WhetherMultipleBroadCategories) and WhetherMultipleSmallCategories and not WhetherIndividualPlots ):
    #when you're making a plot for one category together in a single image. This will plot ZTT or W or HF incl sub categories.
    if((not WhetherMultipleBroadCategories) and not WhetherIndividualPlots ):
            
            if categories[cat]=='HF':
                    lumi = lumi_HF
                    analyzed_lumi = analyzed_lumi_HF
            if categories[cat]=='ZTT':
                    lumi = lumi_ZTT
                    analyzed_lumi = analyzed_lumi_ZTT
            if categories[cat]=='W':
                    lumi = lumi_W
                    analyzed_lumi = analyzed_lumi_W
            
            c = TCanvas("c","c",100,100,W_ref,H_ref)
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
            
            #
            c.SetLogy();
            c.cd()
        
            frame = c.DrawFrame(1.4,0.001, 4.1, 1.2)
            frame.GetYaxis().CenterTitle()
        
            frame.GetYaxis().SetTitleSize(0.05)
            frame.GetXaxis().SetTitleSize(0.05)
            frame.GetXaxis().SetLabelSize(0.04)
            frame.GetYaxis().SetLabelSize(0.04)
            frame.GetYaxis().SetTitleOffset(1.1)
            frame.GetXaxis().SetNdivisions(508)
            frame.GetYaxis().CenterTitle(False)
            frame.GetYaxis().SetTitle("Expected Limit (10^{-7})")
            frame.GetXaxis().SetTitle("Integrated Luminosity fb^{-1}")
        
        
        
        
            #frame.SetMinimum(median[cat].GetHistogram().GetMinimum()*0.8)
            #frame.SetMaximum(median[cat].GetHistogram().GetMinimum()+(5.0/3.0)*(median[cat].GetHistogram().GetMaximum()-median[cat].GetHistogram().GetMinimum()))
            
            frame.SetMinimum(0.01)
            frame.SetMaximum(25.0)
            
            frame.GetXaxis().SetLimits(min(lumi) ,max(lumi)*1.2)
            
            c.SetLogx();
            c.cd()
            
            possible_colors = [1,2,4,6,800,880,416,9]
            
            #for cat in range(Cat_No):
            if categories[0]=='ZTT' or categories[0]=='W':
                    for cat_sub in range(Cat_No_sub[cat]):
                            if subcat[cat][cat_sub] == 'combined':
                                    yellow[cat][cat_sub].SetFillColor(ROOT.kOrange)
                                    yellow[cat][cat_sub].SetLineColor(ROOT.kOrange)
                                    yellow[cat][cat_sub].SetFillStyle(1001)
                                    yellow[cat][cat_sub].Draw('F')
                                 
                                    green[cat][cat_sub].SetFillColor(ROOT.kGreen+1)
                                    green[cat][cat_sub].SetLineColor(ROOT.kGreen+1)
                                    green[cat][cat_sub].SetFillStyle(1001)
                                    green[cat][cat_sub].Draw('Fsame')
                                 
                                    median[cat][cat_sub].SetLineColor(possible_colors[cat_sub])
                                    median[cat][cat_sub].SetLineWidth(2)
                                    median[cat][cat_sub].SetLineStyle(2)
                                    median[cat][cat_sub].Draw('Lsame')
                                    median[cat][cat_sub].Draw()
                            else:
                                    #yellow[cat][cat_sub].SetFillColor(ROOT.kOrange)
                                    #yellow[cat][cat_sub].SetLineColor(ROOT.kOrange)
                                    #yellow[cat][cat_sub].SetFillStyle(1001)
                                    #yellow[cat][cat_sub].Draw('F')
                                 
                                    #green[cat][cat_sub].SetFillColor(ROOT.kGreen+1)
                                    #green[cat][cat_sub].SetLineColor(ROOT.kGreen+1)
                                    #green[cat][cat_sub].SetFillStyle(1001)
                                    #green[cat][cat_sub].Draw('Fsame')
                                 
                                    median[cat][cat_sub].SetLineColor(possible_colors[cat_sub])
                                    median[cat][cat_sub].SetLineWidth(2)
                                    median[cat][cat_sub].SetLineStyle(2)
                                    median[cat][cat_sub].Draw('Lsame')
                                    median[cat][cat_sub].Draw()
            else:
                    for cat_sub in range(Cat_No_sub[cat]):
                            yellow[cat][cat_sub].SetFillColor(ROOT.kOrange)
                            yellow[cat][cat_sub].SetLineColor(ROOT.kOrange)
                            yellow[cat][cat_sub].SetFillStyle(1001)
                            yellow[cat][cat_sub].Draw('F')
                         
                            green[cat][cat_sub].SetFillColor(ROOT.kGreen+1)
                            green[cat][cat_sub].SetLineColor(ROOT.kGreen+1)
                            green[cat][cat_sub].SetFillStyle(1001)
                            green[cat][cat_sub].Draw('Fsame')
                         
                            median[cat][cat_sub].SetLineColor(possible_colors[cat_sub])
                            median[cat][cat_sub].SetLineWidth(2)
                            median[cat][cat_sub].SetLineStyle(2)
                            median[cat][cat_sub].Draw('Lsame')
                            median[cat][cat_sub].Draw()
         
        #    CMS_lumi.CMS_lumi(c,13,11)
            ROOT.gPad.SetTicks(1,1)
        #    CMS_lumi(ROOT.gPad, 5, 0)
            ROOT.gPad.Update()
            frame.Draw('sameaxis')
            
            
            if categories[cat] == 'HF':
                    x, y = ctypes.c_double(0), ctypes.c_double(0)  # Create c_double variables for x and y
                    median[0][0].GetPoint(0, x, y)  # Get the point at index i
                    # Convert ctypes values to Python floats
                    x_val = x.value
                    y_val = y.value
                    
                    # Create and draw the scaled 1/sqrt(L) line
                    y_values = [y_val*math.sqrt(analyzed_lumi / l) if l > 0 else 0 for l in lumi]
                    y_values_L = [y_val*(analyzed_lumi / l) if l > 0 else 0 for l in lumi]
                    scaled_graph = ROOT.TGraph(len(lumi))
                    scaled_graph_L = ROOT.TGraph(len(lumi))
                    for i, (lx, ly) in enumerate(zip(lumi, y_values)):
                        scaled_graph.SetPoint(i, lx, ly)
                    for i, (lx, ly) in enumerate(zip(lumi, y_values_L)):
                        scaled_graph_L.SetPoint(i, lx, ly)
                
                    scaled_graph.SetLineColor(ROOT.kBlue)
                    scaled_graph.SetLineWidth(2)
                    scaled_graph.SetLineStyle(3)  # Dashed line
                    scaled_graph.Draw("Lsame")
                    
                    scaled_graph_L.SetLineColor(ROOT.kRed)
                    scaled_graph_L.SetLineWidth(2)
                    scaled_graph_L.SetLineStyle(3)  # Dashed line
                    scaled_graph_L.Draw("Lsame")
                    
            
            x1 = 0.65
            x2 = x1 + 0.24
            y2 = 0.86
            y1 = 0.70
            legend = TLegend(x1,y1,x2,y2)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetTextSize(0.041)
            legend.SetTextFont(42)
            
            for cat_sub in range(Cat_No_sub[cat]):
                    legend.AddEntry(median[0][cat_sub], subcat_label[0][cat_sub],'L')
            
            if categories[cat] == 'HF':
                    legend.AddEntry(scaled_graph, "sqrt(1/L)", 'L')
                    legend.AddEntry(scaled_graph_L, "(1/L)", 'L')
            
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
            
            CMSStyle.CMS_lumi(c, 5, 0)
        
            latex.SetTextAlign(1)
            latex.DrawLatex(0.15, 0.85, Text)
            latex.Draw('same') 
            print(" ")
            if len(categories) < 2:
                    c.SaveAs("Limit_scan_"+categories[0]+".png")
                    c.SaveAs("Limit_scan_"+categories[0]+".pdf")
            else:
                    c.SaveAs("Limit_scan.png")
            c.Close()
            
            
            
            
            
            
            
            for cat in range(Cat_No):
                    if categories[cat] == 'HF':
                        lumi_vals = lumi_HF
                    elif categories[cat] == 'ZTT':
                        lumi_vals = lumi_ZTT
                    elif categories[cat] == 'W':
                        lumi_vals = lumi_W
                    else:
                        lumi_vals = lumi
                
                    if categories[cat] in ['ZTT', 'W']:
                        for cat_sub in range(Cat_No_sub[cat]):
                            graph_med = median[cat][cat_sub]
                            graph_1s = green[cat][cat_sub]
                            graph_2s = yellow[cat][cat_sub]
                
                            N = graph_med.GetN()
                            if N * 2 != graph_1s.GetN() or N * 2 != graph_2s.GetN():
                                print(f"Warning: Mismatch in points for {categories[cat]}_{subcat_label[cat][cat_sub]}")
                                continue
                
                            output_filename = f"median_limits_{categories[cat]}_{subcat_label[cat][cat_sub]}.txt"
                            with open(output_filename, "w") as f:
                                f.write("# lumi [fb^-1]   expected_limit (10^-7)   -2sigma   -1sigma   +1sigma   +2sigma\n")
                                for i in range(N):
                                    lumi_val = graph_med.GetX()[i]
                                    med = graph_med.GetY()[i]
                
                                    plus1sigma = graph_1s.GetY()[i]
                                    minus1sigma = graph_1s.GetY()[2*N - 1 - i]
                
                                    plus2sigma = graph_2s.GetY()[i]
                                    minus2sigma = graph_2s.GetY()[2*N - 1 - i]
                
                                    f.write(f"{lumi_val:.4f}   {med:.6f}   {minus2sigma:.6f}   {minus1sigma:.6f}   {plus1sigma:.6f}   {plus2sigma:.6f}\n")
                
                    else:
                        graph_med = median[cat][0]
                        graph_1s = green[cat][0]
                        graph_2s = yellow[cat][0]
                
                        N = graph_med.GetN()
                        if N * 2 != graph_1s.GetN() or N * 2 != graph_2s.GetN():
                            print(f"Warning: Mismatch in points for {categories[cat]}")
                            continue
                
                        output_filename = f"median_limits_{categories[cat]}.txt"
                        with open(output_filename, "w") as f:
                            f.write("# lumi [fb^-1]   expected_limit (10^-7)   -2sigma   -1sigma   +1sigma   +2sigma\n")
                            for i in range(N):
                                lumi_val = graph_med.GetX()[i]
                                med = graph_med.GetY()[i]
                
                                plus1sigma = graph_1s.GetY()[i]
                                minus1sigma = graph_1s.GetY()[2*N - 1 - i]
                
                                plus2sigma = graph_2s.GetY()[i]
                                minus2sigma = graph_2s.GetY()[2*N - 1 - i]
                
                                f.write(f"{lumi_val:.4f}   {med:.6f}   {minus2sigma:.6f}   {minus1sigma:.6f}   {plus1sigma:.6f}   {plus2sigma:.6f}\n")
            
            
            
            
            
    
    #To plot multiple sub-categories separately as separate images
    if((not WhetherMultipleBroadCategories) and WhetherMultipleSmallCategories and WhetherIndividualPlots):
            
            if categories[cat]=='HF':
                    lumi = lumi_HF
                    analyzed_lumi = analyzed_lumi_HF
            if categories[cat]=='ZTT':
                    lumi = lumi_ZTT
                    analyzed_lumi = analyzed_lumi_ZTT
            if categories[cat]=='W':
                    lumi = lumi_W
                    analyzed_lumi = analyzed_lumi_W
            
            for cat_sub in range(Cat_No_sub[0]):
            
                    c = TCanvas("c","c",100,100,W_ref,H_ref)
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
                    frame.GetYaxis().SetTitleOffset(1.1)
                    frame.GetXaxis().SetNdivisions(508)
                    frame.GetYaxis().CenterTitle(False)
                    frame.GetYaxis().SetTitle("Expected Limit (10^{-7})")
                    frame.GetXaxis().SetTitle("Integrated Luminosity fb^{-1}")
                
                
                
                
                    #frame.SetMinimum(median[cat].GetHistogram().GetMinimum()*0.8)
                    #frame.SetMaximum(median[cat].GetHistogram().GetMinimum()+(5.0/3.0)*(median[cat].GetHistogram().GetMaximum()-median[cat].GetHistogram().GetMinimum()))
                    
                    frame.SetMinimum(0.01)
                    frame.SetMaximum(25.0)
                    
                    frame.GetXaxis().SetLimits(min(lumi) ,max(lumi)*1.2)
                    c.SetLogx();
                    c.cd()
                
                
                    #for cat in range(Cat_No):
                    yellow[cat][cat_sub].SetFillColor(ROOT.kOrange)
                    yellow[cat][cat_sub].SetLineColor(ROOT.kOrange)
                    yellow[cat][cat_sub].SetFillStyle(1001)
                    yellow[cat][cat_sub].Draw('F')
                 
                    green[cat][cat_sub].SetFillColor(ROOT.kGreen+1)
                    green[cat][cat_sub].SetLineColor(ROOT.kGreen+1)
                    green[cat][cat_sub].SetFillStyle(1001)
                    green[cat][cat_sub].Draw('Fsame')
                 
                    median[cat][cat_sub].SetLineColor(1)
                    median[cat][cat_sub].SetLineWidth(2)
                    median[cat][cat_sub].SetLineStyle(2)
                    median[cat][cat_sub].Draw('Lsame')
                    median[cat][cat_sub].Draw()
                    
                    
                    x, y = ctypes.c_double(0), ctypes.c_double(0)  # Create c_double variables for x and y
                    median[cat][cat_sub].GetPoint(0, x, y)  # Get the point at index i
                    # Convert ctypes values to Python floats
                    x_val = x.value
                    y_val = y.value
                    
                    
                    # Create and draw the scaled 1/sqrt(L) line
                    y_values = [y_val*math.sqrt(analyzed_lumi / l) if l > 0 else 0 for l in lumi]
                    y_values_L = [y_val*(analyzed_lumi / l) if l > 0 else 0 for l in lumi]
                    scaled_graph = ROOT.TGraph(len(lumi))
                    scaled_graph_L = ROOT.TGraph(len(lumi))
                    for i, (lx, ly) in enumerate(zip(lumi, y_values)):
                        scaled_graph.SetPoint(i, lx, ly)
                    for i, (lx, ly) in enumerate(zip(lumi, y_values_L)):
                        scaled_graph_L.SetPoint(i, lx, ly)
                
                    scaled_graph.SetLineColor(ROOT.kBlue)
                    scaled_graph.SetLineWidth(2)
                    scaled_graph.SetLineStyle(3)  # Dashed line
                    scaled_graph.Draw("Lsame")
                    
                    scaled_graph_L.SetLineColor(ROOT.kRed)
                    scaled_graph_L.SetLineWidth(2)
                    scaled_graph_L.SetLineStyle(3)  # Dashed line
                    scaled_graph_L.Draw("Lsame")
                    
                    
                 
                #    CMS_lumi.CMS_lumi(c,13,11)
                    ROOT.gPad.SetTicks(1,1)
                #    CMS_lumi(ROOT.gPad, 5, 0)
                    ROOT.gPad.Update()
                    frame.Draw('sameaxis')
                 
                    x1 = 0.65
                    x2 = x1 + 0.24
                    y2 = 0.86
                    y1 = 0.70
                    legend = TLegend(x1,y1,x2,y2)
                    legend.SetFillStyle(0)
                    legend.SetBorderSize(0)
                    legend.SetTextSize(0.041)
                    legend.SetTextFont(42)
                    legend.AddEntry(median[0][cat_sub], subcat_label[0][cat_sub],'L')
                    legend.AddEntry(scaled_graph, "sqrt(1/L)", 'L')
                    legend.AddEntry(scaled_graph_L, "(1/L)", 'L')
                    
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
                    
                    CMSStyle.CMS_lumi(c, 5, 0)
                
                    latex.SetTextAlign(1)
                    latex.DrawLatex(0.15, 0.85, Text)
                    latex.Draw('same') 
                    print(" ")
                    c.SaveAs("Limit_scan_"+categories[0]+subcat[0][cat_sub]+".png")
                    c.SaveAs("Limit_scan_"+categories[0]+subcat[0][cat_sub]+".png")
                    
                    c.Close()
            
            
    
    
    #For multiple categories in a single plot
    if(WhetherMultipleBroadCategories):
            
            print("WhetherMultipleBroadCategories 2: ",WhetherMultipleBroadCategories)
            
            c = TCanvas("c","c",100,100,W_ref,H_ref)
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
            frame.GetYaxis().SetTitleOffset(1.1)
            frame.GetXaxis().SetNdivisions(508)
            frame.GetYaxis().CenterTitle(False)
            frame.GetYaxis().SetTitle("Expected Limit (10^{-7})")
            frame.GetXaxis().SetTitle("Integrated Luminosity fb^{-1}")
        
        
        
        
            #frame.SetMinimum(median[cat].GetHistogram().GetMinimum()*0.8)
            #frame.SetMaximum(median[cat].GetHistogram().GetMinimum()+(5.0/3.0)*(median[cat].GetHistogram().GetMaximum()-median[cat].GetHistogram().GetMinimum()))
            
            frame.SetMinimum(0.01)
            frame.SetMaximum(25.0)
            
            frame.GetXaxis().SetLimits(min(lumi) ,max(lumi)*1.2)
            
            c.SetLogx();
            c.cd()
        
        
            for cat in range(Cat_No):
                    #yellow[cat][cat_sub].SetFillColor(ROOT.kOrange)
                    #yellow[cat][cat_sub].SetLineColor(ROOT.kOrange)
                    #yellow[cat][cat_sub].SetFillStyle(1001)
                    #yellow[cat][cat_sub].Draw('F')
                 
                    #green[cat][cat_sub].SetFillColor(ROOT.kGreen+1)
                    #green[cat][cat_sub].SetLineColor(ROOT.kGreen+1)
                    #green[cat][cat_sub].SetFillStyle(1001)
                    #green[cat][cat_sub].Draw('Fsame')
                 
                    median[cat][0].SetLineColor(cat+1)
                    median[cat][0].SetLineWidth(2)
                    median[cat][0].SetLineStyle(2)
                    median[cat][0].Draw('Lsame')
                    median[cat][0].Draw()
         
        #    CMS_lumi.CMS_lumi(c,13,11)
            ROOT.gPad.SetTicks(1,1)
        #    CMS_lumi(ROOT.gPad, 5, 0)
            ROOT.gPad.Update()
            frame.Draw('sameaxis')
         
            x1 = 0.65
            x2 = x1 + 0.24
            y2 = 0.86
            y1 = 0.70
            legend = TLegend(x1,y1,x2,y2)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetTextSize(0.041)
            legend.SetTextFont(42)
            for cat in range(Cat_No):
                    legend.AddEntry(median[cat][0], label[cat],'L')
            
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
            
            CMSStyle.CMS_lumi(c, 5, 0)
        
            latex.SetTextAlign(1)
            latex.DrawLatex(0.15, 0.85, Text)
            latex.Draw('same') 
            print(" ")
            c.SaveAs("Limit_scan.png")
            c.Close()
            
            
            for cat in range(Cat_No):
                    graph = median[cat][0]
                    N = graph.GetN()
                    x_vals = graph.GetX()
                    y_vals = graph.GetY()
                    
                    # Save to a text file named after the category
                    output_filename = f"median_limits_cat_{label[cat]}.txt"
                    with open(output_filename, "w") as f:
                        f.write("# lumi [fb^-1]   expected_limit (10^-7)\n")
                        for i in range(N):
                            f.write(f"{x_vals[i]:.4f}   {y_vals[i]:.6f}\n")
            
        


# MAIN
def main():


        
#    categories = ['ZTT','W']
#    categories = ['W']
    categories = ['ZTT']
#    categories = ['HF']
#    categories = ['ZTT','HF','W']
    
    Whether_Hybrid=False
    if(args.Method == 'H'):
            Whether_Hybrid=True
    
    #actual categories
    #categories=[args.category]
    
    #print(lumi_ZTT)
    
    if(args.runType == 'run'):
            executeDataCards_onCondor(lumi_W,lumi_HF,lumi_ZTT,categories,Whether_Hybrid)
    if(args.runType == 'plot'):
            plotUpperLimits(lumi_W,lumi_HF,lumi_ZTT,analyzed_lumi_W,analyzed_lumi_HF,analyzed_lumi_ZTT,categories,Whether_Hybrid,False)

if __name__ == '__main__':
    main()
