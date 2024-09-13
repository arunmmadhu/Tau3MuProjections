#need to cmsenv somewhere first

import os
import numpy as np

lim_mid = []
lim_up = []
lim_up_2 = []
lim_low = []
lim_low_2 = []
lim = [lim_low_2,lim_low,lim_mid,lim_up,lim_up_2]

#lumi = [59.0,97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0]
lumi = [59.0,97.7, 129.0]

def executeDataCards_onCondor(lumi):
        for lu_no in range(len(lumi)):
                
            lu = lumi[lu_no]
            print(lu)
            Whether_Hybrid=True
                
            if(Whether_Hybrid):
                    
                    print("Running Sigma -2")
                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.025 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu),"datacards_modified/dc_"+str(lu)+".txt",str(lu))
                    os.system(command_run)
                    print("Running Sigma -1")
                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.16 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu),"datacards_modified/dc_"+str(lu)+".txt",str(lu))
                    os.system(command_run)
                    print("Running Sigma Median")
                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.5 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu),"datacards_modified/dc_"+str(lu)+".txt",str(lu))
                    os.system(command_run)
                    print("Running Sigma +1")
                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.84 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu),"datacards_modified/dc_"+str(lu)+".txt",str(lu))
                    os.system(command_run)
                    print("Running Sigma +2")
                    command_run = "combineTool.py -M HybridNew --LHCmode LHC-limits  -n %s -d %s --rMin 0 --rMax 50 --cl 0.90 -t 10 --expectedFromGrid 0.975 --job-mode condor --sub-opts='+JobFlavour=\"workday\"'  --task-name HybridTest%s " % (str(lu),"datacards_modified/dc_"+str(lu)+".txt",str(lu))
                    os.system(command_run)
            
        #lim = [lim_mid,lim_low_2,lim_low,lim_up,lim_up_2]
        #print(lim)

# PLOT upper limits
def plotUpperLimits(lumi):
 
    N = len(lumi)
    yellow = TGraph(2*N)   
    green = TGraph(2*N)    
    median = TGraph(N)  
    
    text_limits=open("TextLimits%s"%(prefix)+outputLabel+".txt","w")
    for i in range(N):
        file_name1 = "higgsCombine"+str(lumi[i])+".HybridNew.mH120.123456.quant0.025.root"
        limit1 = getLimits(file_name1)
        file_name2 = "higgsCombine"+str(lumi[i])+".HybridNew.mH120.123456.quant0.160.root"
        limit2 = getLimits(file_name2)
        file_name3 = "higgsCombine"+str(lumi[i])+".HybridNew.mH120.123456.quant0.500.root"
        limit3 = getLimits(file_name3)
        file_name4 = "higgsCombine"+str(lumi[i])+".HybridNew.mH120.123456.quant0.840.root"
        limit4 = getLimits(file_name4)
        file_name5 = "higgsCombine"+str(lumi[i])+".HybridNew.mH120.123456.quant0.975.root"
        limit5 = getLimits(file_name5)
        
        yellow.SetPoint( 2*N-1-i, lumi[i], limit1[2]) # - 2 sigma
        green.SetPoint(  2*N-1-i, lumi[i], limit2[2]) # - 1 sigma
        median.SetPoint(    i,    lumi[i], limit3[2]) #    median
        green.SetPoint(     i,    lumi[i], limit4[2]) # + 1 sigma
        yellow.SetPoint(    i,    lumi[i], limit5[2]) # + 2 sigma
        
        text_limits.write("bdt %.2f     median exp %.2f\n"%(lumi[i],limit3[2]))

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
    frame.GetYaxis().SetTitle("Expected Limit ($10^{-7}$)")
    frame.GetXaxis().SetTitle("Integrated Luminosity $fb^{-1}$")




    #frame.SetMinimum(median.GetHistogram().GetMinimum()*0.8)
    #frame.SetMaximum(median.GetHistogram().GetMinimum()+(5.0/3.0)*(median.GetHistogram().GetMaximum()-median.GetHistogram().GetMinimum()))
    
    frame.SetMinimum(0.0)
    frame.SetMaximum(25.0)

#    frame.GetXaxis().SetLimits(min(values),max(values)*1.2)
    frame.GetXaxis().SetLimits(min(values) ,max(values)*1.2)


 
    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
#    yellow.Draw('F')
 
    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
#    green.Draw('Fsame')
 
    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')
    median.Draw()
 
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
    legend.AddEntry(median, "Asymptotic CL_{s} expected upper limit",'L')
#    legend.AddEntry(median, "HybridNew CL_{s} expected upper limit",'L')
#    legend.AddEntry(green, "#pm 1 std. deviation",'f')
#    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')


    legend.Draw()
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    Text = ''
    if prefix=='taue':
        Text = 'Category: Z#rightarrow#tau_{e}#tau_{3#mu}'

    if prefix=='taumu':
        Text = 'Category: Z#rightarrow#tau_{#mu}#tau_{3#mu}'

    if prefix=='tauhA':
        Text = 'Category: Z#rightarrow#tau_{h,1-prong}#tau_{3#mu}'

    if prefix=='tauhB':
        Text = 'Category: Z#rightarrow#tau_{h,3-prong}#tau_{3#mu}'
        
    if prefix=='all':
        Text = 'Category: Z#rightarrow#tau#tau_{3#mu}'

    latex.SetTextAlign(1)
    latex.DrawLatex(0.15, 0.85, Text)
    latex.Draw('same') 
    print " "
    c.SaveAs("Limit_scan_Category_"+prefix+outputLabel+".png")
    c.Close()
    
#    executeDataCards_onCondor
        
    plotUpperLimits(lumi)
