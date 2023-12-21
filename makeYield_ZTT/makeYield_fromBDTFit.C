#include "RooGlobalFunc.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
//#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TH1D.h"
#include "TH2D.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>

using namespace RooFit;

void makeYield_fromBDTFit () 
{
    //Category names: 0=tauh, 1=taumu, 2=taue
    TString cat_name[3];
    cat_name[0] = "ztau3mutauh_default_";
    cat_name[1] = "ztau3mutaumu_default_";
    cat_name[2] = "ztau3mutaue_default_";
    
    TString cat_base[3];
    cat_base[0] = "ztau3mutauh";
    cat_base[1] = "ztau3mutaumu";
    cat_base[2] = "ztau3mutaue";
    
    TString cat_label[3];
    cat_label[0] = "{h}";
    cat_label[1] = "{#mu}";
    cat_label[2] = "{e}";
    
    TString print_label[3];
    print_label[0] = "h";
    print_label[1] = "mu";
    print_label[2] = "e";
    
    TString hname;
    
    float signal_region_min(1.4);
    float signal_region_max(2.1);
    
    float signal_peak_region_min(1.73);
    float signal_peak_region_max(1.82);
    
    float Loose_BDT_Cut(-0.4);
    
    //Filename and histograms
    TFile * file_tau[3];
    
    TH1D  * tau_T3Mu[3];
    TH1D  * tau_T3Mu_Dat[3];
    TH1D  * tau_BDT_Output_Data[3];
    TH1D  * tau_BDT_Output_MC[3];
    TH2D  * tau_T3Mu_vs_BDT_Output_Data[3];
    TH1D  * tau_T3Mu_vs_BDT_Output_Data_Projection[3];
    
    TH2D  * tau_cut1_vs_cut2_vs_sig[3];
    TH2D  * tau_cut1_vs_cut2_vs_limit[3];
    TH1D  * tau_loose_cut_vs_ratio[3];
    TH1D  * tau_loose_cut_vs_ratio_num[3];
    TH1D  * tau_loose_cut_vs_ratio_den[3];
    
    
    
    
    
    /*
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      file_tau[i] = TFile::Open("LOCAL_COMBINED_"+cat_name[i]+"LumiScaled_Lim.root");
      
      tau_T3Mu[i] = (TH1D*)file_tau[i]->Get(cat_name[i]+"PostBDT_TripletMass_VeryLooseCutMC"+hname);
      tau_T3Mu_Dat[i] = (TH1D*)file_tau[i]->Get(cat_name[i]+"PostBDT_TripletMass_VeryLooseCutData");
      tau_BDT_Output_Data[i] = (TH1D*)file_tau[i]->Get(cat_name[i]+"PostSelection_BDT_OutputData");
      tau_BDT_Output_MC[i] = (TH1D*)file_tau[i]->Get(cat_name[i]+"PostSelection_BDT_OutputMC"+hname);
      tau_T3Mu_vs_BDT_Output_Data[i] = (TH2D*)file_tau[i]->Get(cat_name[i]+"BDT_2Dscan_TripletMassData");
      
      
    }
    */
    
    TFile *TreeFile = new TFile("Combine_Tree_ztau3mutau.root","READ");
    TTree *tree[3];
    
    Float_t tripletMass;
    Float_t bdt_cv;
    Float_t weight;
    Float_t isMC;
    Float_t ifCommonCV;
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      tree[i] = (TTree *) TreeFile->Get(cat_base[i]);
      
      tau_T3Mu[i] = new TH1D("tau_T3Mu","tau_T3Mu_"+hname,40,signal_region_min,signal_region_max);
      tau_T3Mu_Dat[i] = new TH1D("tau_T3Mu_Dat","tau_T3Mu_Dat_"+hname,40,signal_region_min,signal_region_max);
      tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,100,-0.9,0.9);
      tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,100,-0.9,0.9);
      
      tree[i]->SetBranchAddress("tripletMass",&tripletMass);
      tree[i]->SetBranchAddress("bdt_cv",&bdt_cv);
      tree[i]->SetBranchAddress("weight",&weight);
      tree[i]->SetBranchAddress("isMC",&isMC);
      tree[i]->SetBranchAddress("ifCommonCV",&ifCommonCV);
      
      Long64_t nentries = tree[i]->GetEntries();
      for (Long64_t j=0;j<nentries;j++) {
        tree[i]->GetEntry(j);
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max&&ifCommonCV){
        //if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max){
                if(isMC>0){
                  if(bdt_cv>Loose_BDT_Cut){
                    tau_T3Mu[i]->Fill(tripletMass,weight);
                  }
                  tau_BDT_Output_MC[i]->Fill(bdt_cv,weight);
                }
                if(isMC==0 && (tripletMass<=signal_peak_region_min || tripletMass>=signal_peak_region_max) ){//blinded
                  if(bdt_cv>Loose_BDT_Cut){
                    tau_T3Mu_Dat[i]->Fill(tripletMass);
                  }
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      //tau_T3Mu_vs_BDT_Output_Data[i] = (TH2D*)file_tau[i]->Get(cat_name[i]+"BDT_2Dscan_TripletMassData");
      
      
    }
    
    for(int i=0; i<3; i++){
            tau_T3Mu[i]->Sumw2();
            tau_T3Mu_Dat[i]->Sumw2();
            tau_BDT_Output_Data[i]->Sumw2();
            tau_BDT_Output_MC[i]->Sumw2();
    }
    
    
    
    
    //double sig_norm = 0.0000283 * 0.215; //average normalization factor for the three signal samples
    //tau_BDT_Output_MC[1]->Scale(sig_norm);
    
    
    
    
    /*
    //For fitting BDT Output: Using crystal ball
    //General
    RooRealVar * BDTOutput_x[3];
    RooRealVar * cbmean[3];
    RooRealVar * cbsigma[3];
    RooRealVar * n[3];
    RooRealVar * alpha[3];
    
    RooCBShape * cball[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, #tau_"+cat_label[i],-0.9,0.9);
      BDTOutput_x[i]->setRange("R1",-0.9,0.9);
      cbmean[i] = new RooRealVar("cbmean"+hname, "cbmean" , -0.5, -0.9,0.0) ;
      cbsigma[i] = new RooRealVar("cbsigma"+hname, "cbsigma" , 0.2, 0.000001, 1.0) ;
      n[i] = new RooRealVar("n"+hname, "n", 15.0, 0.5, 20);
      alpha[i] = new RooRealVar("alpha"+hname,"alpha value CB",5.0,0.01,10);
      
      cball[i] = new RooCBShape("cball"+hname, "crystal ball", *BDTOutput_x[i], *cbmean[i], *cbsigma[i], *alpha[i], *n[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*cball[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("R1"), Save());
      
    }
    */
    
    
    
    
    
    
    //For fitting BDT Output in Data: Using Bifurgauss
    //General
    RooRealVar * BDTOutput_x[3];
    RooRealVar * bgausmean[3];
    RooRealVar * bgaussigma_a[3];
    RooRealVar * bgaussigma_b[3];
    
    
    RooBifurGauss * bgaus_dist[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, #tau_"+cat_label[i],-0.9,0.9);
      BDTOutput_x[i]->setRange("BDT_Range",-0.9,0.9);
      bgausmean[i] = new RooRealVar("bgausmean"+hname, "bgausmean" , -0.5, -0.9,0.0) ;
      bgaussigma_a[i] = new RooRealVar("bgaussigma_a"+hname, "bgaussigma_a" , 0.2, 0.000001, 1.0) ;
      bgaussigma_b[i] = new RooRealVar("bgaussigma_b"+hname, "bgaussigma_b" , 0.2, 0.000001, 1.0);
      
      bgaus_dist[i] = new RooBifurGauss("bgaus_dist"+hname, "bgaus dist", *BDTOutput_x[i], *bgausmean[i], *bgaussigma_a[i], *bgaussigma_b[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*bgaus_dist[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    
    
    
    /*
    //BDT Data Fit Plots
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1800, 600);
    canvas2->Divide(3, 1);
    
    RooPlot * xFrame_bdt[3];
    float BDT_Data_Max[3];
    BDT_Data_Max[0] = 500.0;
    BDT_Data_Max[1] = 3000.0;
    BDT_Data_Max[2] = 350.0;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas2->cd( i+1 );
      xFrame_bdt[i] = BDTOutput_x[i]->frame();
      bdt_data[i]->plotOn(xFrame_bdt[i]);
      BDT_distribution[i]->plotOn(xFrame_bdt[i],LineColor(4),LineWidth(2), Normalization(bdt_data[i]->sumEntries("1", "R1"), RooAbsReal::NumEvent),ProjectionRange("R1"));
      xFrame_bdt[i]->SetTitle("BDT Output, #tau_"+cat_label[i]);
      xFrame_bdt[i]->SetXTitle("BDT Score");
      xFrame_bdt[i]->SetYTitle("Events");
      xFrame_bdt[i]->GetXaxis()->SetRangeUser(-0.9,0.9);
      xFrame_bdt[i]->GetYaxis()->SetRangeUser(0.0001,BDT_Data_Max[i]);
      canvas2->GetPad(i+1)->SetLogy();
      xFrame_bdt[i]->Draw();
    }
    */
    
    
    
    
    
    
    
    
    //For fitting BDT Output in MC: Using Bifurgauss
    //General
    RooRealVar * bgausmeanMC[3];
    RooRealVar * bgaussigmaMC_a[3];
    RooRealVar * bgaussigmaMC_b[3];
    
    
    RooBifurGauss * bgaus_distMC[3];
    RooDataHist * bdt_MC[3];
    RooRealVar * BDTNormMC[3];
    RooAddPdf * BDT_distributionMC[3];
    RooFitResult * fitresult_bdtMC[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      bgausmeanMC[i] = new RooRealVar("bgausmeanMC"+hname, "bgausmeanMC" , 0.5, 0.0,0.9) ;
      bgaussigmaMC_a[i] = new RooRealVar("bgaussigmaMC_a"+hname, "bgaussigmaMC_a" , 0.2, 0.000001, 1.0) ;
      bgaussigmaMC_b[i] = new RooRealVar("bgaussigmaMC_b"+hname, "bgaussigmaMC_b" , 0.2, 0.000001, 1.0);
      
      bgaus_distMC[i] = new RooBifurGauss("bgaus_distMC"+hname, "bgaus dist MC", *BDTOutput_x[i], *bgausmeanMC[i], *bgaussigmaMC_a[i], *bgaussigmaMC_b[i]);
      bdt_MC[i] = new RooDataHist("bdt_MC"+hname, "bdt_MC", *BDTOutput_x[i], Import(*tau_BDT_Output_MC[i]));
      BDTNormMC[i] = new RooRealVar("BDTNormMC"+hname, "BDTNormMC", 5.0,0.0,50);
      BDT_distributionMC[i] = new RooAddPdf("BDT_distributionMC"+hname, "BDT_distributionMC", RooArgList(*bgaus_distMC[i]), RooArgList(*BDTNormMC[i]));
      fitresult_bdtMC[i] = BDT_distributionMC[i]->fitTo(*bdt_MC[i], Range("BDT_Range"), Save());
      
    }
    
    
    
    
    /*
    //MC BDT Fit Plots
    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 1800, 600);
    canvas3->Divide(3, 1);
    
    RooPlot * xFrame_bdtMC[3];
    //float BDT_MC_Max[3];
    //BDT_MC_Max[0] = 100.0;
    //BDT_MC_Max[1] = 400.0;
    //BDT_MC_Max[2] = 50.0;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas3->cd( i+1 );
      xFrame_bdtMC[i] = BDTOutput_x[i]->frame();
      bdt_MC[i]->plotOn(xFrame_bdtMC[i]);
      BDT_distributionMC[i]->plotOn(xFrame_bdtMC[i],LineColor(4),LineWidth(2), Normalization(bdt_MC[i]->sumEntries("1", "R1"), RooAbsReal::NumEvent),ProjectionRange("R1"));
      xFrame_bdtMC[i]->SetTitle("BDT Output, #tau_"+cat_label[i]);
      xFrame_bdtMC[i]->SetXTitle("BDT Score");
      xFrame_bdtMC[i]->SetYTitle("Events");
      xFrame_bdtMC[i]->GetXaxis()->SetRangeUser(-0.9,0.9);
      //xFrame_bdtMC[i]->GetYaxis()->SetRangeUser(0.0000001,BDT_Data_Max[i]);
      //canvas2->GetPad(i+1)->SetLogy();
      xFrame_bdtMC[i]->Draw();
    }
    */
    
    
    
    
    
    
    
    

    //cout << "tauh MC count: " << tauh_Vis_Mass->Integral() << " tauh Data count: " << tauh_Vis_Mass_Dat->Integral() << endl;
    //cout << "taumu MC count: " << taumu_Vis_Mass->Integral() << " taumu Data count: " << taumu_Vis_Mass_Dat->Integral() << endl;
    //cout << "taue MC count: " << taue_Vis_Mass->Integral() << " taue Data count: " << taue_Vis_Mass_Dat->Integral() << endl;
    
    
    
    
    
    //Triplet Mass Fits
    RooRealVar * InvMass[3];
    
    RooPolynomial * poly[3];
    RooDataHist * data[3];
    RooRealVar * LineNorm[3];
    RooAddPdf * pdf[3];
    RooFitResult * fitresult[3];
    
    RooRealVar * mean[3];
    RooRealVar * sigma[3];
    RooGaussian * Gauss[3];
    RooDataHist * mc[3];
    RooRealVar * GaussNorm[3];
    RooAddPdf * mc_pdf[3];
    RooFitResult * mc_fitresult[3];
    
    RooExponential * Expo[3];
    RooRealVar * lambda[3];
    RooRealVar * ExpNorm[3];
    RooAddPdf * exp_pdf[3];
    RooFitResult * exp_fitresult[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      InvMass[i] = new RooRealVar("InvMass"+hname,"InvMass, #tau_"+cat_label[i],signal_region_min,signal_region_max);
      InvMass[i]->setRange("R1",signal_region_min,signal_peak_region_min); //background   
      InvMass[i]->setRange("R2",signal_peak_region_max,signal_region_max); //background
      InvMass[i]->setRange("R3",1.72,1.81); //signal range for fitting
      InvMass[i]->setRange("R4",signal_peak_region_min,signal_peak_region_max); //signal range for yield
      
      
      //Flat fit for data
      poly[i] = new RooPolynomial("poly"+hname, "poly dist", *InvMass[i]);
      data[i] = new RooDataHist("data"+hname, "data", *InvMass[i], Import(*tau_T3Mu_Dat[i]));
      LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,15);
      pdf[i] = new RooAddPdf("pdf"+hname, "pdf", RooArgList(*poly[i]), RooArgList(*LineNorm[i]));
      fitresult[i] = pdf[i]->fitTo(*data[i], Range("R1,R2"), Save());
      
      //Testing Exponential Fit
      lambda[i] = new RooRealVar("lambda"+hname,"lambda",0.1, -1.5, 1.5);
      Expo[i] = new RooExponential("Expo"+hname, "Exponential PDF", *InvMass[i],  *lambda[i]);
      ExpNorm[i] = new RooRealVar("ExpNorm"+hname, "ExpNorm",  2.0,0.001,15);
      exp_pdf[i] = new RooAddPdf("exp_pdf"+hname, "exp_pdf", RooArgList(*Expo[i]), RooArgList(*ExpNorm[i]));
      exp_fitresult[i] = exp_pdf[i]->fitTo(*data[i], Range("R1,R2"), Save());
      
      
      //Gaussian fit for MC
      mean[i] = new RooRealVar("mean"+hname, "mean" , 1.776,0.,5.) ;
      sigma[i] = new RooRealVar("sigma"+hname, "sigma" , 0.5,0.001,10) ;
      
      Gauss[i] = new RooGaussian("Gauss"+hname, "Gauss dist", *InvMass[i], *mean[i], *sigma[i]);
      mc[i] = new RooDataHist("mc"+hname, "mc", *InvMass[i], Import(*tau_T3Mu[i]));
      GaussNorm[i] = new RooRealVar("GaussNorm"+hname, "GaussNorm",  0.5,0.001,1.0);
      mc_pdf[i] = new RooAddPdf("mc_pdf"+hname, "mc_pdf", RooArgList(*Gauss[i]), RooArgList(*GaussNorm[i]));
      mc_fitresult[i] = mc_pdf[i]->fitTo(*mc[i], Range("R3"), Save());
      
    }
    
    
    
    
    
    /*
    // This gives the integral from the fits.
    double pdf1_integral_restricted = pdf1.createIntegral(InvMass1,NormSet(InvMass1),Range("R4"))->getVal();
    double mc_pdf1_integral_restricted = mc_pdf1.createIntegral(InvMass1,NormSet(InvMass1),Range("R4"))->getVal();
    double pdf2_integral_restricted = pdf2.createIntegral(InvMass2,NormSet(InvMass2),Range("R4"))->getVal();
    double mc_pdf2_integral_restricted = mc_pdf2.createIntegral(InvMass2,NormSet(InvMass2),Range("R4"))->getVal();
    double pdf3_integral_restricted = pdf3.createIntegral(InvMass3,NormSet(InvMass3),Range("R4"))->getVal();
    double mc_pdf3_integral_restricted = mc_pdf3.createIntegral(InvMass3,NormSet(InvMass3),Range("R4"))->getVal();
    
    // Normalizations need to be added manually. The pdfs are normalized to 1 and scaled to the data plotted. nData1 and nSignal1 are for normalization (same region as fit). "R4" is for yields.
    const double nData1 = data1.sumEntries("1", "R1,R2");
    const double nSignal1 = mc1.sumEntries("1", "R3");
    const double nSignal1_restricted = mc1.sumEntries("1", "R4");// used a separate range for getting the yield and a different range for fitting
    const double nData2 = data2.sumEntries("1", "R1,R2");
    const double nSignal2 = mc2.sumEntries("1", "R3");
    const double nSignal2_restricted = mc2.sumEntries("1", "R4");
    const double nData3 = data3.sumEntries("1", "R1,R2");
    const double nSignal3 = mc3.sumEntries("1", "R3");
    const double nSignal3_restricted = mc3.sumEntries("1", "R4");
    
    cout << "nData1: " << nData1 << " nSignal1: " << nSignal1_restricted << endl;
    cout << "nData2: " << nData2 << " nSignal2: " << nSignal2_restricted << endl;
    cout << "nData3: " << nData3 << " nSignal3: " << nSignal3_restricted << endl;
    
    cout << " LineNorm1: " << LineNorm1.getValV() << " guessed data in signal region: " << (pdf1_integral_restricted/(1-pdf1_integral_restricted))*nData1 << endl;
    cout << " LineNorm1: " << LineNorm2.getValV() << " guessed data in signal region: " << (pdf2_integral_restricted/(1-pdf2_integral_restricted))*nData2 << endl;
    cout << " LineNorm1: " << LineNorm3.getValV() << " guessed data in signal region: " << (pdf3_integral_restricted/(1-pdf3_integral_restricted))*nData3 << endl;
    
    cout << "mc_pdf1_integral: " << GaussNorm1.getValV()*mc_pdf1_integral_restricted << " pdf1_integral: " << LineNorm1.getValV()*pdf1_integral_restricted << endl;
    cout << "mc_pdf2_integral: " << GaussNorm2.getValV()*mc_pdf2_integral_restricted << " pdf2_integral: " << LineNorm2.getValV()*pdf2_integral_restricted << endl;
    cout << "mc_pdf3_integral: " << GaussNorm3.getValV()*mc_pdf3_integral_restricted << " pdf3_integral: " << LineNorm3.getValV()*pdf3_integral_restricted << endl;
    
    double scaling1 = mc1.sumEntries("1")/tauh_T3Mu->GetEntries();
    double scaling2 = mc2.sumEntries("1")/taumu_T3Mu->GetEntries();
    double scaling3 = mc3.sumEntries("1")/taue_T3Mu->GetEntries();
    
    cout << "Unscaled mc1: " << nSignal1_restricted/scaling1 << " scaling 1: " << scaling1 << endl;// Getting the unweighted content of the signal histogram
    cout << "Unscaled mc2: " << nSignal2_restricted/scaling2 << " scaling 2: " << scaling2 << endl;
    cout << "Unscaled mc3: " << nSignal3_restricted/scaling3 << " scaling 3: " << scaling3 << endl;
    */
    
    
    double pdf_integral_restricted[3];
    double mc_pdf_integral_restricted[3];
    double nData[3];
    double nSignal[3];
    double nSignal_restricted[3];
    double scaling[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      // This gives the integral from the fits.
      pdf_integral_restricted[i] = pdf[i]->createIntegral(*InvMass[i],NormSet(*InvMass[i]),Range("R4"))->getVal();
      mc_pdf_integral_restricted[i] = mc_pdf[i]->createIntegral(*InvMass[i],NormSet(*InvMass[i]),Range("R4"))->getVal();
      
      // Normalizations need to be added manually. The pdfs are normalized to 1 and scaled to the data plotted. nData1 and nSignal1 are for normalization (same region as fit). "R4" is for yields.
      nData[i] = data[i]->sumEntries("1", "R1,R2");
      nSignal[i] = mc[i]->sumEntries("1", "R3");
      nSignal_restricted[i] = mc[i]->sumEntries("1", "R4");// used a separate range for getting the yield and a different range for fitting
      
      
      
      cout << "  " << endl;
      
      cout << "nData "+print_label[i]+" : " << nData[i] << " nSignal "+print_label[i]+" : " << nSignal_restricted[i] << endl;
      cout << " LineNorm "+print_label[i]+" : " << LineNorm[i]->getValV() << " guessed data in signal region: " << (pdf_integral_restricted[i]/(1-pdf_integral_restricted[i]))*nData[i] << endl;
      cout << "mc_pdf "+print_label[i]+" _integral: " << GaussNorm[i]->getValV()*mc_pdf_integral_restricted[i] << " pdf "+print_label[i]+" _integral: " << LineNorm[i]->getValV()*pdf_integral_restricted[i] << endl;
      
      scaling[i] = mc[i]->sumEntries("1")/(tau_T3Mu[i]->GetEntries());
      
      cout << "Unscaled mc "+print_label[i]+" : " << nSignal_restricted[i]/scaling[i] << " scaling  "+print_label[i]+" : " << scaling[i] << endl;// Getting the unweighted content of the signal histogram
      
      cout << "  " << endl;
      
    }
    
    
    
    
    /*
    //Triplet Mass Fit Plots
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 1800, 600);
    canvas1->Divide(3, 1);
    
    RooPlot * xFrame[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas1->cd( i+1 );
      xFrame[i] = InvMass[i]->frame();
      data[i]->plotOn(xFrame[i]);
      pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
      //exp_pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));//plots exponential fit
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), #tau_"+cat_label[i]);
      xFrame[i]->SetXTitle("3#mu inv. mass (GeV)");
      xFrame[i]->SetYTitle("Events");
      xFrame[i]->Draw();
    }
    */
    
    
    
    
    
    
    
    
    
    /*
    
    
    double TightBDTCutFraction1 = BDT_distribution1.createIntegral(BDTOutput_x1,NormSet(BDTOutput_x1),Range("R_Small"))->getVal() / BDT_distribution1.createIntegral(BDTOutput_x1,NormSet(BDTOutput_x1),Range("R_Large"))->getVal();
    double TightBDTCutFraction2 = BDT_distribution2.createIntegral(BDTOutput_x2,NormSet(BDTOutput_x2),Range("R_Small"))->getVal() / BDT_distribution2.createIntegral(BDTOutput_x2,NormSet(BDTOutput_x2),Range("R_Large"))->getVal();
    double TightBDTCutFraction3 = BDT_distribution3.createIntegral(BDTOutput_x3,NormSet(BDTOutput_x3),Range("R_Small"))->getVal() / BDT_distribution3.createIntegral(BDTOutput_x3,NormSet(BDTOutput_x3),Range("R_Large"))->getVal();
    
    cout << "Yield mc1: " << BDT_distributionMC1.createIntegral(BDTOutput_x1,NormSet(BDTOutput_x1),Range("R_Small"))->getVal() * BDTNormMC1.getValV() << " Yield Background 1: " << (pdf1_integral_restricted/(1-pdf1_integral_restricted))*nData1*TightBDTCutFraction1 << endl;
    cout << "Yield mc2: " << BDT_distributionMC2.createIntegral(BDTOutput_x2,NormSet(BDTOutput_x2),Range("R_Small"))->getVal() * BDTNormMC2.getValV() << " Yield Background 2: " << (pdf2_integral_restricted/(1-pdf2_integral_restricted))*nData2*TightBDTCutFraction2 << endl;
    cout << "Yield mc3: " << BDT_distributionMC3.createIntegral(BDTOutput_x3,NormSet(BDTOutput_x3),Range("R_Small"))->getVal() * BDTNormMC3.getValV() << " Yield Background 3: " << (pdf3_integral_restricted/(1-pdf3_integral_restricted))*nData3*TightBDTCutFraction3 << endl;
    */
    
    
    
    
    
    
    
    //From BDT Output:
    BDTOutput_x[0]->setRange("R_Small",0.309056,100); BDTOutput_x[0]->setRange("R_Large",Loose_BDT_Cut,100);
    BDTOutput_x[1]->setRange("R_Small",0.329983,100); BDTOutput_x[1]->setRange("R_Large",Loose_BDT_Cut,100);
    BDTOutput_x[2]->setRange("R_Small",0.333186,100); BDTOutput_x[2]->setRange("R_Large",Loose_BDT_Cut,100);
    
    double TightBDTCutFraction[3];
    
    RooAddPdf BDT_distributionTest = *BDT_distribution[0];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      TightBDTCutFraction[i] = BDT_distribution[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Small"))->getVal() / BDT_distribution[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Large"))->getVal();
      cout << "Yield mc "+print_label[i]+" : " << BDT_distributionMC[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Small"))->getVal() * (bdt_MC[i]->sumEntries("1", "R1")) << " Yield Background "+print_label[i]+" : " << (pdf_integral_restricted[i]/(1-pdf_integral_restricted[i]))*nData[i]*TightBDTCutFraction[i] << endl;
      
    }
    
    
    
        //Calculate significance
        double a[3],b[3];
        double N_s_1[3], N_s_2[3];
        double N_b_1[3], N_b_2[3];
        double S1[3], S2[3], S[3];
        std::vector<double> S1_list[3], S2_list[3], S_list[3], a_list[3], b_list[3], N_s_1_yield_list[3], N_s_2_yield_list[3], N_b_1_yield_list[3], N_b_2_yield_list[3];
        
        //double X_min = std::min(tau_BDT_Output_Data[0]->GetXaxis()->GetXmin(), tau_BDT_Output_MC[0]->GetXaxis()->GetXmin());
        //double X_max = std::max(tau_BDT_Output_Data[0]->GetXaxis()->GetXmax(), tau_BDT_Output_MC[0]->GetXaxis()->GetXmax());
        
        double X_min = 0.2;
        double X_max = 0.6;
        
        //Loop on both cuts in [X_min;X_max]
        Int_t dim = 0;
        //Increase N to increase (a,b) scan granularity!
        Int_t N = 30; double step = (X_max - X_min)/N;
        
        for(int i=0; i<3; i++){
          hname=to_string(i+1);
          
          //tau_cut1_vs_cut2_vs_sig[i] = new TH2D("tau_cut1_vs_cut2_vs_sig","tau_cut1_vs_cut2_vs_sig_"+hname,N,X_min,X_max,N,X_min,X_max,"a","b");
          tau_cut1_vs_cut2_vs_sig[i] = new TH2D("tau_cut1_vs_cut2_vs_sig","tau_cut1_vs_cut2_vs_sig_"+hname,N,X_min,X_max,N,X_min,X_max);
              
        }
        
        for(int i=0; i<N; i++){
                for(int j=0; j<N; j++){
                        
                        for(int k=0; k<3; k++){
                                
                                a[k] = X_min + i * step;
                                b[k] = X_min + j * step;
                                if(a[k]<b[k]) continue;
                                
                                
                                
                                
                                BDTOutput_x[k]->setRange("R_a",a[k],100); BDTOutput_x[k]->setRange("R_b",b[k],a[k]);
                                
                                
                                //computing areas in range [a;X_max]
                                N_s_1[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() * (BDTNormMC[k]->getValV());
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   ((BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal()) / (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()));
                                
                                //computing areas in range [b;a]
                                N_s_2[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() * (BDTNormMC[k]->getValV());
                                N_b_2[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   ((BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal()) / (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()));
                                
                                //cout<< "k is: "<< k <<endl;
                                
                                /*
                                cout<< "Count a small: "<< bdt_data[k]->sumEntries("1", "R_a") << " Count a large: "<< bdt_data[k]->sumEntries("1", "R_Large") <<endl;
                                cout<< "Count b small: "<< bdt_data[k]->sumEntries("1", "R_b") << " Count b large: "<< bdt_data[k]->sumEntries("1", "R_Large") <<endl;
                                
                                cout<< "Integral a small: "<< (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal())*BDTNorm[k]->getValV() << " Integral a large: "<< (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal())*BDTNorm[k]->getValV() <<endl;
                                cout<< "Integral b small: "<< (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal())*BDTNorm[k]->getValV() << " Integral b large: "<< (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal())*BDTNorm[k]->getValV() <<endl;
                                */
                                
                                /*
                                N_s_1[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() * (BDTNormMC[k]->getValV()) * 4500.0/59.0;
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * 4500.0/59.0;
                                N_s_2[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() * (BDTNormMC[k]->getValV()) * 4500.0/59.0;
                                N_b_2[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * 4500.0/59.0;
                                */
                                
                                
                                if(N_b_1[k]>0.0&&N_b_2[k]>0.0){
                                        //S1[k] = N_s_1[k] / sqrt(N_s_1[k] + N_b_1[k]);
                                        //S2[k] = N_s_2[k] / sqrt(N_s_2[k] + N_b_2[k]);
                                        
                                        S1[k] = sqrt( (2 * (N_s_1[k] + N_b_1[k]) * log(1 + (N_s_1[k]/N_b_1[k])) ) - 2 * N_s_1[k] );
                                        S2[k] = sqrt( (2 * (N_s_2[k] + N_b_2[k]) * log(1 + (N_s_2[k]/N_b_2[k])) ) - 2 * N_s_2[k] );
                                        
                                        //Combined significance
                                        S[k] = sqrt(S1[k]*S1[k] + S2[k]*S2[k]);
                                        
                                        if(!isnan(S[k])){
                                        
                                                //cout<<"N_s_1[k]: "<<N_s_1[k]<<" N_b_1[k]: "<<N_b_1[k]<<" N_s_2[k]: "<<N_s_2[k]<<" N_b_2[k]: "<<N_b_2[k]<<endl;
                                                //cout<<"S: "<<S[k]<<endl;
                                                
                                                S1_list[k].push_back(S1[k]);
                                                S2_list[k].push_back(S2[k]);
                                                a_list[k].push_back(a[k]);
                                                b_list[k].push_back(b[k]);
                                                S_list[k].push_back(S[k]);
                                                
                                                
                                                N_s_1_yield_list[k].push_back(N_s_1[k]);
                                                N_s_2_yield_list[k].push_back(N_s_2[k]);
                                                N_b_1_yield_list[k].push_back(N_b_1[k]);
                                                N_b_2_yield_list[k].push_back(N_b_2[k]);
                                                
                                                dim++;
                                                tau_cut1_vs_cut2_vs_sig[k]->Fill(a[k]+0.00000001,b[k]+0.00000001,S[k] );
                                        
                                        }
                                }
                        
                        }
                        
                }
        }
        
    //Taking absolute maximum of the combined significance
    
    double S_max[3];
    int S_maxIndex[3];
    float a_max[3];
    float b_max[3];
    double * a_array[3];
    double * sig_a_array[3];
    
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      if(S_list[i].size()>0){
              S_max[i] = *max_element(S_list[i].begin(), S_list[i].end());
              S_maxIndex[i] = std::max_element(S_list[i].begin(),S_list[i].end()) - S_list[i].begin();
              cout<<"S_max[i]: "<<S_max[i]<<" S_maxIndex[i]: "<<S_maxIndex[i]<<endl;
              
              a_max[i] = a_list[i].at(S_maxIndex[i]);
              b_max[i] = b_list[i].at(S_maxIndex[i]);
              
              cout<<"b cut: "<<b_max[i]<<" a cut: "<<a_max[i]<<endl;
              cout<<"a signal yield: "<< N_s_1_yield_list[i].at(S_maxIndex[i]) <<" a bkg yield: "<< N_b_1_yield_list[i].at(S_maxIndex[i]) <<endl;
              cout<<"b,a signal yield: "<< N_s_2_yield_list[i].at(S_maxIndex[i]) <<" b,a bkg yield: "<< N_b_2_yield_list[i].at(S_maxIndex[i]) <<endl;
              
              a_array[i] = &a_list[i][0];
              sig_a_array[i] = &S1_list[i][0];
      }
    }
    
    
    /*
    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 1800, 600);
    canvas4->Divide(3, 1);
    
    TGraph * Sig_Plot[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas4->cd( i+1 );
      Sig_Plot[i] = new TGraph(a_list[i].size(),a_array[i],sig_a_array[i]);
      Sig_Plot[i]->SetTitle("Significance vs. BDT Cut Value, #tau_"+cat_label[i]+";BDT Cut Value;Significance");
      Sig_Plot[i]->Draw();
      
    }
    */
    /*
    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 1800, 600);
    canvas5->Divide(3, 1);
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas5->cd( i+1 );
      tau_cut1_vs_cut2_vs_sig[i]->SetTitle("Significance vs. BDT Cut Values, #tau_"+cat_label[i]+";a;b");
      tau_cut1_vs_cut2_vs_sig[i]->SetStats(0);
      tau_cut1_vs_cut2_vs_sig[i]->Draw("colz");
      
    }
    */
    
    
    
    
    
    /*
    //Testing Symmetry of triplet mass
    double Loose_BDT_min = -0.9;
    double Loose_BDT_max = 0.3;
    Int_t Loose_BDT_N = 12; double Loose_BDT_step = (Loose_BDT_max - Loose_BDT_min)/Loose_BDT_N;
    double Loose_BDT_val(0.0);//x - val
    double Ratio_val(0.0);//y - val
    double binError(0.0);
    std::vector<double> Loose_BDT_list[3], Ratio_list[3];
    
    
    TH1F *myh = new TH1F("myh","T3M Distribution", 40, signal_region_min, signal_region_max);// used to temporarily store histograms
    for(int i=0; i<3; i++){
            hname=to_string(i+1);
            tau_loose_cut_vs_ratio[i] = new TH1D("tau_loose_cut_vs_ratio","tau_loose_cut_vs_ratio_"+hname,Loose_BDT_N,Loose_BDT_min,Loose_BDT_max);
            tau_loose_cut_vs_ratio[i]->Sumw2();
            
            tau_loose_cut_vs_ratio_num[i] = new TH1D("tau_loose_cut_vs_ratio_num","tau_loose_cut_vs_ratio_num_"+hname,Loose_BDT_N,Loose_BDT_min,Loose_BDT_max);
            tau_loose_cut_vs_ratio_den[i] = new TH1D("tau_loose_cut_vs_ratio_den","tau_loose_cut_vs_ratio_den_"+hname,Loose_BDT_N,Loose_BDT_min,Loose_BDT_max);
            tau_loose_cut_vs_ratio_num[i]->Sumw2();
            tau_loose_cut_vs_ratio_den[i]->Sumw2();
            
    }
    
    
    for(int i=0; i<Loose_BDT_N; i++){
            Loose_BDT_val = Loose_BDT_min + i * Loose_BDT_step + Loose_BDT_step/2.0;
            for(int k=0; k<3; k++){
                    
                    //tree[k]->Draw("tripletMass>>myh", const_cast<char*>(("tripletMass >= "+std::to_string(signal_region_min)+" && tripletMass <= "+std::to_string(signal_region_max)+" && isMC==0 && (tripletMass<="+std::to_string(signal_peak_region_min)+" || tripletMass>="+std::to_string(signal_peak_region_max)+") && bdt_cv>"+std::to_string(Loose_BDT_val)).c_str()),"goff");// This is useless. Used for testing validity of LooseBDT Cuts.
                    
                    
                    
                    
                    tree[k]->Draw("tripletMass>>myh", const_cast<char*>(("tripletMass >= "+std::to_string(signal_region_min)+" && tripletMass <= "+std::to_string(signal_region_max)+" && isMC==0 && (tripletMass<="+std::to_string(signal_peak_region_min)+" || tripletMass>="+std::to_string(signal_peak_region_max)+") && bdt_cv>"+std::to_string(Loose_BDT_val)).c_str()),"goff");
                    
                    //tree[k]->Draw("tripletMass>>myh", const_cast<char*>(("tripletMass >= "+std::to_string(signal_region_min)+" && tripletMass <= "+std::to_string(signal_region_max)+" && isMC==0 && (tripletMass<="+std::to_string(signal_peak_region_min)+" || tripletMass>="+std::to_string(signal_peak_region_max)+") && bdt_cv>"+std::to_string(Loose_BDT_val-0.05)+" && bdt_cv<"+std::to_string(Loose_BDT_val+0.05)).c_str()),"goff");
                    
                    data[k]->reset();
                    data[k] = new RooDataHist("data"+hname, "data", *InvMass[k], Import(*myh));
                    
                    //BDTOutput_x[k]->setRange("R_Large",Loose_BDT_val,100);
                    
                    //Ratio_val = data[k]->sumEntries("1", "R1")  /   data[k]->sumEntries("1", "R2") * ( (signal_region_max-signal_peak_region_max) / (signal_peak_region_min-signal_region_min) );
                    
                    //Ratio_val = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  data[k]->sumEntries("1", "R1,R2")  /   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal());
                    
                    
                    //Using exponential
                    exp_fitresult[k] = exp_pdf[k]->fitTo(*data[k], Range("R1,R2"), Save());
                    Ratio_val = exp_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getVal() / exp_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R1,R2"))->getVal();
                    //binError = exp_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getPropagatedError(*exp_fitresult[k]) * ExpNorm[k]->getValV();
                    binError = exp_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getPropagatedError(*exp_fitresult[k]);
                    
                    
                    
                    
                    if(Ratio_val>0.0){
                            Loose_BDT_list[k].push_back(Loose_BDT_val);
                            Ratio_list[k].push_back(Ratio_val);
                            //tau_loose_cut_vs_ratio_num[k]->SetBinContent( tau_loose_cut_vs_ratio_num[k]->FindBin(Loose_BDT_val+0.000001), exp_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getVal() );
                            //tau_loose_cut_vs_ratio_den[k]->SetBinContent( tau_loose_cut_vs_ratio_den[k]->FindBin(Loose_BDT_val+0.000001), exp_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R1,R2"))->getVal() );
                            //tau_loose_cut_vs_ratio_num[k]->Divide( tau_loose_cut_vs_ratio_den[k]);
                            tau_loose_cut_vs_ratio[k]->SetBinContent( tau_loose_cut_vs_ratio[k]->FindBin(Loose_BDT_val), Ratio_val );
                            tau_loose_cut_vs_ratio[k]->SetBinError( tau_loose_cut_vs_ratio[k]->FindBin(Loose_BDT_val), binError );
                            //tau_loose_cut_vs_ratio[k]->Fill( Loose_BDT_val+0.000001, Ratio_val );
                    }
            }
    }
    
    TCanvas *canvas6 = new TCanvas("canvas6", "canvas6", 1800, 600);
    canvas6->Divide(3, 1);
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas6->cd( i+1 );
      //tau_loose_cut_vs_ratio[i]->SetTitle("Ratio vs. BDT Cut Values, #tau_"+cat_label[i]+";BDT Cut;Ratio");
      tau_loose_cut_vs_ratio[i]->SetTitle("SB Symmetry, #tau_"+cat_label[i]+";BDT Cut;Ratio");
      tau_loose_cut_vs_ratio[i]->SetStats(0);
      tau_loose_cut_vs_ratio[i]->SetMaximum(0.2);
      tau_loose_cut_vs_ratio[i]->SetMinimum(0.1);
      tau_loose_cut_vs_ratio[i]->Draw("E1");
      
    }
    */
        
    
    
    
    
    
}