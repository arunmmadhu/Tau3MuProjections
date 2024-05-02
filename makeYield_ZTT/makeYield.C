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
#include "TH1F.h"
#include "TH2F.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>

using namespace RooFit;

void makeYield () 
{
    system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/ZTT/CMSSW_11_3_4/src/; cmsenv; cd -");
    
    //Category names: 0=tauh_A, 1=tauh_B, 2=taumu, 3=taue
    
    TString cat_base[4];
    cat_base[0] = "ztau3mutauh_A";
    cat_base[1] = "ztau3mutauh_B";
    cat_base[2] = "ztau3mutaumu";
    cat_base[3] = "ztau3mutaue";
    
    TString cat_label[4];
    cat_label[0] = "{h,A}";
    cat_label[1] = "{h,B}";
    cat_label[2] = "{#mu}";
    cat_label[3] = "{e}";
    
    TString print_label[4];
    print_label[0] = "h_A";
    print_label[1] = "h_B";
    print_label[2] = "mu";
    print_label[3] = "e";
    
    TString hname;
    
    //BDT Cuts a
    float BDT_Cut_a[4];
    
    //Cuts from limits
    
    BDT_Cut_a[0] = 0.345;
    BDT_Cut_a[1] = 0.12;
    BDT_Cut_a[2] = 0.35;
    BDT_Cut_a[3] = 0.175;
    
    
    //Cuts from significance
    
    BDT_Cut_a[0] = 0.375;
    BDT_Cut_a[1] = 0.18;
    BDT_Cut_a[2] = 0.44;
    BDT_Cut_a[3] = 0.1925;
    
    
    float signal_region_min(1.6);
    float signal_region_max(2.0);
    
    float signal_peak_region_min[4];
    float signal_peak_region_max[4];
    
    float BDT_Score_Min(-0.9);
    
    int triplet_mass_bins(40);
    
    float Gaussian_Sigma_From_Loose_BDT_Cut[4];
    Gaussian_Sigma_From_Loose_BDT_Cut[0] = 0.0170238;
    Gaussian_Sigma_From_Loose_BDT_Cut[1] = 0.0167713;
    Gaussian_Sigma_From_Loose_BDT_Cut[2] = 0.0172921;
    Gaussian_Sigma_From_Loose_BDT_Cut[3] = 0.0172922;
    
    //sigma_scale from training
    
    float sigma_scale[4];
    
    sigma_scale[0] = 2.0;
    sigma_scale[1] = 2.0;
    sigma_scale[2] = 2.0;
    sigma_scale[3] = 2.8;
    
    for(int i=0; i<4; i++){
      hname=to_string(i+1);
      
      signal_peak_region_min[i]=1.77686-sigma_scale[i]*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      signal_peak_region_max[i]=1.77686+sigma_scale[i]*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
      
      int val_1 = (int)(((signal_peak_region_min[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      int val_2 = (int)(((signal_peak_region_max[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      
      signal_peak_region_min[i]=signal_region_min + val_1 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      signal_peak_region_max[i]=signal_region_min + val_2 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
    }
    
    //Filename and histograms
    TFile * file_tau[4];
    
    TH1D  * tau_T3Mu[4];
    TH1D  * tau_T3Mu_Dat[4];
    TH1D  * tau_BDT_Output_Data[4];
    TH1D  * tau_BDT_Output_MC[4];
    TH2D  * tau_T3Mu_vs_BDT_Output_Data[4];
    TH1D  * tau_T3Mu_vs_BDT_Output_Data_Projection[4];
    
    TH2D  * tau_cut1_vs_cut2_vs_sig[4];
    TH2D  * tau_cut1_vs_cut2_vs_limit[4];
    
    TFile *TreeFile = new TFile("Combine_Tree_ztau3mutau.root","READ");
    TTree *tree[4];
    
    Float_t tripletMass;
    Float_t bdt_cv;
    Float_t weight;
    Float_t isMC;
    for(int i=0; i<4; i++){
      hname=to_string(i+1);
      
      tree[i] = (TTree *) TreeFile->Get(cat_base[i]);
      
      tau_T3Mu[i] = new TH1D("tau_T3Mu","tau_T3Mu_"+hname,40,signal_region_min,signal_region_max);
      tau_T3Mu_Dat[i] = new TH1D("tau_T3Mu_Dat","tau_T3Mu_Dat_"+hname,40,signal_region_min,signal_region_max);
      tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,100,BDT_Score_Min,1.0);
      tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,100,BDT_Score_Min,1.0);
      
      tree[i]->SetBranchAddress("tripletMass",&tripletMass);
      tree[i]->SetBranchAddress("bdt_cv",&bdt_cv);
      tree[i]->SetBranchAddress("weight",&weight);
      tree[i]->SetBranchAddress("isMC",&isMC);
      
      Long64_t nentries = tree[i]->GetEntries();
      for (Long64_t j=0;j<nentries;j++) {
        tree[i]->GetEntry(j);
        
        if(tripletMass>=signal_region_min && tripletMass<=signal_region_max && bdt_cv>=BDT_Cut_a[i]){
        //if(tripletMass>=signal_region_min && tripletMass<=signal_region_max && bdt_cv >= BDT_Cut_b[i] && bdt_cv < BDT_Cut_a[i]){
                if(isMC>0){
                  tau_T3Mu[i]->Fill(tripletMass,weight);
                  tau_BDT_Output_MC[i]->Fill(bdt_cv,weight);
                }
                if(isMC==0 && (tripletMass<=signal_peak_region_min[i] || tripletMass>=signal_peak_region_max[i]) ){//blinded
                  tau_T3Mu_Dat[i]->Fill(tripletMass);
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      
      
    }
    
    
    
    
    
    
    //Triplet Mass Fits
    RooRealVar * InvMass[4];
    
    RooPolynomial * poly[4];
    RooDataHist * data[4];
    RooRealVar * LineNorm[4];
    RooAddPdf * pdf[4];
    RooFitResult * fitresult[4];
    
    RooRealVar * mean[4];
    RooRealVar * sigma[4];
    RooGaussian * Gauss[4];
    RooDataHist * mc[4];
    RooRealVar * GaussNorm[4];
    RooAddPdf * mc_pdf[4];
    RooFitResult * mc_fitresult[4];
    
    RooExponential * Expo[4];
    RooRealVar * lambda[4];
    RooRealVar * ExpNorm[4];
    RooAddPdf * exp_pdf[4];
    RooFitResult * exp_fitresult[4];
    
    for(int i=0; i<4; i++){
      hname=to_string(i+1);
      
      InvMass[i] = new RooRealVar("InvMass"+hname,"InvMass, #tau_"+cat_label[i],signal_region_min,signal_region_max);
      InvMass[i]->setRange("R1",signal_region_min,signal_peak_region_min[i]); //background   
      InvMass[i]->setRange("R2",signal_peak_region_max[i],signal_region_max); //background
      InvMass[i]->setRange("R3",1.73,1.82); //signal range for fitting
      InvMass[i]->setRange("R4",signal_peak_region_min[i],signal_peak_region_max[i]); //signal range for yield
      
      
      //Flat fit for data
      poly[i] = new RooPolynomial("poly"+hname, "poly dist", *InvMass[i]);
      data[i] = new RooDataHist("data"+hname, "data", *InvMass[i], Import(*tau_T3Mu_Dat[i]));
      LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,20000);
      pdf[i] = new RooAddPdf("pdf"+hname, "pdf", RooArgList(*poly[i]), RooArgList(*LineNorm[i]));
      fitresult[i] = pdf[i]->fitTo(*data[i], Range("R1,R2"), Save());
      
      
      //Testing Exponential Fit
      lambda[i] = new RooRealVar("lambda"+hname,"lambda",0.1, -1.5, 1.5);
      Expo[i] = new RooExponential("Expo"+hname, "Exponential PDF", *InvMass[i],  *lambda[i]);
      ExpNorm[i] = new RooRealVar("ExpNorm"+hname, "ExpNorm",  2.0,0.001,10000);
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
    
    
    double pdf_integral_restricted[4];
    double mc_pdf_integral_restricted[4];
    double nData[4];
    double nSignal[4];
    double nSignal_restricted[4];
    double scaling[4];
    
    for(int i=0; i<4; i++){
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
    
    
    
    //Triplet Mass Fit Plots
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 2400, 600);
    canvas1->Divide(4, 1);
    
    RooPlot * xFrame[4];
    
    for(int i=0; i<4; i++){
      hname=to_string(i+1);
      
      canvas1->cd( i+1 );
      xFrame[i] = InvMass[i]->frame();
      data[i]->plotOn(xFrame[i]);
      pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
      //exp_pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));//plots exponential fit
      //mc[i]->plotOn(xFrame[i]);
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), #tau_"+cat_label[i]);
      xFrame[i]->SetXTitle("3#mu inv. mass (GeV)");
      xFrame[i]->SetYTitle("Events");
      xFrame[i]->Draw();
    }
    
    
    
    
}