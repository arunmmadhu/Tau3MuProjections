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
    
    //BDT Cuts b
    float BDT_Cut_b[3];
    BDT_Cut_b[0] = 0.275;
    BDT_Cut_b[1] = 0.325;
    BDT_Cut_b[2] = 0.25;
    
    //BDT Cuts a
    float BDT_Cut_a[3];
    BDT_Cut_a[0] = 0.425;
    BDT_Cut_a[1] = 0.45;
    BDT_Cut_a[2] = 0.4;
    
    float signal_region_min(1.4);
    float signal_region_max(2.1);
    
    float signal_peak_region_min(1.73);
    float signal_peak_region_max(1.82);
    
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
    
    TFile *TreeFile = new TFile("Combine_Tree_ztau3mutau.root","READ");
    TTree *tree[3];
    
    Float_t tripletMass;
    Float_t bdt_cv;
    Float_t weight;
    Float_t isMC;
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      tree[i] = (TTree *) TreeFile->Get(cat_base[i]);
      
      tau_T3Mu[i] = new TH1D("tau_T3Mu","tau_T3Mu_"+hname,70,signal_region_min,signal_region_max);
      tau_T3Mu_Dat[i] = new TH1D("tau_T3Mu_Dat","tau_T3Mu_Dat_"+hname,70,signal_region_min,signal_region_max);
      tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,100,-0.9,0.9);
      tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,100,-0.9,0.9);
      
      tree[i]->SetBranchAddress("tripletMass",&tripletMass);
      tree[i]->SetBranchAddress("bdt_cv",&bdt_cv);
      tree[i]->SetBranchAddress("weight",&weight);
      tree[i]->SetBranchAddress("isMC",&isMC);
      
      Long64_t nentries = tree[i]->GetEntries();
      for (Long64_t j=0;j<nentries;j++) {
        tree[i]->GetEntry(j);
        
        //if(tripletMass>=signal_region_min && tripletMass<=signal_region_max && bdt_cv>=BDT_Cut_a[i]){
        if(tripletMass>=signal_region_min && tripletMass<=signal_region_max && bdt_cv >= BDT_Cut_b[i] && bdt_cv < BDT_Cut_a[i]){
                if(isMC>0){
                  tau_T3Mu[i]->Fill(tripletMass,weight);
                  tau_BDT_Output_MC[i]->Fill(bdt_cv,weight);
                }
                if(isMC==0 && (tripletMass<=signal_peak_region_min || tripletMass>=signal_peak_region_max) ){//blinded
                  tau_T3Mu_Dat[i]->Fill(tripletMass);
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      //tau_T3Mu_vs_BDT_Output_Data[i] = (TH2D*)file_tau[i]->Get(cat_name[i]+"BDT_2Dscan_TripletMassData");
      
      
    }
    
    
    
    
    
    
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
      InvMass[i]->setRange("R3",1.70,1.85); //signal range for fitting
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
      //mc[i]->plotOn(xFrame[i]);
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), #tau_"+cat_label[i]);
      xFrame[i]->SetXTitle("3#mu inv. mass (GeV)");
      xFrame[i]->SetYTitle("Events");
      xFrame[i]->Draw();
    }
    
    
    
    
}