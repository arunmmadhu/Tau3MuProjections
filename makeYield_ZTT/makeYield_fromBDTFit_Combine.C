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
using std::ofstream;

void makeYield_fromBDTFit_Combine () 
{
    system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv");
    
    //specify the luminosities here
    //double lumi_values[] = {97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0};
    double lumi_values[] = {59.0};
    int lumi_size = sizeof(lumi_values)/ sizeof(double);
    
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
    
    TString card_modifier_name[3];
    card_modifier_name[0] = "ZTT_tauh_test";
    card_modifier_name[1] = "ZTT_taumu_test";
    card_modifier_name[2] = "ZTT_taue_test";
    
    TString combined_card_name[3];
    combined_card_name[0] = "Cat_1_Mod";
    combined_card_name[1] = "Cat_2_Mod";
    combined_card_name[2] = "Cat_3_Mod";
    
    TString hname;
    
    float signal_region_min(1.6);
    float signal_region_max(2.0);
    
    //float signal_peak_region_min(1.73);
    //float signal_peak_region_max(1.82);
    
    float signal_peak_region_min[3];
    float signal_peak_region_max[3];
    
    float BDT_Score_Min(-0.9);
    
    int triplet_mass_bins(40);
    
    float Loose_BDT_Cut(-0.4);
    
    float Gaussian_Sigma_From_Loose_BDT_Cut[3];
    Gaussian_Sigma_From_Loose_BDT_Cut[0] = 0.0165883;
    Gaussian_Sigma_From_Loose_BDT_Cut[1] = 0.0169745;
    Gaussian_Sigma_From_Loose_BDT_Cut[2] = 0.0170951;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      signal_peak_region_min[i]=1.77686-2.3*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      signal_peak_region_max[i]=1.77686+2.3*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
      
      int val_1 = (int)(((signal_peak_region_min[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      int val_2 = (int)(((signal_peak_region_max[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      
      signal_peak_region_min[i]=signal_region_min + val_1 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      signal_peak_region_max[i]=signal_region_min + val_2 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
    }
    
    //Filename and histograms
    TFile * file_tau[3];
    
    TH1D  * tau_T3Mu[3];
    TH1D  * tau_T3Mu_Dat[3];
    TH1D  * tau_BDT_Output_Data[3];
    TH1D  * tau_BDT_Output_MC[3];
    TH2D  * tau_T3Mu_vs_BDT_Output_Data[3];
    TH1D  * tau_T3Mu_vs_BDT_Output_Data_Projection[3];
    
    TH2D  * tau_cut1_vs_cut2_vs_sig[lumi_size][3];
    TH1D  * bdt_cut_vs_limit[lumi_size][3];
    
    
    TFile *TreeFile = new TFile("Combine_Tree_ztau3mutau.root","READ");
    TTree *tree[3];
    
    Float_t tripletMass;
    Float_t bdt_cv;
    Float_t weight;
    Float_t isMC;
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      tree[i] = (TTree *) TreeFile->Get(cat_base[i]);
      
      tau_T3Mu[i] = new TH1D("tau_T3Mu","tau_T3Mu_"+hname,40,signal_region_min,signal_region_max);
      tau_T3Mu_Dat[i] = new TH1D("tau_T3Mu_Dat","tau_T3Mu_Dat_"+hname,40,signal_region_min,signal_region_max);
      tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,100,BDT_Score_Min,0.9);
      tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,100,BDT_Score_Min,0.9);
      
      tree[i]->SetBranchAddress("tripletMass",&tripletMass);
      tree[i]->SetBranchAddress("bdt_cv",&bdt_cv);
      tree[i]->SetBranchAddress("weight",&weight);
      tree[i]->SetBranchAddress("isMC",&isMC);
      
      Long64_t nentries = tree[i]->GetEntries();
      for (Long64_t j=0;j<nentries;j++) {
        tree[i]->GetEntry(j);
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max){
                if(isMC>0){
                  if(bdt_cv>Loose_BDT_Cut){
                    tau_T3Mu[i]->Fill(tripletMass,weight);
                  }
                  tau_BDT_Output_MC[i]->Fill(bdt_cv,weight);
                }
                if(isMC==0 && (tripletMass<=signal_peak_region_min[i] || tripletMass>=signal_peak_region_max[i]) ){//blinded
                  if(bdt_cv>Loose_BDT_Cut){
                    tau_T3Mu_Dat[i]->Fill(tripletMass);
                  }
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      //tau_T3Mu_vs_BDT_Output_Data[i] = (TH2D*)file_tau[i]->Get(cat_name[i]+"BDT_2Dscan_TripletMassData");
      
      
    }
    
    
    
    
    
    
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
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, #tau_"+cat_label[i],BDT_Score_Min,0.9);
      BDTOutput_x[i]->setRange("BDT_Range",BDT_Score_Min,0.9);
      bgausmean[i] = new RooRealVar("bgausmean"+hname, "bgausmean" , -0.5, BDT_Score_Min,0.0) ;
      bgaussigma_a[i] = new RooRealVar("bgaussigma_a"+hname, "bgaussigma_a" , 0.2, 0.000001, 1.0) ;
      bgaussigma_b[i] = new RooRealVar("bgaussigma_b"+hname, "bgaussigma_b" , 0.2, 0.000001, 1.0);
      
      bgaus_dist[i] = new RooBifurGauss("bgaus_dist"+hname, "bgaus dist", *BDTOutput_x[i], *bgausmean[i], *bgaussigma_a[i], *bgaussigma_b[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*bgaus_dist[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    
    
    
    /*
    //BDT Fit Plots
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1800, 600);
    canvas2->Divide(3, 1);
    
    RooPlot * xFrame_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas2->cd( i+1 );
      xFrame_bdt[i] = BDTOutput_x[i]->frame();
      bdt_data[i]->plotOn(xFrame_bdt[i]);
      BDT_distribution[i]->plotOn(xFrame_bdt[i],LineColor(4),LineWidth(2), Normalization(bdt_data[i]->sumEntries("1", "BDT_Range"), RooAbsReal::NumEvent),ProjectionRange("BDT_Range"));
      xFrame_bdt[i]->SetTitle("BDT Output, #tau_"+cat_label[i]);
      xFrame_bdt[i]->SetXTitle("BDT Score");
      xFrame_bdt[i]->SetYTitle("Events");
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
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas3->cd( i+1 );
      xFrame_bdtMC[i] = BDTOutput_x[i]->frame();
      bdt_MC[i]->plotOn(xFrame_bdtMC[i]);
      BDT_distributionMC[i]->plotOn(xFrame_bdtMC[i],LineColor(4),LineWidth(2), Normalization(bdt_MC[i]->sumEntries("1", "BDT_Range"), RooAbsReal::NumEvent),ProjectionRange("BDT_Range"));
      xFrame_bdtMC[i]->SetTitle("BDT Output, #tau_"+cat_label[i]);
      xFrame_bdtMC[i]->SetXTitle("BDT Score");
      xFrame_bdtMC[i]->SetYTitle("Events");
      xFrame_bdtMC[i]->Draw();
    }
    */
    
    
    
    
    
    
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
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      InvMass[i] = new RooRealVar("InvMass"+hname,"InvMass, #tau_"+cat_label[i],signal_region_min,signal_region_max);
      InvMass[i]->setRange("R1",signal_region_min,signal_peak_region_min[i]); //background   
      InvMass[i]->setRange("R2",signal_peak_region_max[i],signal_region_max); //background
      InvMass[i]->setRange("R3",1.73,1.82); //signal range for fitting
      InvMass[i]->setRange("R4",signal_peak_region_min[i],signal_peak_region_max[i]); //signal range for yield
      
      
      //Flat fit for data
      poly[i] = new RooPolynomial("poly"+hname, "poly dist", *InvMass[i]);
      data[i] = new RooDataHist("data"+hname, "data", *InvMass[i], Import(*tau_T3Mu_Dat[i]));
      LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,200.0);
      pdf[i] = new RooAddPdf("pdf"+hname, "pdf", RooArgList(*poly[i]), RooArgList(*LineNorm[i]));
      fitresult[i] = pdf[i]->fitTo(*data[i], Range("R1,R2"), Save());
      
      
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
      
      cout << "For category, " << cat_label[i] << ", sigma is: " << sigma[i]->getValV() << endl;
      
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
    /*
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 1800, 600);
    canvas1->Divide(3, 1);
    
    RooPlot * xFrame[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas1->cd( i+1 );
      xFrame[i] = InvMass[i]->frame();
      data[i]->plotOn(xFrame[i]);
      pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
      mc[i]->plotOn(xFrame[i]);
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), #tau_"+cat_label[i]);
      xFrame[i]->SetXTitle("3#mu inv. mass (GeV)");
      xFrame[i]->SetYTitle("Events");
      xFrame[i]->Draw();
    }
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
    
    
    
    
    
    
        //Calculate limits
        double a[3],b[3];// BDT score cut
        double N_s_1[3], N_s_2[3];
        double N_b_1[3], N_b_2[3];
        double S1[3], S2[3], S[3];
        std::vector<double> S1_list[3], S2_list[3], S_list[lumi_size][3], a_list[lumi_size][3], b_list[lumi_size][3], N_s_1_yield_list[lumi_size][3], N_s_2_yield_list[lumi_size][3], N_b_1_yield_list[lumi_size][3], N_b_2_yield_list[lumi_size][3];
        
        TString command_a[3];
        TString command_b[3];
        TString command_a_and_b[3];
        TString command_run[3];
        TString command_copy[3];
        TString command_recreate_workdir;
        TString command_recreate_lumidir;
        TString command_recreate_lumi_scan_dir;
        TString command_copy_datacard[3];
        
        
        
        
        double Xa_min[3];
        double Xa_max[3];
        
        Xa_min[0] = 0.34;
        Xa_min[1] = 0.36;
        Xa_min[2] = 0.30;
        Xa_max[0] = 0.52;
        Xa_max[1] = 0.54;
        Xa_max[2] = 0.48;
        
        //Loop on both cuts in [X_min;X_max]
        Int_t dim = 0;
        //Increase N to increase (a,b) scan granularity!
        Int_t N_a = 18; //Int_t N_b = 5;//ideal: N = 16
        double step_a[3];
        double step_b[3];
        
        
        command_recreate_workdir = "rm -r workdir/; mkdir workdir/";
        system(command_recreate_workdir);
        
    
    for(int iter_lumi=0; iter_lumi<lumi_size; iter_lumi++){//iterating over different values of lumis
        
        command_recreate_lumidir = "rm -r workdir/lumi_" + to_string((int)lumi_values[iter_lumi])+"; mkdir workdir/lumi_" + to_string((int)lumi_values[iter_lumi])+"";
        system(command_recreate_lumidir);
        
        for(int i=0; i<3; i++){// stores limit as function of bdt_cut
          hname=to_string(i+1);
          
          step_a[i]=(Xa_max[i] - Xa_min[i])/N_a;
          
          bdt_cut_vs_limit[iter_lumi][i] = new TH1D("bdt_cut_vs_limit","bdt_cut_vs_limit_"+hname,N_a,Xa_min[i],Xa_max[i]);
              
        }
        
        for(int i=0; i<N_a; i++){//For loop for bdt cuts in range [X_min;X_max]
                
                command_recreate_lumi_scan_dir = "rm -r workdir/lumi_" + to_string((int)lumi_values[iter_lumi])+"/scan_"+to_string(i+1)+"; mkdir workdir/lumi_" + to_string((int)lumi_values[iter_lumi])+"/scan_"+to_string(i+1)+"/";
                system(command_recreate_lumi_scan_dir);
                
                        for(int k=0; k<3; k++){// for tauh/taumu/taue or A/B/C
                                
                                a[k] = Xa_min[k] + i * step_a[k];// bdt cut
                                
                                cout<<"i: "<<i<<" lumi: "<<lumi_values[iter_lumi]<<endl;
                                
                                
                                //Calculate signal and background yields from BDT output plots
                                
                                BDTOutput_x[k]->setRange("R_a",a[k],100.0);
                                N_s_1[k] = bdt_MC[k]->sumEntries("1", "R_a") * lumi_values[iter_lumi]/59.0;
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * lumi_values[iter_lumi]/59.0;
                                
                                if( (std::round(N_b_1[k] / 0.00001) * 0.00001) > 0.0){
                                        
                                        //Create datacard from signal and bkg values
                                        command_a[k] = "python card_modifiers/"+ card_modifier_name[k] +".py --luminosity " + std::to_string(59.0) + " --s " + std::to_string(N_s_1[k]) + " --b " + std::to_string(N_b_1[k]) + " --cuttype a";
                                        system(command_a[k]);
                                        
                                        
                                        //Copy datacards to workdir
                                        command_copy_datacard[k] = "cp "+combined_card_name[k]+"_a.txt workdir/lumi_" + to_string((int)lumi_values[iter_lumi])+"/scan_"+to_string(i+1)+"/";
                                        system(command_copy_datacard[k]);
                                        
                                        
                                        //Whether to use HybridNew or AsymptoticLimits for calculating limits
                                        bool Whether_Hybrid(true);
                                        
                                        
                                        
                                        
                                        //HybridNew
                                        if(Whether_Hybrid){
                                        
                                        command_run[k] = "combine -M HybridNew "+combined_card_name[k]+"_a.txt --cl 0.9 -t -1  --expectedFromGrid=0.5 --rMin 0.0 --rMax 15.0  > out_mid_"+ to_string(k+1) +".txt";
                                        system(command_run[k]);
                                        
                                        std::ifstream f1("out_mid_"+ to_string(k+1) +".txt");
                                        std::string line;
                                        while (std::getline(f1, line)) {
                                          if (line.find("Limit: r <") != std::string::npos) {
                                            std::vector<std::string> linsp;
                                            std::istringstream iss(line);
                                            for (std::string s; iss >> s;) linsp.push_back(s);
                                            S[k] = std::stod(linsp[3]);
                                          }
                                        }
                                        
                                        }
                                        
                                        
                                        //AsymptoticLimits
                                        else{
                                        
                                        command_run[k] = "combine -M AsymptoticLimits "+combined_card_name[k]+"_a.txt --cl 0.9 -t -1  > out"+ to_string(k+1) +".txt";
                                        system(command_run[k]);
                                        
                                        std::ifstream f1("out"+ to_string(k+1) +".txt");
                                        std::string line;
                                        while (std::getline(f1, line)) {
                                          if (line.find("50.0%") != std::string::npos) {
                                            std::vector<std::string> linsp;
                                            std::istringstream iss(line);
                                            for (std::string s; iss >> s;) linsp.push_back(s);
                                            S[k] = std::stod(linsp.back());
                                          }
                                        }
                                        
                                        }
                                        
                                        
                                        
                                        
                                        if(!isnan(S[k])){
                                        
                                                a_list[iter_lumi][k].push_back(a[k]);
                                                S_list[iter_lumi][k].push_back(S[k]);
                                                
                                                
                                                N_s_1_yield_list[iter_lumi][k].push_back(N_s_1[k]);
                                                N_b_1_yield_list[iter_lumi][k].push_back(N_b_1[k]);
                                                
                                                dim++;
                                                bdt_cut_vs_limit[iter_lumi][k]->Fill(a[k]+0.00000001,S[k] );
                                        
                                        }//isnan check
                                }//bkg non-zero check
                        
                        }//end k
                        
        }//end i
    }//end iter_lumi
    
    
    
    //Taking absolute minimum of the limits
    
    
    double S_max[3];
    int S_maxIndex[3];
    float a_max[3];
    float b_max[3];
    double * a_array[3];
    double * sig_a_array[3];
    
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      S_max[i] = *min_element(S_list[0][i].begin(), S_list[0][i].end());
      S_maxIndex[i] = std::min_element(S_list[0][i].begin(),S_list[0][i].end()) - S_list[0][i].begin();
      cout<<"S_min[i]: "<<S_max[i]<<" S_minIndex[i]: "<<S_maxIndex[i]<<endl;
      
      a_max[i] = a_list[0][i].at(S_maxIndex[i]);
      
      cout<<" a cut: "<<a_max[i]<<endl;
      cout<<"a signal yield: "<< N_s_1_yield_list[0][i].at(S_maxIndex[i]) <<" a bkg yield: "<< N_b_1_yield_list[0][i].at(S_maxIndex[i]) <<endl;
      
    }
    
    
    //Store lumis, sig/bkg yield in a text file. If this doesn't work, run line #33 first
    
    ofstream myfile;
    myfile.open("ZTT_Lumi_Limit.txt", ios::out | ios::trunc);
    
    
    for(int iter_lumi=0; iter_lumi<lumi_size; iter_lumi++){
            for(int i=0; i<3; i++){
              hname=to_string(i+1);
              
              double min_limit_val = *min_element(S_list[iter_lumi][i].begin(), S_list[iter_lumi][i].end());
              int min_limit_idx = std::min_element(S_list[iter_lumi][i].begin(),S_list[iter_lumi][i].end()) - S_list[iter_lumi][i].begin();
              
              myfile << "Luminosity "<< lumi_values[iter_lumi]<<" limit "<< min_limit_val <<" category "<< "Cat_"+print_label[i] <<" a_cut "<<a_list[iter_lumi][i].at(min_limit_idx)<<" a_signal_yield "<< N_s_1_yield_list[iter_lumi][i].at(min_limit_idx) <<" a_bkg_yield "<< N_b_1_yield_list[iter_lumi][i].at(min_limit_idx) << "\n";
              
            }
    }
    myfile.close();
    
    
    
    
    // Plot of limit as function of bdt_cut
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
      //bdt_cut_vs_limit[i]->SetTitle("Limit vs. BDT Cut Values, #tau_"+cat_label[i]+";a;b");
      bdt_cut_vs_limit[i]->SetTitle("Limit vs. BDT Cut Values, #tau_"+cat_label[i]+";a;b");
      bdt_cut_vs_limit[i]->SetStats(0);
      bdt_cut_vs_limit[i]->SetMinimum(S_max[i]);
      bdt_cut_vs_limit[i]->Draw("colz");
      
    }
    */
    
    
}