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

void makeYield_fromBDTFit_Combine_W () 
{
    system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv");
    
    //specify the luminosities here
    double lumi_values[] = {97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0};
    //double lumi_values[] = {59.0};
    int lumi_size = sizeof(lumi_values)/ sizeof(double);
    
    //Category names: 0=tauh, 1=taumu, 2=taue
    TString cat_name[3];
    cat_name[0] = "ztau3mutauh_default_";
    cat_name[1] = "ztau3mutaumu_default_";
    cat_name[2] = "ztau3mutaue_default_";
    
    TString cat_base[3];
    cat_base[0] = "tree";
    cat_base[1] = "tree";
    cat_base[2] = "tree";
    
    TString cat_label[3];
    cat_label[0] = "A";
    cat_label[1] = "B";
    cat_label[2] = "C";
    
    TString print_label[3];
    print_label[0] = "A";
    print_label[1] = "B";
    print_label[2] = "C";
    
    TString card_modifier_name[3];
    card_modifier_name[0] = "ZTT_tauh_test";
    card_modifier_name[1] = "ZTT_taumu_test";
    card_modifier_name[2] = "ZTT_taue_test";
    
    TString combined_card_name[3];
    combined_card_name[0] = "ZTT_tauh_Combined_Mod";
    combined_card_name[1] = "ZTT_taumu_Combined_Mod";
    combined_card_name[2] = "ZTT_taue_Combined_Mod";
    
    TString hname;
    
    float signal_region_min(1.6);
    float signal_region_max(2.0);
    
    //float signal_peak_region_min(1.73);
    //float signal_peak_region_max(1.82);
    
    float signal_peak_region_min[3];
    float signal_peak_region_max[3];
    
    float BDT_Score_Min(0.98);
    
    int triplet_mass_bins(40);
    
    float Loose_BDT_Cut[3];
    Loose_BDT_Cut[0] = 0.98;
    Loose_BDT_Cut[1] = 0.98;
    Loose_BDT_Cut[2] = 0.98;
    
    float Gaussian_Sigma_From_Loose_BDT_Cut[3];
    Gaussian_Sigma_From_Loose_BDT_Cut[0] = 0.0115631;
    Gaussian_Sigma_From_Loose_BDT_Cut[1] = 0.016846;
    Gaussian_Sigma_From_Loose_BDT_Cut[2] = 0.0216494;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      signal_peak_region_min[i]=1.77686-2*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      signal_peak_region_max[i]=1.77686+2*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
      
      int val_1 = (int)(((signal_peak_region_min[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      int val_2 = (int)(((signal_peak_region_max[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      
      signal_peak_region_min[i]=signal_region_min + val_1 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      signal_peak_region_max[i]=signal_region_min + val_2 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
    }
    
    double MC_NORM17 = 30541./(500e+3)*(8580+11370)*0.1138/0.1063*1E-7;
    double MC_NORM18 = 59828./(492e+3)*(8580+11370)*0.1138/0.1063*1E-7;
    
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
    
    
    TFile *TreeFile_sig = new TFile("luca_root/signal_threeMedium_weighted_16Mar2022.root","READ");
    TFile *TreeFile_bkg = new TFile("luca_root/background_threeMedium-UNBLINDED.root","READ");
    TTree *tree_MC[3];
    TTree *tree_bkg[3];
    
    Double_t tripletMass;
    Double_t bdt_cv;
    Double_t mcweight;
    Double_t weight;
    Double_t isMC;
    Double_t cand_refit_tau_massE;
    Float_t year;
    Double_t tau_sv_ls;
    
    bool categ[3];//Categ A B or C
    
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      tree_MC[i] = (TTree *) TreeFile_sig->Get(cat_base[i]);
      tree_bkg[i] = (TTree *) TreeFile_bkg->Get(cat_base[i]);
      
      tau_T3Mu[i] = new TH1D("tau_T3Mu","tau_T3Mu_"+hname,40,signal_region_min,signal_region_max);
      tau_T3Mu_Dat[i] = new TH1D("tau_T3Mu_Dat","tau_T3Mu_Dat_"+hname,40,signal_region_min,signal_region_max);
      tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,160,BDT_Score_Min,1.0);
      tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,20,BDT_Score_Min,1.0);
      
      
      //MC
      tree_MC[i]->SetBranchAddress("cand_refit_tau_mass",&tripletMass);
      tree_MC[i]->SetBranchAddress("cand_refit_tau_massE",&cand_refit_tau_massE);
      tree_MC[i]->SetBranchAddress("bdt",&bdt_cv);
      tree_MC[i]->SetBranchAddress("mcweight",&mcweight);
      tree_MC[i]->SetBranchAddress("weight",&weight);
      tree_MC[i]->SetBranchAddress("year",&year);
      tree_MC[i]->SetBranchAddress("tau_sv_ls",&tau_sv_ls);
      
      
      
      Long64_t nentries1 = tree_MC[i]->GetEntries();
      for (Long64_t j=0;j<nentries1;j++) {
        tree_MC[i]->GetEntry(j);
        
        //if(j<100){
        //        cout<<"For MC: mcweight: "<< mcweight << " weight: "<< weight <<endl;
        //}
        
        categ[0]=sqrt(cand_refit_tau_massE)/tripletMass >=0.000 && sqrt(cand_refit_tau_massE)/tripletMass <0.007;
        categ[1]=sqrt(cand_refit_tau_massE)/tripletMass >=0.007 && sqrt(cand_refit_tau_massE)/tripletMass <0.012;
        categ[2]=sqrt(cand_refit_tau_massE)/tripletMass >=0.012 && sqrt(cand_refit_tau_massE)/tripletMass <9999.;
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max && categ[i] && year==18 && tau_sv_ls>2.0 ){
                //MC
                  if(bdt_cv>Loose_BDT_Cut[i]){
                    tau_T3Mu[i]->Fill(tripletMass,mcweight*MC_NORM18);
                  }
                  tau_BDT_Output_MC[i]->Fill(bdt_cv,mcweight*MC_NORM18);
        }
        
      }
      
      //Bkg
      tree_bkg[i]->SetBranchAddress("cand_refit_tau_mass",&tripletMass);
      tree_bkg[i]->SetBranchAddress("cand_refit_tau_massE",&cand_refit_tau_massE);
      tree_bkg[i]->SetBranchAddress("bdt",&bdt_cv);
      //tree_bkg[i]->SetBranchAddress("mcweight",&mcweight);
      //tree_bkg[i]->SetBranchAddress("weight",&weight);
      tree_bkg[i]->SetBranchAddress("year",&year);
      tree_bkg[i]->SetBranchAddress("tau_sv_ls",&tau_sv_ls);
      
      Long64_t nentries2 = tree_bkg[i]->GetEntries();
      for (Long64_t j=0;j<nentries2;j++) {
        tree_bkg[i]->GetEntry(j);
        
        //if(j<100){
        //        cout<<"For Data: mcweight: "<< mcweight << " weight: "<< weight <<endl;
        //}
        
        categ[0]=sqrt(cand_refit_tau_massE)/tripletMass >=0.000 && sqrt(cand_refit_tau_massE)/tripletMass <0.007;
        categ[1]=sqrt(cand_refit_tau_massE)/tripletMass >=0.007 && sqrt(cand_refit_tau_massE)/tripletMass <0.012;
        categ[2]=sqrt(cand_refit_tau_massE)/tripletMass >=0.012 && sqrt(cand_refit_tau_massE)/tripletMass <9999.;
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max && categ[i] && year==18 && tau_sv_ls>2.0 ){
                //Bkg
                if( (tripletMass<=signal_peak_region_min[i] || tripletMass>=signal_peak_region_max[i]) ){//blinded
                  if(bdt_cv>Loose_BDT_Cut[i]){
                    tau_T3Mu_Dat[i]->Fill(tripletMass);
                  }
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      //tau_T3Mu_vs_BDT_Output_Data[i] = (TH2D*)file_tau[i]->Get(cat_name[i]+"BDT_2Dscan_TripletMassData");
      
      
    }
    
    
    
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
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, Cat "+cat_label[i],BDT_Score_Min,1.0);
      BDTOutput_x[i]->setRange("BDT_Range",BDT_Score_Min,1.0);
      cbmean[i] = new RooRealVar("cbmean"+hname, "cbmean" , -0.5, -9.9,0.0) ;
      cbsigma[i] = new RooRealVar("cbsigma"+hname, "cbsigma" , 0.2, 0.000001, 10.0) ;
      n[i] = new RooRealVar("n"+hname, "n", 15.0, 0.5, 20);
      alpha[i] = new RooRealVar("alpha"+hname,"alpha value CB",5.0,0.01,10);
      
      cball[i] = new RooCBShape("cball"+hname, "crystal ball", *BDTOutput_x[i], *cbmean[i], *cbsigma[i], *alpha[i], *n[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*cball[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    */
    
    
    /*
    //For fitting BDT Output: Using negative power function
    //General
    RooRealVar * BDTOutput_x[3];
    RooRealVar * cbmean[3];
    RooRealVar * cbsigma[3];
    
    RooRealVar * n_exp[3];
    RooRealVar * x00[3];
    RooRealVar * coeff_1[3];
    RooRealVar * coeff_2[3];
    RooRealVar * coeff_3[3];
    
    RooGenericPdf * cball[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    //RooGenericPdf power_func("power_func", "power_func", "pow((x-x00),n_exp)",RooArgSet(*x00[i], *n_exp[i]));
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, Cat "+cat_label[i],BDT_Score_Min,1.0);
      BDTOutput_x[i]->setRange("BDT_Range",BDT_Score_Min,1.0);
      cbmean[i] = new RooRealVar("cbmean"+hname, "cbmean" , -0.5, -9.9,0.0) ;
      cbsigma[i] = new RooRealVar("cbsigma"+hname, "cbsigma" , 0.2, 0.000001, 10.0) ;
      
      n_exp[i] = new RooRealVar("n_exp"+hname, "n_exp", 15.0, 1.0, 20.0);
      x00[i] = new RooRealVar("x00"+hname,"x00",0.5,-1.0,1.0);
      coeff_1[i] = new RooRealVar("coeff_1"+hname,"coeff_1",-0.01,-0.1,0.0);
      coeff_2[i] = new RooRealVar("coeff_2"+hname,"coeff_2",0.01,0.001,0.1);
      coeff_3[i] = new RooRealVar("coeff_3"+hname,"coeff_3",0.01,0.001,0.1);
      
      //cball[i] = new RooGenericPdf("cball"+hname, "cball", "1/TMath::Power(@0+@2,@1)", RooArgSet(*BDTOutput_x[i], *n_exp[i], *coeff_1[i]));
      
      cball[i] = new RooGenericPdf("cball"+hname, "cball", "1/TMath::Power(abs(@2*@0*@0+@3*@0+@4),@1)", RooArgSet(*BDTOutput_x[i], *n_exp[i], *coeff_1[i], *coeff_2[i], *coeff_3[i]));
      
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,50000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*cball[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    */
    
    
    
    //For fitting BDT Output: Using Argus function
    RooRealVar * BDTOutput_x[3];
    RooArgusBG * cball[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    RooRealVar * m0[3];
    RooRealVar * k0[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, Cat "+cat_label[i],BDT_Score_Min,1.0);
      BDTOutput_x[i]->setRange("BDT_Range",BDT_Score_Min,1.0);
      
      m0[i] = new RooRealVar("m0"+hname,"m0",1.0, 1.0, 1.0);
      k0[i] = new RooRealVar("k0"+hname,"k0",-30, -50, 0.0);
      
      cball[i] = new RooArgusBG("cball"+hname, "cball", *BDTOutput_x[i], *m0[i], *k0[i]);
      
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 10.0,0.5,200);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*cball[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    
    
    
    
    
    
    /*
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
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, Cat "+cat_label[i],BDT_Score_Min,1.0);
      BDTOutput_x[i]->setRange("BDT_Range",BDT_Score_Min,1.0);
      bgausmean[i] = new RooRealVar("bgausmean"+hname, "bgausmean" , -0.99,-10,0.0) ;
      bgaussigma_a[i] = new RooRealVar("bgaussigma_a"+hname, "bgaussigma_a" , 0.2, 0.000001, 10.0) ;
      bgaussigma_b[i] = new RooRealVar("bgaussigma_b"+hname, "bgaussigma_b" , 0.2, 0.000001, 10.0);
      
      bgaus_dist[i] = new RooBifurGauss("bgaus_dist"+hname, "bgaus dist", *BDTOutput_x[i], *bgausmean[i], *bgaussigma_a[i], *bgaussigma_b[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,100000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*bgaus_dist[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    */
    
    
    /*
    //For fitting BDT Output in Data: Using Exponential
    //General
    RooRealVar * BDTOutput_x[3];
    RooRealVar * exp_coeff[3];
    
    
    RooExponential * exp_dist[3];
    RooDataHist * bdt_data[3];
    RooRealVar * BDTNorm[3];
    RooAddPdf * BDT_distribution[3];
    RooFitResult * fitresult_bdt[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i] = new RooRealVar("BDTOutput_x"+hname,"BDT Output, Cat "+cat_label[i],BDT_Score_Min,1.0);
      BDTOutput_x[i]->setRange("BDT_Range",BDT_Score_Min,1.0);
      exp_coeff[i] = new RooRealVar("exp_coeff"+hname, "exp_coeff" , -5, -100, -0.1) ;
      
      exp_dist[i] = new RooExponential("exp_dist"+hname, "exp dist", *BDTOutput_x[i], *exp_coeff[i]);
      bdt_data[i] = new RooDataHist("bdt_data"+hname, "bdt_data", *BDTOutput_x[i], Import(*tau_BDT_Output_Data[i]));
      BDTNorm[i] = new RooRealVar("BDTNorm"+hname, "BDTNorm", 500.0,100.0,100000);
      BDT_distribution[i] = new RooAddPdf("BDT_distribution"+hname, "BDT_distribution", RooArgList(*exp_dist[i]), RooArgList(*BDTNorm[i]));
      fitresult_bdt[i] = BDT_distribution[i]->fitTo(*bdt_data[i], Range("BDT_Range"), Save());
      
    }
    */
    
    
    
    
    //BDT Data Fit Plots
    /*
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1800, 600);
    canvas2->Divide(3, 1);
    
    RooPlot * xFrame_bdt[3];
    float BDT_Data_Max[3];
    BDT_Data_Max[0] = 25.0;
    BDT_Data_Max[1] = 25.0;
    BDT_Data_Max[2] = 25.0;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas2->cd( i+1 );
      xFrame_bdt[i] = BDTOutput_x[i]->frame();
      bdt_data[i]->plotOn(xFrame_bdt[i]);
      BDT_distribution[i]->plotOn(xFrame_bdt[i],LineColor(4),LineWidth(2), Normalization(bdt_data[i]->sumEntries("1", "BDT_Range"), RooAbsReal::NumEvent),ProjectionRange("BDT_Range"));
      xFrame_bdt[i]->SetTitle("BDT Output, Cat "+cat_label[i]);
      xFrame_bdt[i]->SetXTitle("BDT Score");
      xFrame_bdt[i]->SetYTitle("Events");
      xFrame_bdt[i]->GetXaxis()->SetRangeUser(0.0,1.0);
      xFrame_bdt[i]->GetYaxis()->SetRangeUser(0.0001,BDT_Data_Max[i]);
      //canvas2->GetPad(i+1)->SetLogy();
      xFrame_bdt[i]->Draw();
    }
    
    */
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    //For fitting BDT Output in MC: Using Exponential
    //General
    RooRealVar * exp_coeffMC[3];
    
    
    RooExponential * exp_distMC[3];
    RooDataHist * bdt_MC[3];
    RooRealVar * BDTNormMC[3];
    RooAddPdf * BDT_distributionMC[3];
    RooFitResult * fitresult_bdtMC[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      exp_coeffMC[i] = new RooRealVar("exp_coeffMC"+hname, "exp_coeffMC" , 5,1.0,100) ;
      
      exp_distMC[i] = new RooExponential("exp_distMC"+hname, "exp dist MC", *BDTOutput_x[i], *exp_coeffMC[i]);
      bdt_MC[i] = new RooDataHist("bdt_MC"+hname, "bdt_MC", *BDTOutput_x[i], Import(*tau_BDT_Output_MC[i]));
      BDTNormMC[i] = new RooRealVar("BDTNormMC"+hname, "BDTNormMC", 5.0,0.0,100);
      BDT_distributionMC[i] = new RooAddPdf("BDT_distributionMC"+hname, "BDT_distributionMC", RooArgList(*exp_distMC[i]), RooArgList(*BDTNormMC[i]));
      fitresult_bdtMC[i] = BDT_distributionMC[i]->fitTo(*bdt_MC[i], Range("BDT_Range"), Save());
      
    }
    */
    
    /*
    //For fitting BDT Output in MC: Using power function
    //General
    RooRealVar * exp_coeffMC[3];
    
    RooRealVar * n_exp_MC[3];
    RooRealVar * coeff_1_MC[3];
    RooRealVar * coeff_2_MC[3];
    RooRealVar * coeff_3_MC[3];
    
    
    RooGenericPdf * exp_distMC[3];
    RooDataHist * bdt_MC[3];
    RooRealVar * BDTNormMC[3];
    RooAddPdf * BDT_distributionMC[3];
    RooFitResult * fitresult_bdtMC[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i]->setRange("BDT_Range_MC",BDT_Score_Min,1.0);
      
      exp_coeffMC[i] = new RooRealVar("exp_coeffMC"+hname, "exp_coeffMC" , 5,1.0,100) ;
      n_exp_MC[i] = new RooRealVar("n_exp_MC"+hname, "n_exp_MC", 15.0, 1.0, 20.0);
      coeff_1_MC[i] = new RooRealVar("coeff_1_MC"+hname,"coeff_1_MC",0.001,0.0001,0.01);
      coeff_2_MC[i] = new RooRealVar("coeff_2_MC"+hname,"coeff_2_MC",0.01,0.001,0.1);
      coeff_3_MC[i] = new RooRealVar("coeff_3_MC"+hname,"coeff_3_MC",0.01,0.001,0.1);
      
      exp_distMC[i] = new RooGenericPdf("exp_distMC"+hname, "exp dist MC", "1/TMath::Power(@4*(1.0-@0)+@2-@3*@0*@0,@1)", RooArgSet(*BDTOutput_x[i], *n_exp_MC[i], *coeff_3_MC[i], *coeff_1_MC[i], *coeff_2_MC[i]));
      
      bdt_MC[i] = new RooDataHist("bdt_MC"+hname, "bdt_MC", *BDTOutput_x[i], Import(*tau_BDT_Output_MC[i]));
      BDTNormMC[i] = new RooRealVar("BDTNormMC"+hname, "BDTNormMC", 5.0,0.0,100);
      BDT_distributionMC[i] = new RooAddPdf("BDT_distributionMC"+hname, "BDT_distributionMC", RooArgList(*exp_distMC[i]), RooArgList(*BDTNormMC[i]));
      fitresult_bdtMC[i] = BDT_distributionMC[i]->fitTo(*bdt_MC[i], Range("BDT_Range_MC"), Save());
      
    }
    */
    
    
    //For fitting BDT Output in MC: Using power function
    //General
    
    RooRealVar * m0_MC[3];
    RooRealVar * k0_MC[3];
    RooRealVar * p0_MC[3];
    
    
    RooArgusBG * exp_distMC[3];
    RooDataHist * bdt_MC[3];
    RooRealVar * BDTNormMC[3];
    RooAddPdf * BDT_distributionMC[3];
    RooFitResult * fitresult_bdtMC[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      BDTOutput_x[i]->setRange("BDT_Range_MC",BDT_Score_Min,1.0);
      
      m0_MC[i] = new RooRealVar("m0_MC"+hname,"m0_MC",0.9999, 0.999, 1.0);
      k0_MC[i] = new RooRealVar("k0_MC"+hname,"k0_MC",-200.0, -5000.0, -100.0);
      p0_MC[i] = new RooRealVar("p0_MC"+hname,"p0_MC",0.3,0.0,0.99);
      
      exp_distMC[i] = new RooArgusBG("exp_distMC"+hname, "exp dist MC", *BDTOutput_x[i], *m0_MC[i], *k0_MC[i], *p0_MC[i]);
      
      bdt_MC[i] = new RooDataHist("bdt_MC"+hname, "bdt_MC", *BDTOutput_x[i], Import(*tau_BDT_Output_MC[i]));
      BDTNormMC[i] = new RooRealVar("BDTNormMC"+hname, "BDTNormMC", 5.0,0.0,10.0);
      BDT_distributionMC[i] = new RooAddPdf("BDT_distributionMC"+hname, "BDT_distributionMC", RooArgList(*exp_distMC[i]), RooArgList(*BDTNormMC[i]));
      fitresult_bdtMC[i] = BDT_distributionMC[i]->fitTo(*bdt_MC[i], Range("BDT_Range_MC"), Save());
      
    }
    
    
    
    
    
    //MC BDT Fit Plots
    
    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 1800, 600);
    canvas3->Divide(3, 1);
    
    RooPlot * xFrame_bdtMC[3];
    float BDT_MC_Max[3];
    BDT_MC_Max[0] = 0.3;
    BDT_MC_Max[1] = 0.3;
    BDT_MC_Max[2] = 0.08;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas3->cd( i+1 );
      xFrame_bdtMC[i] = BDTOutput_x[i]->frame();
      bdt_MC[i]->plotOn(xFrame_bdtMC[i]);
      //BDT_distributionMC[i]->plotOn(xFrame_bdtMC[i],LineColor(4),LineWidth(2), Normalization(bdt_MC[i]->sumEntries("1", "BDT_Range_MC"), RooAbsReal::NumEvent),ProjectionRange("BDT_Range_MC"));
      xFrame_bdtMC[i]->SetTitle("BDT Output, Cat "+cat_label[i]);
      xFrame_bdtMC[i]->SetXTitle("BDT Score");
      xFrame_bdtMC[i]->SetYTitle("Events");
      xFrame_bdtMC[i]->GetXaxis()->SetRangeUser(0.0,1.0);
      xFrame_bdtMC[i]->GetYaxis()->SetRangeUser(0.0000001,BDT_MC_Max[i]);
      //canvas3->GetPad(i+1)->SetLogy();
      xFrame_bdtMC[i]->Draw();
    }
    
    
    

    return;
    
    

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
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      InvMass[i] = new RooRealVar("InvMass"+hname,"InvMass, Cat "+cat_label[i],signal_region_min,signal_region_max);
      InvMass[i]->setRange("R1",signal_region_min,signal_peak_region_min[i]); //background   
      InvMass[i]->setRange("R2",signal_peak_region_max[i],signal_region_max); //background
      InvMass[i]->setRange("R3",1.73,1.82); //signal range for fitting
      InvMass[i]->setRange("R4",signal_peak_region_min[i],signal_peak_region_max[i]); //signal range for yield
      
      
      //Flat fit for data
      poly[i] = new RooPolynomial("poly"+hname, "poly dist", *InvMass[i]);
      data[i] = new RooDataHist("data"+hname, "data", *InvMass[i], Import(*tau_T3Mu_Dat[i]));
      LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,15);
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
      //data[i]->plotOn(xFrame[i]);
      //pdf[i]->plotOn(xFrame[i],LineColor(4),LineWidth(2), Normalization(nData[i], RooAbsReal::NumEvent),ProjectionRange("R1,R2"));
      mc[i]->plotOn(xFrame[i]);
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), Cat "+cat_label[i]);
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
    BDTOutput_x[0]->setRange("R_Small",0.309056,1.0); BDTOutput_x[0]->setRange("R_Large",Loose_BDT_Cut[0],1.0);
    BDTOutput_x[1]->setRange("R_Small",0.329983,1.0); BDTOutput_x[1]->setRange("R_Large",Loose_BDT_Cut[1],1.0);
    BDTOutput_x[2]->setRange("R_Small",0.333186,1.0); BDTOutput_x[2]->setRange("R_Large",Loose_BDT_Cut[2],1.0);
    
    /*
    double TightBDTCutFraction[3];
    
    RooAddPdf BDT_distributionTest = *BDT_distribution[0];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      TightBDTCutFraction[i] = BDT_distribution[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Small"))->getVal() / BDT_distribution[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Large"))->getVal();
      cout << "Yield mc "+print_label[i]+" : " << BDT_distributionMC[i]->createIntegral(*BDTOutput_x[i],NormSet(*BDTOutput_x[i]),Range("R_Small"))->getVal() * (bdt_MC[i]->sumEntries("1", "R1")) << " Yield Background "+print_label[i]+" : " << (pdf_integral_restricted[i]/(1-pdf_integral_restricted[i]))*nData[i]*TightBDTCutFraction[i] << endl;
      
    }
    */
    
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
        
        //double X_min = std::min(tau_BDT_Output_Data[0]->GetXaxis()->GetXmin(), tau_BDT_Output_MC[0]->GetXaxis()->GetXmin());
        //double X_max = std::max(tau_BDT_Output_Data[0]->GetXaxis()->GetXmax(), tau_BDT_Output_MC[0]->GetXaxis()->GetXmax());
        
        /*
        double Xb_min[3];
        double Xb_max[3];
        
        Xb_min[0] = 0.965;
        Xb_min[1] = 0.965;
        Xb_min[2] = 0.9575;
        Xb_max[0] = 0.965;
        Xb_max[1] = 0.965;
        Xb_max[2] = 0.9575;
        
        Xa_min[0] = 0.995;
        Xa_min[1] = 0.998;
        Xa_min[2] = 0.994;
        Xa_max[0] = 0.995;
        Xa_max[1] = 0.998;
        Xa_max[2] = 0.994;
        
        Xa_min[0] = 0.994;
        Xa_min[1] = 0.994;
        Xa_min[2] = 0.991;
        Xa_max[0] = 0.994;
        Xa_max[1] = 0.994;
        Xa_max[2] = 0.991;
        
        Xa_min[0] = 0.99;
        Xa_min[1] = 0.99;
        Xa_min[2] = 0.99;
        Xa_max[0] = 1.00;
        Xa_max[1] = 1.00;
        Xa_max[2] = 1.00;
        */
        
        
        double Xa_min[3];
        double Xa_max[3];
        
        
        Xa_min[0] = 0.99;
        Xa_min[1] = 0.99;
        Xa_min[2] = 0.99;
        Xa_max[0] = 1.00;
        Xa_max[1] = 1.00;
        Xa_max[2] = 1.00;
        
        //Loop on both cuts in [X_min;X_max]
        Int_t dim = 0;
        //Increase N to increase (a,b) scan granularity!
        Int_t N_a = 10; //Int_t N_b = 1;//ideal: N = 16
        double step_a[3];
        double step_b[3];
        
    
    for(int iter_lumi=0; iter_lumi<lumi_size; iter_lumi++){//iterating over different values of lumis
        
        
        for(int i=0; i<3; i++){
          hname=to_string(i+1);
          
          step_a[i]=(Xa_max[i] - Xa_min[i])/N_a;
          //step_b[i]=(Xb_max[i] - Xb_min[i])/N_b;
          
          bdt_cut_vs_limit[iter_lumi][i] = new TH1D("bdt_cut_vs_limit","bdt_cut_vs_limit_"+hname,N_a,Xa_min[i],Xa_max[i]);
              
        }
        
        for(int i=0; i<N_a; i++){
                //for(int j=0; j<N_b; j++){
                        
                        for(int k=0; k<3; k++){// for tauh/taumu/taue
                                
                                a[k] = Xa_min[k] + i * step_a[k];
                                //b[k] = Xb_min[k] + j * step_b[k];
                                //if(a[k]<b[k]) continue;
                                
                                cout<<"i: "<<i<<" lumi: "<<lumi_values[iter_lumi]<<endl;
                                
                                
                                BDTOutput_x[k]->setRange("R_a",a[k],1.0); //BDTOutput_x[k]->setRange("R_b",b[k],a[k]);
                                
                                /*
                                //computing areas in range [a;X_max]
                                N_s_1[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() * (BDTNormMC[k]->getValV());
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal());
                                
                                //computing areas in range [b;a]
                                N_s_2[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() * (BDTNormMC[k]->getValV());
                                N_b_2[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal());
                                */
                                
                                N_s_1[k] = bdt_MC[k]->sumEntries("1", "R_a") * lumi_values[iter_lumi]/59.0;
                                //cout<<"N_s_1: "<<N_s_1[k]<<endl;
                                N_b_1[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_a"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * lumi_values[iter_lumi]/59.0;
                                //N_s_2[k] = BDT_distributionMC[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() * (BDTNormMC[k]->getValV()) * lumi_values[iter_lumi]/59.0;
                                //N_b_2[k] = (pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]))  *  nData[k]  *   (BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_b"))->getVal() / BDT_distribution[k]->createIntegral(*BDTOutput_x[k],NormSet(*BDTOutput_x[k]),Range("R_Large"))->getVal()) * lumi_values[iter_lumi]/59.0;
                                
                                
                                //cout<<"pdf_integral_restricted: "<<pdf_integral_restricted<<" nData[k]: "<<nData[k]<<endl;
                                //cout<<"N_s_1[k]: "<<N_s_1[k]<<" N_b_1[k]: "<<N_b_1[k]<<" N_s_2[k]: "<<N_s_2[k]<<" N_b_2[k]: "<<N_b_2[k]<<endl;
                                
                                //if(false){
                                if(N_b_1[k]>0.0){
                                        //S1[k] = N_s_1[k] / sqrt(N_s_1[k] + N_b_1[k]);
                                        //S2[k] = N_s_2[k] / sqrt(N_s_2[k] + N_b_2[k]);
                                        
                                        //S1, S2, S is now used to store limits, not significances
                                        S1[k] = sqrt( (2 * (N_s_1[k] + N_b_1[k]) * log(1 + (N_s_1[k]/N_b_1[k])) ) - 2 * N_s_1[k] );
                                        //S2[k] = sqrt( (2 * (N_s_2[k] + N_b_2[k]) * log(1 + (N_s_2[k]/N_b_2[k])) ) - 2 * N_s_2[k] );
                                        
                                        //Combined significance
                                        //S[k] = sqrt(S1[k]*S1[k] + S2[k]*S2[k]);
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        command_a[k] = "python Combine/jupyternb/"+ card_modifier_name[k] +".py --luminosity " + std::to_string(59.0) + " --s " + std::to_string(N_s_1[k]) + " --b " + std::to_string(N_b_1[k]) + " --cuttype a";
                                        //command_b[k] = "python Combine/jupyternb/"+ card_modifier_name[k] +".py --luminosity " + std::to_string(59.0) + " --s " + std::to_string(N_s_2[k]) + " --b " + std::to_string(N_b_2[k]) + " --cuttype b";
                                        system(command_a[k]);
                                        //system(command_b[k]);
                                        
                                        //command_a_and_b[k] = "combineCards.py "+combined_card_name[k]+"_a.txt "+combined_card_name[k]+"_b.txt > "+combined_card_name[k]+".txt";
                                        //system(command_a_and_b[k]);
                                        command_run[k] = "combine -M AsymptoticLimits "+combined_card_name[k]+"_a.txt --cl 0.9 -t -1  > out"+ to_string(k+1) +".txt";
                                        system(command_run[k]);
                                        
                                        
                                        
                                        std::ifstream f1("out"+ to_string(k+1) +".txt");
                                        std::string line;
                                        while (std::getline(f1, line)) {
                                          if (line.find("50.0%") != std::string::npos) {
                                            std::vector<std::string> linsp;
                                            std::istringstream iss(line);
                                            for (std::string s; iss >> s;) linsp.push_back(s);
                                            //cout << "linsp mu: " << linsp.back() << endl;
                                            S[k] = std::stod(linsp.back());
                                          }
                                        }
                                        
                                        //Not useful since we don't know the best one
                                        //command_copy[k] = "cp "+combined_card_name[k]+".txt Output/Lumi_"+combined_card_name[k]+".txt";
                                        //system(command_copy[k]);
                                        
                                        //cout<<"N_s_1[k]: "<<N_s_1[k]<<" N_b_1[k]: "<<N_b_1[k]<<" N_s_2[k]: "<<N_s_2[k]<<" N_b_2[k]: "<<N_b_2[k]<<endl;
                                        //cout<<"S: "<<S[k]<<endl;
                                        
                                        
                                        if(!isnan(S[k])){
                                        
                                                //cout<<"N_s_1[k]: "<<N_s_1[k]<<" N_b_1[k]: "<<N_b_1[k]<<" N_s_2[k]: "<<N_s_2[k]<<" N_b_2[k]: "<<N_b_2[k]<<endl;
                                                //cout<<"S: "<<S[k]<<endl;
                                                
                                                //S1_list[k].push_back(S1[k]);
                                                //S2_list[k].push_back(S2[k]);
                                                a_list[iter_lumi][k].push_back(a[k]);
                                                //b_list[iter_lumi][k].push_back(b[k]);
                                                S_list[iter_lumi][k].push_back(S[k]);
                                                
                                                
                                                N_s_1_yield_list[iter_lumi][k].push_back(N_s_1[k]);
                                                //N_s_2_yield_list[iter_lumi][k].push_back(N_s_2[k]);
                                                N_b_1_yield_list[iter_lumi][k].push_back(N_b_1[k]);
                                                //N_b_2_yield_list[iter_lumi][k].push_back(N_b_2[k]);
                                                
                                                dim++;
                                                bdt_cut_vs_limit[iter_lumi][k]->Fill(a[k]+0.00000001,S[k] );
                                        
                                        }//if!isnan
                                }//if(N_b_1[k]>0.0&&N_b_2[k]>0.0)
                        
                        }
                        
                //}
        }
    }//end iter_lumi
    //Taking absolute maximum of the combined significance
    
    
    
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
      //b_max[i] = b_list[0][i].at(S_maxIndex[i]);
      
      cout<<" a cut: "<<a_max[i]<<endl;
      cout<<"a signal yield: "<< N_s_1_yield_list[0][i].at(S_maxIndex[i]) <<" a bkg yield: "<< N_b_1_yield_list[0][i].at(S_maxIndex[i]) <<endl;
      //cout<<"b,a signal yield: "<< N_s_2_yield_list[0][i].at(S_maxIndex[i]) <<" b,a bkg yield: "<< N_b_2_yield_list[0][i].at(S_maxIndex[i]) <<endl;
      
      //a_array[i] = &a_list[i][0];
      //sig_a_array[i] = &S1_list[i][0];
      
    }
    
    
    // If this doesn't work, run line #33 first
    ofstream myfile;
    myfile.open("W_Lumi_Limit.txt", ios::out | ios::trunc);
    
    
    for(int iter_lumi=0; iter_lumi<lumi_size; iter_lumi++){
            for(int i=0; i<3; i++){
              hname=to_string(i+1);
              
              double min_limit_val = *min_element(S_list[iter_lumi][i].begin(), S_list[iter_lumi][i].end());
              int min_limit_idx = std::min_element(S_list[iter_lumi][i].begin(),S_list[iter_lumi][i].end()) - S_list[iter_lumi][i].begin();
              
              myfile << "Luminosity "<< lumi_values[iter_lumi]<<" limit "<< min_limit_val <<" category "<< "Cat_"+print_label[i] <<" a_cut "<<a_list[iter_lumi][i].at(min_limit_idx)<<" a_signal_yield "<< N_s_1_yield_list[iter_lumi][i].at(min_limit_idx) <<" a_bkg_yield "<< N_b_1_yield_list[iter_lumi][i].at(min_limit_idx) << "\n";
              
              //a_array[i] = &a_list[i][0];
              //sig_a_array[i] = &S1_list[i][0];
              
            }
    }
    myfile.close();
    
    
    
    
    /*
    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 1800, 600);
    canvas4->Divide(3, 1);
    
    TGraph * Sig_Plot[3];
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas4->cd( i+1 );
      Sig_Plot[i] = new TGraph(a_list[i].size(),a_array[i],sig_a_array[i]);
      Sig_Plot[i]->SetTitle("Significance vs. BDT Cut Value, Cat "+cat_label[i]+";BDT Cut Value;Significance");
      Sig_Plot[i]->Draw();
      
    }
    */
    
    
    
    /*
    //std::vector<double> lumi = {97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0};
    std::vector<double> lumi = {97.7, 129.0};
    std::vector<double> ZTT_taumu_lim;
    std::vector<double> ZTT_tauh_lim;
    std::vector<double> ZTT_taue_lim;
    TString command1 = "";
    TString command2 = "";
    TString command3 = "";
    for(int i=0; i<lumi.size(); i++){
      command1 = "python Combine/jupyternb/ZTT_taumu_test.py --luminosity " + std::to_string(lumi[i]);
      command2 = "python Combine/jupyternb/ZTT_tauh_test.py --luminosity " + std::to_string(lumi[i]);
      command3 = "python Combine/jupyternb/ZTT_taue_test.py --luminosity " + std::to_string(lumi[i]);
      system(command1);
      system(command2);
      system(command3);
      system("combine -M AsymptoticLimits ZTT_taumu_Combined_Mod.txt --cl 0.9 -t -1  > out1.txt");
      system("combine -M AsymptoticLimits ZTT_tauh_Combined_Mod.txt --cl 0.9 -t -1 > out2.txt");
      system("combine -M AsymptoticLimits ZTT_taue_Combined_Mod.txt --cl 0.9 -t -1 > out3.txt");
      
      
      std::ifstream f1("out1.txt");
      std::string line_mu;
      while (std::getline(f1, line_mu)) {
        if (line_mu.find("50.0%") != std::string::npos) {
          std::vector<std::string> linsp;
          std::istringstream iss(line_mu);
          for (std::string s; iss >> s;) linsp.push_back(s);
          cout << "linsp mu: " << linsp.back() << endl;
          ZTT_taumu_lim.push_back(std::stod(linsp.back()));
        }
      }
      
      std::ifstream f2("out2.txt");
      std::string line_h;
      while (std::getline(f2, line_h)) {
        if (line_h.find("50.0%") != std::string::npos) {
          std::vector<std::string> linsp;
          std::istringstream iss(line_h);
          for (std::string s; iss >> s;) linsp.push_back(s);
          cout << "linsp h: " << linsp.back() << endl;
          ZTT_tauh_lim.push_back(std::stod(linsp.back()));
        }
      }
      
      std::ifstream f3("out3.txt");
      std::string line_e;
      while (std::getline(f3, line_e)) {
        if (line_e.find("50.0%") != std::string::npos) {
          std::vector<std::string> linsp;
          std::istringstream iss(line_e);
          for (std::string s; iss >> s;) linsp.push_back(s);
          cout << "linsp e: " << linsp.back() << endl;
          ZTT_taue_lim.push_back(std::stod(linsp.back()));
        }
      }
      
    }
    */
    
    /*
    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 1800, 600);
    canvas5->Divide(3, 1);
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      canvas5->cd( i+1 );
      //bdt_cut_vs_limit[i]->SetTitle("Limit vs. BDT Cut Values, Cat "+cat_label[i]+";a;b");
      bdt_cut_vs_limit[i]->SetTitle("Limit vs. BDT Cut Values, Cat "+cat_label[i]+";a;b");
      bdt_cut_vs_limit[i]->SetStats(0);
      bdt_cut_vs_limit[i]->SetMinimum(S_max[i]);
      bdt_cut_vs_limit[i]->Draw("colz");
      
    }
    */
    
    
}