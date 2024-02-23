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
    //system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv");
    
    //system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/ZTT/CMSSW_11_3_4/src/; cmsenv");
    
    //specify the luminosities here
    //double lumi_values[] = {97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0};
    double lumi_values[] = {4500.0};
    int lumi_size = sizeof(lumi_values)/ sizeof(double);
    
    //Category names: 0=A, 1=B, 2=C
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
    
    std::string card_modifier_name[3];
    card_modifier_name[0] = "ZTT_tauh_test";
    card_modifier_name[1] = "ZTT_taumu_test";
    card_modifier_name[2] = "ZTT_taue_test";
    
    std::string combined_card_name[3];
    combined_card_name[0] = "Cat_1_Mod";
    combined_card_name[1] = "Cat_2_Mod";
    combined_card_name[2] = "Cat_3_Mod";
    
    TString hname;
    
    float signal_region_min(1.6);
    float signal_region_max(2.0);
    
    float signal_peak_region_min[3];
    float signal_peak_region_max[3];
    
    float BDT_Score_Min(0.98);
    
    int triplet_mass_bins(40);
    
    float Loose_BDT_Cut[3];
    
    Loose_BDT_Cut[0] = 0.98;
    Loose_BDT_Cut[1] = 0.98;
    Loose_BDT_Cut[2] = 0.98;
    
    
    /*
    Loose_BDT_Cut[0] = 0.995;
    Loose_BDT_Cut[1] = 0.998;
    Loose_BDT_Cut[2] = 0.994;
    
    
    Loose_BDT_Cut[0] = 0.993;
    Loose_BDT_Cut[1] = 0.994;
    Loose_BDT_Cut[2] = 0.991;
    
    */
    
    
    
    float Gaussian_Sigma_From_Loose_BDT_Cut[3];
    Gaussian_Sigma_From_Loose_BDT_Cut[0] = 0.0115727;
    Gaussian_Sigma_From_Loose_BDT_Cut[1] = 0.0174467;
    Gaussian_Sigma_From_Loose_BDT_Cut[2] = 0.025071;
    
    float Sigma_Multiplier[3];
    Sigma_Multiplier[0] = 3.5;
    Sigma_Multiplier[1] = 2.3;
    Sigma_Multiplier[2] = 1.6;
    
    for(int i=0; i<3; i++){
      hname=to_string(i+1);
      
      signal_peak_region_min[i]=1.77686-Sigma_Multiplier[i]*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      signal_peak_region_max[i]=1.77686+Sigma_Multiplier[i]*Gaussian_Sigma_From_Loose_BDT_Cut[i];
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
      
      int val_1 = (int)(((signal_peak_region_min[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      int val_2 = (int)(((signal_peak_region_max[i]-signal_region_min)/((signal_region_max-signal_region_min)/triplet_mass_bins)) + .5);
      
      signal_peak_region_min[i]=signal_region_min + val_1 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      signal_peak_region_max[i]=signal_region_min + val_2 * ((signal_region_max-signal_region_min)/triplet_mass_bins);
      
      cout<<"min: "<< signal_peak_region_min[i] <<" max: "<< signal_peak_region_max[i] <<endl;
      
      // Luca used: min: 1.74 max: 1.82 for all 3 categories. ALso, luca uses the full mass window when computing signal yield
      //signal_peak_region_min[i]=1.6;
      //signal_peak_region_max[i]=2.0;
      
    }
    
    double MC_NORM17 = 30541./(500e+3)*(8580+11370)*0.1138/0.1063*1E-7;
    double MC_NORM18 = 59828./(492e+3)*(8580+11370)*0.1138/0.1063*1E-7;
    
    double phi_veto_min = 1.020 - 0.020;
    double phi_veto_max = 1.020 + 0.020;
    
    double omega_veto_min = 0.782 - 0.020;
    double omega_veto_max = 0.782 + 0.020;
    
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
    
    Double_t cand_refit_mass12;
    Double_t cand_refit_mass13;
    Double_t cand_refit_mass23;
    Int_t cand_charge12;
    Int_t cand_charge13;
    Int_t cand_charge23;
    
    bool categ[3];//Categ A B or C
    
    bool mass_veto(false);
    
    
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
      
      tree_MC[i]->SetBranchAddress("cand_refit_mass12",&cand_refit_mass12);
      tree_MC[i]->SetBranchAddress("cand_refit_mass13",&cand_refit_mass13);
      tree_MC[i]->SetBranchAddress("cand_refit_mass23",&cand_refit_mass23);
      tree_MC[i]->SetBranchAddress("cand_charge12",&cand_charge12);
      tree_MC[i]->SetBranchAddress("cand_charge13",&cand_charge13);
      tree_MC[i]->SetBranchAddress("cand_charge23",&cand_charge23);
      
      
      
      Long64_t nentries1 = tree_MC[i]->GetEntries();
      for (Long64_t j=0;j<nentries1;j++) {
        tree_MC[i]->GetEntry(j);
        
        categ[0]=sqrt(cand_refit_tau_massE)/tripletMass >=0.000 && sqrt(cand_refit_tau_massE)/tripletMass <0.007;
        categ[1]=sqrt(cand_refit_tau_massE)/tripletMass >=0.007 && sqrt(cand_refit_tau_massE)/tripletMass <0.012;
        categ[2]=sqrt(cand_refit_tau_massE)/tripletMass >=0.012 && sqrt(cand_refit_tau_massE)/tripletMass <9999.;
        
        mass_veto = ( (cand_charge12==0)? (cand_refit_mass12>=phi_veto_max||cand_refit_mass12<=phi_veto_min)&&(cand_refit_mass12>=omega_veto_max||cand_refit_mass12<=omega_veto_min) : true ) && ( (cand_charge13==0)? (cand_refit_mass13>=phi_veto_max||cand_refit_mass13<=phi_veto_min)&&(cand_refit_mass13>=omega_veto_max||cand_refit_mass13<=omega_veto_min) : true ) && ( (cand_charge23==0)? (cand_refit_mass23>=phi_veto_max||cand_refit_mass23<=phi_veto_min)&&(cand_refit_mass23>=omega_veto_max||cand_refit_mass23<=omega_veto_min) : true );
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max && categ[i] && year==18 && tau_sv_ls>2.0 && mass_veto ){
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
      tree_bkg[i]->SetBranchAddress("year",&year);
      tree_bkg[i]->SetBranchAddress("tau_sv_ls",&tau_sv_ls);
      
      tree_bkg[i]->SetBranchAddress("cand_refit_mass12",&cand_refit_mass12);
      tree_bkg[i]->SetBranchAddress("cand_refit_mass13",&cand_refit_mass13);
      tree_bkg[i]->SetBranchAddress("cand_refit_mass23",&cand_refit_mass23);
      tree_bkg[i]->SetBranchAddress("cand_charge12",&cand_charge12);
      tree_bkg[i]->SetBranchAddress("cand_charge13",&cand_charge13);
      tree_bkg[i]->SetBranchAddress("cand_charge23",&cand_charge23);
      
      Long64_t nentries2 = tree_bkg[i]->GetEntries();
      for (Long64_t j=0;j<nentries2;j++) {
        tree_bkg[i]->GetEntry(j);
        
        categ[0]=sqrt(cand_refit_tau_massE)/tripletMass >=0.000 && sqrt(cand_refit_tau_massE)/tripletMass <0.007;
        categ[1]=sqrt(cand_refit_tau_massE)/tripletMass >=0.007 && sqrt(cand_refit_tau_massE)/tripletMass <0.012;
        categ[2]=sqrt(cand_refit_tau_massE)/tripletMass >=0.012 && sqrt(cand_refit_tau_massE)/tripletMass <9999.;
        
        mass_veto = ( (cand_charge12==0)? (cand_refit_mass12>=phi_veto_max||cand_refit_mass12<=phi_veto_min)&&(cand_refit_mass12>=omega_veto_max||cand_refit_mass12<=omega_veto_min) : true ) && ( (cand_charge13==0)? (cand_refit_mass13>=phi_veto_max||cand_refit_mass13<=phi_veto_min)&&(cand_refit_mass13>=omega_veto_max||cand_refit_mass13<=omega_veto_min) : true ) && ( (cand_charge23==0)? (cand_refit_mass23>=phi_veto_max||cand_refit_mass23<=phi_veto_min)&&(cand_refit_mass23>=omega_veto_max||cand_refit_mass23<=omega_veto_min) : true );
        
        if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max && categ[i] && year==18 && tau_sv_ls>2.0 && mass_veto ){
                //Bkg
                if( (tripletMass<=signal_peak_region_min[i] || tripletMass>=signal_peak_region_max[i]) ){//blinded
                  if(bdt_cv>Loose_BDT_Cut[i]){
                    tau_T3Mu_Dat[i]->Fill(tripletMass);
                  }
                  tau_BDT_Output_Data[i]->Fill(bdt_cv);
                }
        }
        
      }
      
      
    }
    
    
    
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
    /*
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
      LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,200.0);
      pdf[i] = new RooAddPdf("pdf"+hname, "pdf", RooArgList(*poly[i]), RooArgList(*LineNorm[i]));
      fitresult[i] = pdf[i]->fitTo(*data[i], Range("R1,R2"), Save());
      
      
      //Gaussian fit for MC
      mean[i] = new RooRealVar("mean"+hname, "mean" , 1.776,0.,5.) ;
      sigma[i] = new RooRealVar("sigma"+hname, "sigma" , 0.5,0.001,10) ;
      
      Gauss[i] = new RooGaussian("Gauss"+hname, "Gauss dist", *InvMass[i], *mean[i], *sigma[i]);
      mc[i] = new RooDataHist("mc"+hname, "mc", *InvMass[i], Import(*tau_T3Mu[i]));
      GaussNorm[i] = new RooRealVar("GaussNorm"+hname, "GaussNorm",  0.5,0.001,10.0);
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
      //mc[i]->plotOn(xFrame[i]);
      mc_pdf[i]->plotOn(xFrame[i],LineColor(1),LineWidth(2), Normalization(nSignal[i], RooAbsReal::NumEvent),ProjectionRange("R3"));
      xFrame[i]->SetTitle("3#mu inv. mass (GeV), Cat "+cat_label[i]);
      xFrame[i]->SetXTitle("3#mu inv. mass (GeV)");
      xFrame[i]->SetYTitle("Events");
      xFrame[i]->Draw();
    }
    */
    
    
    
    
    
    
    
    
    //From BDT Output:
    BDTOutput_x[0]->setRange("R_Small",0.309056,1.0); BDTOutput_x[0]->setRange("R_Large",Loose_BDT_Cut[0],1.0);
    BDTOutput_x[1]->setRange("R_Small",0.329983,1.0); BDTOutput_x[1]->setRange("R_Large",Loose_BDT_Cut[1],1.0);
    BDTOutput_x[2]->setRange("R_Small",0.333186,1.0); BDTOutput_x[2]->setRange("R_Large",Loose_BDT_Cut[2],1.0);
    
    
    
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
        
        
        
        Xa_min[0] = 0.985;
        Xa_min[1] = 0.985;
        Xa_min[2] = 0.985;
        Xa_max[0] = 0.998;
        Xa_max[1] = 0.998;
        Xa_max[2] = 0.998;
        
        
        //Xa_min[0] = 0.995;
        //Xa_min[1] = 0.998;
        //Xa_min[2] = 0.994;    
        
        
        
        //For loop on both cuts in [X_min;X_max]
        Int_t dim = 0;
        //Increase N to increase (a,b) scan granularity!
        Int_t N_a = 13; //Int_t N_b = 1;//ideal: N = 10
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
                                
                                BDTOutput_x[k]->setRange("R_a",a[k],1.0);
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
    myfile.open("W_Lumi_Limit.txt", ios::out | ios::trunc);
    
    
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
    //
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