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

void makeYield_fromBDTFit_Combine_PeakWidth () 
{
    //system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/HF_and_W/CMSSW_10_2_13/src; cmsenv; cd -");
    system("cd /afs/cern.ch/work/m/mmadhu/Analysis/combinestats/t3mcombine/ZTT/CMSSW_11_3_4/src/; cmsenv; cd -");
    
    //specify the luminosities here
    //double lumi_values[] = {59.0, 97.7, 129.0, 377.0, 700.0, 1500.0, 2250.0, 3000.0, 3750.0, 4500.0};
    double lumi_values[] = {59.0};
    int lumi_size = sizeof(lumi_values)/ sizeof(double);
    
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
    
    std::string card_modifier_name[4];
    card_modifier_name[0] = "ZTT_tauha_test";
    card_modifier_name[1] = "ZTT_tauhb_test";
    card_modifier_name[2] = "ZTT_taumu_test";
    card_modifier_name[3] = "ZTT_taue_test";
    
    std::string combined_card_name[4];
    combined_card_name[0] = "Cat_1a_Mod";
    combined_card_name[1] = "Cat_1b_Mod";
    combined_card_name[2] = "Cat_2_Mod";
    combined_card_name[3] = "Cat_3_Mod";
    
    TString hname;
    
    float signal_region_min(1.6);
    float signal_region_max(2.0);
    
    //float signal_peak_region_min(1.73);
    //float signal_peak_region_max(1.82);
    
    float signal_peak_region_min[4];
    float signal_peak_region_max[4];
    
    float BDT_Score_Min(-0.9);
    
    int triplet_mass_bins(40);
    
    float Loose_BDT_Cut(-0.4);
    
    
    //BDT Cuts a
    float BDT_Cut_a[4];
    
    //Cuts from limits
    
    BDT_Cut_a[0] = 0.36;
    BDT_Cut_a[1] = 0.14;
    BDT_Cut_a[2] = 0.35;
    BDT_Cut_a[3] = 0.18;
    
    
    //Cuts from significance
    /*
    BDT_Cut_a[0] = 0.3825;
    BDT_Cut_a[1] = 0.26;
    BDT_Cut_a[2] = 0.44;
    BDT_Cut_a[3] = 0.27;
    */
    
    
    
    float Gaussian_Sigma_From_Loose_BDT_Cut[4];
    Gaussian_Sigma_From_Loose_BDT_Cut[0] = 0.0170238;
    Gaussian_Sigma_From_Loose_BDT_Cut[1] = 0.0167713;
    Gaussian_Sigma_From_Loose_BDT_Cut[2] = 0.0172921;
    Gaussian_Sigma_From_Loose_BDT_Cut[3] = 0.0172922;
    
    double Xa_min[4];
    double Xa_max[4];
    
    Xa_min[0] = 1.5;
    Xa_min[1] = 1.5;
    Xa_min[2] = 1.5;
    Xa_min[3] = 1.5;
    Xa_max[0] = 3.5;
    Xa_max[1] = 3.5;
    Xa_max[2] = 3.5;
    Xa_max[3] = 3.5;
    
    Int_t N_a = 20;
    
    double sigma_scale[4];
    double step[4];
    
    std::vector<double> S_list[4], a_list[4];
    
    for(int m=0; m<N_a; m++){
            
            for(int i=0; i<4; i++){
              hname=to_string(i+1);
              
              step[i] = (Xa_max[i] - Xa_min[i])/N_a;
              
              sigma_scale[i] = Xa_min[i] + m * step[i];
              
              cout<< "sigma_scale "<< i << " is: " << sigma_scale[i] <<endl;
              
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
            
            TH1D  * tau_T3Mu[4];
            TH1D  * tau_T3Mu_Dat[4];
            TH1D  * tau_BDT_Output_Data[4];
            TH1D  * tau_BDT_Output_MC[4];
            TH2D  * tau_T3Mu_vs_BDT_Output_Data[4];
            TH1D  * tau_T3Mu_vs_BDT_Output_Data_Projection[4];
            
            TH1D  * bdt_cut_vs_sig[lumi_size][4];
            TH1D  * bdt_cut_vs_limit[lumi_size][4];
            
            
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
              tau_BDT_Output_MC[i] = new TH1D("tau_BDT_Output_MC","tau_BDT_Output_MC_"+hname,100,BDT_Score_Min,0.999);
              tau_BDT_Output_Data[i] = new TH1D("tau_BDT_Output_Data","tau_BDT_Output_Data_"+hname,100,BDT_Score_Min,0.999);
              
              tree[i]->SetBranchAddress("tripletMass",&tripletMass);
              tree[i]->SetBranchAddress("bdt_cv",&bdt_cv);
              tree[i]->SetBranchAddress("weight",&weight);
              tree[i]->SetBranchAddress("isMC",&isMC);
              
              Long64_t nentries = tree[i]->GetEntries();
              for (Long64_t j=0;j<nentries;j++) {
                tree[i]->GetEntry(j);
                
                //if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max){
                if(tripletMass>=signal_region_min && tripletMass<=signal_region_max && bdt_cv>=BDT_Cut_a[i]){
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
              LineNorm[i] = new RooRealVar("LineNorm"+hname, "LineNorm", 2.0,0.001,20000.0);
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
            
            
            double pdf_integral_restricted[4];
            
            double N_s_1[4], N_b_1[4];
            double S[4];
            
            TString command_a[4];
            TString command_Bayesian[4];
            TString command_UL_Calc[4];
            
            for(int k=0; k<4; k++){// for tauh A B/taumu/taue or A/B/C
                    
                    cout<< "sigma_scale "<< k << " is: " << sigma_scale[k] <<endl;
                    
                    pdf_integral_restricted[k] = pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getVal();
                    
                    //N_s_1[k] = GaussNorm[k]->getValV()* mc_pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getVal();
                    //N_b_1[k] = LineNorm[k]->getValV()* pdf[k]->createIntegral(*InvMass[k],NormSet(*InvMass[k]),Range("R4"))->getVal();
                    
                        N_s_1[k] = 0.0;
                        N_b_1[k] = 0.0;
                        
                        //Getting signal and bkg from the ttree rather than the histograms (histograms are still not fully correct, even though they aren't coming from a fit)
                        Long64_t nentries = tree[k]->GetEntries();
                        for (Long64_t j=0;j<nentries;j++) {
                                tree[k]->GetEntry(j);
                                
                                if(tripletMass>=signal_region_min&&tripletMass<=signal_region_max){
                                        if(isMC>0){
                                          if(bdt_cv>BDT_Cut_a[k] && (tripletMass>=signal_peak_region_min[k] && tripletMass<=signal_peak_region_max[k]) ){
                                            N_s_1[k]+=weight;
                                          }
                                        }
                                        if(isMC==0){
                                          if(bdt_cv>BDT_Cut_a[k] && (tripletMass<=signal_peak_region_min[k] || tripletMass>=signal_peak_region_max[k]) ){//blinded
                                            N_b_1[k]+=1.0;
                                          }
                                        }
                                }
                                
                        }
                        
                        N_b_1[k] = N_b_1[k]*(pdf_integral_restricted[k]/(1-pdf_integral_restricted[k]));
                    
                    cout<<"signal: "<<N_s_1[k]<<" bkg: "<< N_b_1[k] <<endl;
                    
                    if( (std::round(N_b_1[k] / 0.00001) * 0.00001) > 0.0){
                    
                            //Create datacard from signal and bkg values
                            command_a[k] = "python card_modifiers/"+ card_modifier_name[k] +".py --luminosity " + std::to_string(59.0) + " --s " + std::to_string(N_s_1[k]) + " --b " + std::to_string(N_b_1[k]) + " --cuttype a";
                            system(command_a[k]);
                            
                            
                            bool Whether_Bayesian(true);
                            
                            //BayesianSimple
                            if(Whether_Bayesian){
                                    
                                    command_Bayesian[k] = "combine -M BayesianSimple "+combined_card_name[k]+"_a.txt --cl 0.9 -t 100  > Test_Bayesian_"+ to_string(k+1) +".txt";
                                    system(command_Bayesian[k]);
                            
                                    std::ifstream bayes("Test_Bayesian_"+ to_string(k+1) +".txt");
                                    std::string line2;
                                    while (std::getline(bayes, line2)) {
                                            if (line2.find("median expected limit") != std::string::npos) {
                                            //if (line2.find("90% CL") != std::string::npos) {
                                                    std::vector<std::string> linsp;
                                                    std::istringstream iss(line2);
                                                    for (std::string s; iss >> s;) linsp.push_back(s);
                                                    //cout<<"Limit Bayesian: "<<linsp[3]<<endl;
                                                    cout<<"Limit Bayesian: "<<linsp[5]<<endl;
                                                    S[k] = std::stod(linsp[5]);
                                            }
                                    }
                            }
                            
                            command_UL_Calc[k] = "python3 ../CLs_UL_Calculator.py " + std::to_string(N_s_1[k]) + " " + std::to_string(N_b_1[k]) + " > out_UL_Calc_"+ to_string(k+1) +".txt";
                            system(command_UL_Calc[k]);
                            std::ifstream UL_Cal("out_UL_Calc_"+ to_string(k+1) +".txt");
                            std::string line1;
                            while (std::getline(UL_Cal, line1)) {
                                  if (line1.find("upper limit") != std::string::npos) {
                                    std::vector<std::string> linsp;
                                    std::istringstream iss(line1);
                                    for (std::string s; iss >> s;) linsp.push_back(s);
                                    cout<<"Limit UL_Calc: "<<linsp.back()<<endl;
                                    //S[k] = std::stod(linsp.back());
                                  }
                                  if (line1.find("best value of CLs") != std::string::npos) {
                                    std::vector<std::string> linsp;
                                    std::istringstream iss(line1);
                                    for (std::string s; iss >> s;) linsp.push_back(s);
                                    //cout<<"Value of CLs: "<<linsp.back()<<endl;
                                  }
                            }
                            
                            
                            if(!isnan(S[k])){
                                                
                                    a_list[k].push_back(sigma_scale[k]);
                                    S_list[k].push_back(S[k]);
                                                
                            }//isnan check
                    
                    }//bkg non-zero check
                    
                    
                    
            }
            
            
            
            
    }// end m
    
    
    //Taking absolute minimum of the limits
    double S_min[4];
    int S_minIndex[4];
    float a_min[4];
    
    cout<<" ---Results--- "<<endl;
    
    for(int i=0; i<4; i++){
      hname=to_string(i+1);
      
      for(int k=0; k<S_list[i].size(); k++){
              cout<<"Limit Bayesian: " << S_list[i][k]<< " at sigma scale: "<< a_list[i][k] <<endl;
      }
      
      S_min[i] = *min_element(S_list[i].begin(), S_list[i].end());
      S_minIndex[i] = std::min_element(S_list[i].begin(),S_list[i].end()) - S_list[i].begin();
      cout<<"S_min[i]: "<<S_min[i]<<" S_minIndex[i]: "<<S_minIndex[i]<<endl;
      
      a_min[i] = a_list[i].at(S_minIndex[i]);
      
      cout<<"Best sigma_scale: "<<a_min[i]<<endl;
      
    }
    
    
    
}