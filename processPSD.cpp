#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TPaveLabel.h"
#include <TApplication.h>

using namespace std;

int main( int argc, char* argv[]){
  //TApplication * myapp = new TApplication("myapp",0,0);
  char *resdir = argv[1];
  char *tier2_name = argv[2];
  int chn = atoi(argv[3]);
  
  double max_spectrum = 10000;
  
  TFile *tier2_file = TFile::Open(tier2_name);
  TTree *tier2 = (TTree*)tier2_file->Get("tier2");
  double min_AoE;
  double min_AoE_corr;
  
  cout << "Analizing channel " << chn << " using Standard Analysis ..." <<endl;
  
  TCanvas *c = new TCanvas("c","energy");
  gStyle->SetPalette(1);
  //c->SetGrid();
  gPad->SetLogy();
  
  //------- draw the spectrum to find the maximum -------------
  
  TH1D *h = new TH1D("h","", 1000, max_spectrum/4,max_spectrum);
  
  tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>h",chn),Form("GEMDFADC_eventType[%d]==0",chn));
  c->Update();
  
  //--------- find parameters from spectrum -------
  int maxSpectrum = h->GetMaximumBin();
  if (maxSpectrum > 10){
    double mean = h->GetBinCenter(maxSpectrum);
    double calib = 2614.511/mean;
    cout << "Calibration = " << calib << endl;
    //DEP+FEP region
    double min = 1570 / calib;
    double max = 1640 / calib;
    int bin = 280;//( max - min ) * 1 ;
    
    //----------- draw the spectrum in the DEP+FEP range -----------
    
    TH1D *h_ene = new TH1D("h_ene","DEP region; Uncalibrated Energy [chn]; counts",bin,min,max);
    h_ene->GetXaxis()->CenterTitle();
    h_ene->GetYaxis()->CenterTitle();
    tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>h_ene",chn));
    c->Update();
    
    //----- fit of the energy peaks with a double-gaussian --------
    
    TF1 *f_2gaus = new TF1("f_2gaus","gaus(0)+gaus(3)+pol0(6)",min,max);
    
    f_2gaus->SetParName( 0, "Amp. - DEP" );
    f_2gaus->SetParName( 1, "Mean - DEP" );
    f_2gaus->SetParName( 2, "Sigma - DEP" );
    f_2gaus->SetParName( 3, "Amp. - FEP" );
    f_2gaus->SetParName( 4, "Mean - FEP" );
    f_2gaus->SetParName( 5, "Sigma - FEP" );
    f_2gaus->SetParName( 6, "Flat Bkg" );
    
    //--------- initialization of parameters of the fit ----------
    
    double A_d = h_ene->GetMaximum();
    double mu_d = h_ene->GetBinCenter( h_ene->GetMaximumBin());
    h_ene->GetXaxis()->SetRangeUser((max+min)/2,max);
    double A_f = h_ene->GetMaximum();
    double mu_f = h_ene->GetBinCenter( h_ene->GetMaximumBin());
    h_ene->GetXaxis()->SetRangeUser(min,max);
    
    f_2gaus->SetParameter(0, A_d);
    f_2gaus->SetParameter(1, mu_d);
    f_2gaus->SetParameter(2, mu_d * 0.001);
    f_2gaus->SetParLimits(2, 0, 100);
    f_2gaus->SetParameter(3, A_f);
    f_2gaus->SetParameter(4, mu_f);
    f_2gaus->SetParameter(5, mu_f * 0.001);
    f_2gaus->SetParLimits(5, 0, 100);
    f_2gaus->SetParameter(6, 10);
    
    gStyle->SetOptFit(111);
    
    h_ene->Fit("f_2gaus", "L", "", min, max);
    c->Update();
    
    double A_dep = f_2gaus->GetParameter(0);
    double mu_dep = f_2gaus->GetParameter(1);
    double sigma_dep = f_2gaus->GetParameter(2);
    double A_fep = f_2gaus->GetParameter(3);
    double mu_fep = f_2gaus->GetParameter(4);
    double sigma_fep = f_2gaus->GetParameter(5);
    double bkg = f_2gaus->GetParameter(6);
    
    //------------ calculation of Area of DEP, FEP -----------
    
    int count = 0;
    int sumCount = 0;
    double n_sigma = 4.5;//number of sigma to calculate area
    
    //----- DEP Area ------
    int min_bin = h_ene->GetXaxis()->FindBin(mu_dep - n_sigma*sigma_dep);
    int max_bin = h_ene->GetXaxis()->FindBin(mu_dep + n_sigma*sigma_dep);
    int min_bin_1 = h_ene->GetXaxis()->FindBin(mu_dep - 2*n_sigma*sigma_dep);
    int max_bin_1 = h_ene->GetXaxis()->FindBin(mu_dep + 2*n_sigma*sigma_dep);
    int delta = max_bin - min_bin;
    
    for (Int_t j = min_bin; j <= max_bin; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;
    }
    double Area_d_tot = sumCount;//total DEP area
    //cout << "AreaDEP_tot = " << Area_d << endl;
    
    //---background
    sumCount = 0;
    for (Int_t j =  min_bin_1; j <= min_bin; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;
    }
    for (Int_t j =  max_bin; j <= max_bin_1; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;
    }
    double b1_d = sumCount;
    double Area_d = Area_d_tot - b1_d; //subtraction of background
    //cout << "AreaDEP = " << Area_d << endl;
    
    //----- FEP Area ------
    min_bin = h_ene->GetXaxis()->FindBin(mu_fep - n_sigma*sigma_fep);
    max_bin = h_ene->GetXaxis()->FindBin(mu_fep + n_sigma*sigma_fep);
    min_bin_1 = h_ene->GetXaxis()->FindBin(mu_fep - 2*n_sigma*sigma_fep);
    max_bin_1 = h_ene->GetXaxis()->FindBin(mu_fep + 2*n_sigma*sigma_fep);
    delta = max_bin - min_bin;
    //cout << min_bin << " "<<  max_bin << endl;
    
    sumCount = 0;
    for (Int_t j = min_bin; j <= max_bin; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;}
    double Area_f_tot = sumCount;//total FEP area
    //cout << "AreaFEP_tot = " << Area_f << endl;
    
    //---background
    sumCount = 0;
    for (Int_t j =  min_bin_1; j <= min_bin; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;
    }
    for (Int_t j =  max_bin; j <= max_bin_1; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;
    }
    double b1_f = sumCount;
    double Area_f = Area_f_tot - b1_f;//subtraction of background
    //cout << "AreaFEP = " << Area_f << endl;
    
    //---- DEP interval ----
    double min_dep = mu_dep - n_sigma * sigma_dep;
    double max_dep = mu_dep + n_sigma * sigma_dep;
    double min_dep2 = mu_dep - 2*n_sigma * sigma_dep;
    double max_dep2 = mu_dep + 2*n_sigma * sigma_dep;
    
    ////////////////////////////////////////////////////
    //----------------- AoE_analysis ------------------
    ////////////////////////////////////////////////////
    
    TCanvas *c1 = new TCanvas("c1","AoE");
    //c1->SetGrid();
    c1->cd();
    gStyle->SetPalette(1);
    
    TCut cut_e = Form("GEMDEnergyGauss_energy[%d] > %g && GEMDEnergyGauss_energy[%d] < %g", chn, min_dep, chn, max_dep);
    TCut cut_e1 = Form("GEMDEnergyGauss_energy[%d] > %g && GEMDEnergyGauss_energy[%d] < %g",chn, min_dep2, chn, min_dep);
    TCut cut_e2 = Form("GEMDEnergyGauss_energy[%d] > %g && GEMDEnergyGauss_energy[%d] < %g",chn, max_dep, chn, max_dep2);
    //--------- find the value of the A/E ----------
    
    TH1D *AoE = new TH1D("AoE","A over E; A/E; counts",1000,0.02,0.06);
    tier2->Draw(Form("GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d]>>AoE",chn,chn),cut_e);
    c1->Update();
    
    double mean_aoe = AoE->GetBinCenter( AoE->GetMaximumBin());
    cout << "mean_aoe = " << mean_aoe << endl;
    min_AoE = mean_aoe * 0.94;
    double max_AoE = mean_aoe * 1.06;
    
    //--------------- draw the A/E of the DEP --------------
    
    TH1D *h_AoE = new TH1D("h_AoE",Form("chn%d - AoE - DEP-SEP-FEP-2039keV;A/E;counts",chn), 300, 0.7, 1.1);
    TH1D *h_AoE1 = new TH1D("h_AoE1",Form("chn%d - AoE - DEP;A/E;counts",chn), 300, 0.7, 1.1);
    TH1D *h_AoE2 = new TH1D("h_AoE2",Form("chn%d - AoE - DEP;A/E;counts",chn), 300, 0.7, 1.1);
    tier2->Draw(Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g>>h_AoE", chn, chn, mean_aoe),cut_e);	
    tier2->Draw(Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g>>h_AoE1", chn, chn, mean_aoe),cut_e1);
    tier2->Draw(Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g>>h_AoE2", chn, chn, mean_aoe),cut_e2);
    //c1->Update();
    h_AoE1->Add(h_AoE2);
    h_AoE->Add(h_AoE1,-1);
    h_AoE->Draw();
    
    h_AoE->Scale(1./h_AoE->GetBinContent(h_AoE->GetMaximumBin()));
    h_AoE->GetXaxis()->CenterTitle();
    h_AoE->GetYaxis()->CenterTitle();
    h_AoE->SetLineWidth(2);
    h_AoE->Draw();
    c1->Update();
    
    //---------------- fit the A/E of the DEP ------------
    
    TF1 *f_Agaus = new TF1("f_Agaus","gaus(0) + pol1(3)", 0.7, 1.1);
    
    double A_AoE = h_AoE->GetMaximum();
    double mu_AoE = h_AoE->GetBinCenter( h_AoE->GetMaximumBin());
    double sigma_AoE = mu_AoE * 0.01;
    
    f_Agaus->SetParName( 0, "Amp. - A/E" );
    f_Agaus->SetParName( 1, "Mean - A/E" );
    f_Agaus->SetParName( 2, "Sigma - A/E" );
    f_Agaus->SetParName( 3, "Flat Bkg" );
    f_Agaus->SetParName( 4, "Linear Bkg" );
    
    f_Agaus->SetParameter(0, 1.);//A_AoE);
    f_Agaus->SetParameter(1, 1.);//mu_AoE);
    f_Agaus->SetParameter(2, 0.01);//sigma_AoE);
    f_Agaus->SetParLimits(2, 0, 1.);
    
    gStyle->SetOptFit(111);
    
    h_AoE->Fit("f_Agaus","L0","", 0.7, 1.1);
    
    A_AoE = f_Agaus->GetParameter(0);
    mu_AoE = f_Agaus->GetParameter(1);
    sigma_AoE = f_Agaus->GetParameter(2);
    double mu_AoE_err = f_Agaus->GetParError(1);
    double sigma_AoE_err = f_Agaus->GetParError(2);
    
    double FWHM_AoE = 2.35*sigma_AoE/mu_AoE;
    double FWHM_AoE_err = (sigma_AoE_err/sigma_AoE+mu_AoE_err/mu_AoE)*FWHM_AoE;
    
    //------------ print the results of A/E ---------
    
    //cout<<"FWHM = "<< FWHM_AoE <<" +/- "<<FWHM_AoE_err<<" %"<< endl;
    
    TString resultAoE = Form("FWHM = ( %4.2f +/- %4.2f ) %%",FWHM_AoE*100,FWHM_AoE_err*100);
    TPaveLabel *label = new TPaveLabel(0.75,0.38,0.98,0.46,resultAoE,"NDC");
    label->Draw();
    c1->Update();
    
    ////////////////////////////////////////////////////
    //---------------- PSA analysis --------------------
    ////////////////////////////////////////////////////
    
    //------------------ AoE interval ------------------- 

    double totIntegral = h_AoE->Integral(h_AoE->GetXaxis()->FindBin(0.), h_AoE->GetXaxis()->FindBin(1.8));
    cout << "Integral = " << totIntegral << endl;
    TCut cut_aoe;
    double AoE_lowCut = 0.901;
    double cutIntegral;
    double SurvFrac_d;
    for (int i = 0; i < 100; i++){
    	cutIntegral = h_AoE->Integral(h_AoE->GetXaxis()->FindBin(AoE_lowCut), h_AoE->GetXaxis()->FindBin(1.07));
    	SurvFrac_d = (double)cutIntegral/totIntegral;
    	
    	if(SurvFrac_d < 0.9) break;
    	else AoE_lowCut += 0.001;
    }
    cout << "AoE_cut = " << AoE_lowCut << endl;      
    cout << "CutIntegral = " << cutIntegral << endl;
    cout << "SF_dep = " << SurvFrac_d << endl;
    
    TString cutAoE = Form("A/E_cut = %5.3f ", AoE_lowCut);
    
    TPaveLabel *label_ = new TPaveLabel(0.75, 0.26, 0.98, 0.34, cutAoE,"NDC");
    label_->Draw();
    c1->Update();
    
    double n_sigma_aoe = 2.2;
    
    TH1D *hene_psa = new TH1D("hene_psa",Form("chn%d - DEP region; uncalibrated energy; counts",chn), bin, min, max);
    
    TCanvas *c2 = new TCanvas("c2","PSA analysis");
    //c2->SetGrid();
    c2->cd();
    gStyle->SetPalette(1);
    
    tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>h_ene",chn));
    c2->Update();
    TCut cut_aoe1 = Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g > %g && (GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g <1.07", chn, chn, mean_aoe, AoE_lowCut, chn,  chn, mean_aoe);
    
    tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>hene_psa",chn),cut_aoe1,"sames");
    
    hene_psa->GetXaxis()->CenterTitle();
    hene_psa->GetYaxis()->CenterTitle();
    c2->Update();
    
    //------------- fit DEP region after PSD ----------------
    
    f_2gaus->SetParameter(0, A_dep*0.9);
    f_2gaus->SetParameter(1, mu_dep);
    f_2gaus->SetParameter(2, sigma_dep);
    f_2gaus->SetParLimits(2, 0, 100);
    f_2gaus->SetParameter(3, A_fep*0.15);
    f_2gaus->SetParameter(4, mu_fep);
    f_2gaus->SetParameter(5, sigma_fep);
    f_2gaus->SetParLimits(5, 0, 100);
    f_2gaus->SetParameter(6, bkg*0.5);
    
    gStyle->SetOptFit(111);
    
    hene_psa->Fit("f_2gaus", "L+", "sames", min, max);
    
    //-------------- get the results from fit --------------
    
    double bkg_psa = f_2gaus->GetParameter(6);
    
    //------------ calculation of Area of DEP,FEP -----------
    
    //----- DEP area -----
    sumCount = 0;
    min_bin = hene_psa->GetXaxis()->FindBin(mu_dep - n_sigma*sigma_dep);
    max_bin = hene_psa->GetXaxis()->FindBin(mu_dep + n_sigma*sigma_dep);
    min_bin_1 = hene_psa->GetXaxis()->FindBin(mu_dep - 2*n_sigma*sigma_dep);
    max_bin_1 = hene_psa->GetXaxis()->FindBin(mu_dep + 2*n_sigma*sigma_dep);
    delta = max_bin - min_bin;
    
    for (Int_t j = min_bin; j < max_bin + 1; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;}
    double Area_d_psa_tot = sumCount;//total DEP area
    //cout << "Area_DEP_tot_psa = " << Area_d_psa_tot << endl;
    //---background
    sumCount = 0;
    for (Int_t j =  min_bin_1; j <= min_bin; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;}
    for (Int_t j =  max_bin; j <= max_bin_1; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;}
    double b2_d = sumCount;
    double Area_d_psa = Area_d_psa_tot - b2_d;//subtraction of background
    //cout << "Area_DEP_psa = " << Area_d_psa << endl;
    
    //----- FEP area -----
    min_bin = hene_psa->GetXaxis()->FindBin(mu_fep - n_sigma*sigma_fep);
    max_bin = hene_psa->GetXaxis()->FindBin(mu_fep + n_sigma*sigma_fep);
    min_bin_1 = hene_psa->GetXaxis()->FindBin(mu_fep - 2*n_sigma*sigma_fep);
    max_bin_1 = hene_psa->GetXaxis()->FindBin(mu_fep + 2*n_sigma*sigma_fep);
    delta = max_bin - min_bin;
    
    sumCount = 0;
    for (Int_t j = min_bin; j <= max_bin; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;
    }
    double Area_f_psa_tot = sumCount;//total FEP area
    //cout << "Area_FEP_tot_psa = " << Area_f_psa_tot << endl;

    //---background
    sumCount = 0;
    for (Int_t j =  min_bin_1; j <= min_bin; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;
    }
    for (Int_t j =  max_bin; j <= max_bin_1; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;
    }
    double b2_f = sumCount;
    double Area_f_psa = Area_f_psa_tot - b2_f;//subtraction of background
    //cout << "Area_FEP_psa = " << Area_f_psa << endl;
    
    //------- Accentance calculation --------
    
    double Acc_d = Area_d_psa/Area_d;
    double Acc_f = Area_f_psa/Area_f;
    double Acc_bkg = bkg_psa/bkg;
    //double Acc_bkg = b2_f/b1_f;
    double Acc_d_e = pow( ( Area_d_psa_tot + b2_d ) * ( Area_d*Area_d - 2*Area_d*Area_d_psa) + Area_d_psa*Area_d_psa*(Area_d_tot + b1_d), 0.5)/(Area_d*Area_d);
    double Acc_f_e = pow( ( Area_f_psa_tot + b2_f ) * ( Area_f*Area_f - 2*Area_f*Area_f_psa) + (Area_f_psa*Area_f_psa)*(Area_f_tot + b1_f), 0.5)/(Area_f*Area_f);
    //double Acc_bkg_e = (pow(((Acc_bkg/100)*(1-(Acc_bkg/100)))/bkg,0.5));
    double Acc_bkg_e = pow( b2_f * ( b1_f*b1_f - 2*b1_f*b2_f) + (b2_f*b2_f)*b1_f, 0.5)/(b1_f*b1_f);
    
    //-------- print the results of PSD ----------  
    
    TString result_psa_dep = Form("SF_dep = ( %4.2f +/- %4.2f ) %% ",Acc_d*100,Acc_d_e*100);
    TString result_psa_fep = Form("SF_fep = ( %4.2f +/- %4.2f ) %% ",Acc_f*100,Acc_f_e*100);
    //TString result_psa_bkg = Form("SF_bkg = ( %4.2f +/- %4.2f ) %% ",Acc_bkg*100,Acc_bkg_e*100);
    
    TPaveLabel *label1 = new TPaveLabel(0.7,0.5,0.98,0.58,result_psa_dep,"NDC");
    TPaveLabel *label2 = new TPaveLabel(0.7,0.4,0.98,0.48,result_psa_fep,"NDC");
    //TPaveLabel *label3 = new TPaveLabel(0.7,0.3,0.98,0.38,result_psa_bkg,"NDC");
    label1->Draw();
    label2->Draw();
    //label3->Draw();
    c2->Update();
    
    
    
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    //------------------- AoE correction ---------------------------
    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    
    //vector<double> AoE_corr(TH1D hist, int chn, double min_aoe_dep,  double calib, int shapingMethod, int shap){
    cout << "Correction of AoE..." << endl;
    
    //------ energy to analize -----
    double ene[20], ene_err[20], mu_AoE_[20], mu_AoE_err_[20];
    
    ene[0] = 1375.;
    ene[1] = 1425.;
    ene[2] = 1550.;
    ene[3] = 1720.;
    ene[4] = 1845.;
    ene[5] = 1900.;
    ene[6] = 1955.;
    ene[7] = 2010.;
    //ene[8] = 2065.;
    ene[8] = 2155.;
    ene[9] = 2245.;
    ene[10] = 2300.;
    ene[11] = 2355.5;
    int nPeaks = 12;
    
    TCanvas *c3 = new TCanvas("c3","AoE-correction");
    //c3->SetGrid();
    c3->Divide(4,3);

    for(int j = 0; j < nPeaks; j++){
    	//--------- find parameters from spectrum -------
    	ene_err[j]=0;
    	//energy region
    	double peak = ene[j] / calib;
    	double min = peak - 20/calib;
    	double max = peak + 20/calib;
    	int bin = ( max - min ) * 2 ;
      
      ////////////////////////////////////////////////////
    	//----------------- AoE_analysis ------------------
    	////////////////////////////////////////////////////
    	
    	c3->cd(j+1);
    	gStyle->SetPalette(1);
    	
    	TCut cut_e = Form("GEMDEnergyGauss_energy[%d] > %g && GEMDEnergyGauss_energy[%d] < %g",chn,min,chn,max);
    	
    	//--------------- draw the A/E --------------
    	
    	TH1D *h_AoE = new TH1D("h_AoE",Form("chn%d - AoE - region %d;A/E;counts",chn,j), 200, 0.7, 1.1);
    	tier2->Draw(Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g>>h_AoE",chn, chn, mean_aoe),cut_e);
    	h_AoE->GetXaxis()->CenterTitle();
    	h_AoE->GetYaxis()->CenterTitle();
    	c3->Update();
    	
    	//---------------- fit the A/E ------------
    	
    	//TF1 *f_Agaus = new TF1("f_Agaus","gaus(0)+pol1(3)",min_AoE,max_AoE);
    	
    	double fMin = 0;
    	double fMax = 3;
    	
    	TF1 *FitAoE = new TF1( "FitAoE","[0]*TMath::Exp( -TMath::Power( (x-[1])/[2],2. )/2. ) +[3]*((TMath::Exp([4]*(x-[5])))+[6])/((TMath::Exp((x-[5])/[7]))+[5])",fMin, fMax );
    	
    	double A_AoE = h_AoE->GetMaximum();
    	double mu_AoE = h_AoE->GetBinCenter( h_AoE->GetMaximumBin());
    	double sigma_AoE = mu_AoE*0.01;
    	
    	FitAoE->SetParName( 0, "Amp. - A/E" );
    	FitAoE->SetParName( 1, "Mean - A/E" );
    	FitAoE->SetParName( 2, "Sigma - A/E" );
    	//f_Agaus->SetParName( 3, "Flat Bkg" );
    	//f_Agaus->SetParName( 4, "Linear Bkg" );
    	//f_Agaus->SetParameter(0, A_AoE);
    	//f_Agaus->SetParameter(1, mu_AoE);
    	//f_Agaus->SetParameter(2, sigma_AoE);
    	FitAoE->SetParameter(0, A_AoE*0.9);
    	FitAoE->SetParLimits(0, A_AoE*0.6, A_AoE*1.2);
    	FitAoE->SetParameter(1, mu_AoE);
    	FitAoE->SetParLimits(1, 0.98, 1.02);
    	FitAoE->SetParameter(2, 0.01);
    	FitAoE->SetParLimits(2, 0.005, 0.03);
    	FitAoE->SetParameter( 3, A_AoE/3. );
    	FitAoE->SetParameter( 4, 0.1 );
    	//FitAoE->SetParLimits( 4, 0, 1);
    	FitAoE->SetParameter( 5, 0.98 );
    	//FitAoE->SetParLimits( 5, 0.7, 1. );
    	FitAoE->SetParameter( 6, -0.8 );
    	FitAoE->SetParLimits( 6, -1, 0 );
    	FitAoE->SetParameter( 7, 0.02);
    	FitAoE->SetParLimits( 7, 0.005, 0.1 );
    	
    	gStyle->SetOptFit(111);
    	
    	h_AoE->Fit("FitAoE","L","", 0.7, 1.1);
    	TF1 *gaussian = new TF1("gaussian","[0]*TMath::Exp( -TMath::Power( (x-[1])/[2],2. )/2. )",0.7,1.1);
    	gaussian->SetParameters(FitAoE->GetParameter(0),FitAoE->GetParameter(1),FitAoE->GetParameter(2));

    	TF1 *step = new TF1("step", "[0]*((TMath::Exp([1]*(x-[2])))+[3])/((TMath::Exp((x-[2])/[4]))+[2])",fMin, fMax );
    	step->SetParameters(FitAoE->GetParameter(3), FitAoE->GetParameter(4), FitAoE->GetParameter(5), FitAoE->GetParameter(6), FitAoE->GetParameter(7));

    	gaussian->SetLineColor(kAzure);
    	step->SetLineColor(kSpring+2);
    	gaussian->SetLineStyle(2);
    	step->SetLineStyle(2);
    	gaussian->Draw("same");
    	step->Draw("same");
    	c3->Update();
      
    	mu_AoE_[j] = FitAoE->GetParameter(1);
    	mu_AoE_err_[j] = FitAoE->GetParError(1);
    }
    c3->Print(Form("%s/chn%d_AoE_fits.pdf",resdir,chn));
    
    //----- AoE versus E ----
    
    TCanvas *c_3 = new TCanvas("c_3","AoE vs E");
    //c_3->SetGrid();
    
    TGraph *g_e = new TGraphErrors (nPeaks,ene,mu_AoE_,ene_err,mu_AoE_err_);
    
    g_e->SetTitle(Form("AoE vs E - chn %d",chn));
    g_e->GetXaxis()->SetTitle("Energy [keV]");
    g_e->GetYaxis()->SetTitle("AoE peak");
    g_e->GetXaxis()->CenterTitle();
    g_e->GetYaxis()->CenterTitle();
    g_e->SetMarkerColor(4);
    g_e->SetLineColor(4);
    g_e->SetMarkerStyle(21);
    g_e->Draw("AP");
    
    TF1 *pol1 = new TF1("pol1","pol1",1350,2400);
    
    gStyle->SetOptFit(111);
    
    g_e->Fit("pol1","","",1350,2400);
    c_3->Update();
    c_3->Print(Form("%s/chn%d_AoE_vs_E.pdf",resdir,chn));
    
    double par1 = pol1->GetParameter(0);
    double par2 = pol1->GetParameter(1);
    double errpar1 = pol1->GetParError(0);
    double errpar2 = pol1->GetParError(1);
    vector<double> b;
    b.push_back(par1);
    b.push_back(par2);
    b.push_back(errpar1);
    b.push_back(errpar2);
    b.push_back(min_AoE);
    //return(b);
    for (int i = 0; i < 5; i++)
    	cout << b.at(i) << endl;
    
    /////////////////////////////////////////////////////////
    //-------- Survival Fraction for SEP & 2614 ------------
    /////////////////////////////////////////////////////////
    
    TCanvas *c_4 = new TCanvas("c_4","energy");
    //c_4->SetGrid();
    //gPad->SetLogy();

    TH1D *h_sep = new TH1D("h_sep","", 100, 0.7, 1.1);;
    TH1D *h_2614 = new TH1D("h_2614","", 100, 0.7, 1.1);
    
    double peak[2] = {2104.,2614.5};
    double Acc[2] = {0.};
    double Acc_e[2] = {0.};
    
    for (int k = 0; k < 2; k++){
          
      cout << "Calculation of survival fraction for " << peak[k] << endl;
      
      min_AoE_corr = AoE_lowCut + b.at(1)*(peak[k] - 1592.5);
      
      //------- peak region ----------
      double min = (peak[k] - 45) / calib;
      double max = (peak[k] + 45) / calib;
      int bin = ( max - min ) * 2 ;
      
      
      //----------- draw the spectrum in the peak range -----------
      
      h_ene = new TH1D("h_ene","Energy; Uncalibrated Energy [chn]; counts",bin,min,max);
      h_ene->GetXaxis()->CenterTitle();
      h_ene->GetYaxis()->CenterTitle();
      c_4->cd();
      tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>h_ene",chn));
      c_4->Update();
      
      //----- fit of the energy peaks with a double-gaussian --------
      
      TF1 *f_gaus = new TF1("f_gaus","gaus(0)+pol0(3)",min,max);
      
      f_gaus->SetParName( 0, "Amp." );
      f_gaus->SetParName( 1, "Mean" );
      f_gaus->SetParName( 2, "Sigma" );
      f_gaus->SetParName( 3, "Flat Bkg" );
      //f_gaus->SetParName( 4, "Linear Bkg" );
      
      //--------- initialization of parameters of the fit ----------
      
      double A = h_ene->GetMaximum();
      double mu = h_ene->GetBinCenter( h_ene->GetMaximumBin());
      
      f_gaus->SetParameter(0,A);
      f_gaus->SetParameter(1,mu);
      f_gaus->SetParameter(2,mu*0.0005);
      f_gaus->SetParLimits(2,0,100);
      f_gaus->SetParameter(3,10);
      
      gStyle->SetOptFit(111);
      
      h_ene->Fit("f_gaus","L","",min,max);
      c_4->Update();
      
      double A_ = f_gaus->GetParameter(0);
      double mu_ = f_gaus->GetParameter(1);//mu_err=f_gaus->GetParError(1);
      double sigma_ = f_gaus->GetParameter(2);
      
      double n_sigma = 4.5;//number of sigma to calculate area
      
      //----draw the A/E----
      
      TCut cut_e = Form("GEMDEnergyGauss_energy[%d] > %g && GEMDEnergyGauss_energy[%d] < %g",chn, mu_ - n_sigma*sigma_, chn, mu_ + n_sigma*sigma_);
      
      c1->cd();
      if (k == 0){
	tier2->Draw(Form("((GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g)>>h_sep", chn, chn, mean_aoe), cut_e, "same");
	h_sep->Scale(1./h_sep->GetBinContent(h_sep->GetMaximumBin()));
	h_sep->SetLineColor(3);
	h_sep->SetLineWidth(2);
	h_sep->Draw("same");
    	}
      if (k == 1){
	tier2->Draw(Form("((GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g)>>h_2614", chn, chn, mean_aoe), cut_e, "same");
	h_2614->Scale(1./h_2614->GetBinContent( h_2614->GetMaximumBin()));
	h_2614->SetLineColor(2);
	h_2614->SetLineWidth(2);
	h_2614->Draw("same");
      }
      c1->Update();
      //------------ calculation of Area of peak -----------
      
      int count, sumCount = 0;
      
      //----- peak Area ------
      int min_bin = h_ene->GetXaxis()->FindBin(mu_ - n_sigma * sigma_);
      int max_bin = h_ene->GetXaxis()->FindBin(mu_ + n_sigma * sigma_);
      int min_bin_1 = h_ene->GetXaxis()->FindBin(mu_ - 2*n_sigma * sigma_);
      int max_bin_1 = h_ene->GetXaxis()->FindBin(mu_ + 2*n_sigma * sigma_);
      int delta = max_bin - min_bin;
      
      for (Int_t j = min_bin; j <= max_bin; j++){
	count = h_ene->GetBinContent(j);
	sumCount += count;
      }
      double Area_tot = sumCount;//total area
      //cout << "Area_tot = " << Area_tot << endl;
      
      //---background
      sumCount = 0;
      for (Int_t j =  min_bin_1; j <= min_bin; j++){
	count = h_ene->GetBinContent(j);
	sumCount += count;}
      for (Int_t j =  max_bin; j <= max_bin_1; j++){
	count = h_ene->GetBinContent(j);
	sumCount += count;}
      double b1 = sumCount;
      double Area = Area_tot - b1;//subtraction of background
      //cout << "Area = " << Area << endl;
      
      ////////////////////////////////////////////////////
      //---------------- PSA analysis --------------------
      
      hene_psa = new TH1D("hene_psa","Energy Spectrum; energy (FADC chn); counts",bin,min,max);
      c_4->cd();
      //	TCut cut_aoe;
      cut_aoe = Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g > %g && (GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g < 1.07", chn, chn, mean_aoe, min_AoE_corr, chn, chn, mean_aoe);
      
      tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>hene_psa",chn),cut_aoe,"sames");
      
      hene_psa->GetXaxis()->CenterTitle();
      hene_psa->GetYaxis()->CenterTitle();
      
      c_4->Update();
      
      //------------- fit DEP region after PSD ----------------
      
      f_gaus->SetParameter(0, A_*0.3);
      f_gaus->SetParameter(1, mu_);
      f_gaus->SetParameter(2, sigma_);
      f_gaus->SetParLimits(2, 0, 100);
      f_gaus->SetParameter(3, 1);
      
      gStyle->SetOptFit(111);
      
      hene_psa->Fit("f_gaus","L","sames",min,max);
      
      //-------------- get the results from fit --------------
      
      double mu_psa = f_gaus->GetParameter(1);
      double sigma_psa = f_gaus->GetParameter(2);
      
      //------------ calculation of Area of peak -----------
      
      //----- area -----
      sumCount = 0;
      min_bin = hene_psa->GetXaxis()->FindBin(mu_ - n_sigma * sigma_);
      max_bin = hene_psa->GetXaxis()->FindBin(mu_ + n_sigma * sigma_);
      min_bin_1 = hene_psa->GetXaxis()->FindBin(mu_ - 2*n_sigma * sigma_);
      max_bin_1 = hene_psa->GetXaxis()->FindBin(mu_ + 2*n_sigma * sigma_);
      delta = max_bin - min_bin;
      
      for (Int_t j = min_bin; j < max_bin + 1; j++){
	count = hene_psa->GetBinContent(j);
	sumCount += count;}
      double Area_psa = sumCount;//total area
      //cout << "Area_tot_psa = " << Area_psa << endl;
      
      //---background
      sumCount = 0;
      for (Int_t j =  min_bin_1; j <= min_bin; j++){
	count = hene_psa->GetBinContent(j);
	sumCount += count;
      }
      for (Int_t j =  max_bin; j <= max_bin_1; j++){
	count = hene_psa->GetBinContent(j);
	sumCount += count;
      }
      Area_psa -= sumCount;//subtraction of background
      //cout << "Area_psa = " << Area_psa << endl;
      
      //------- Accentance calculation --------
      
      Acc[k] = 1.*Area_psa/Area;
      //double Acc_e = Acc * ( pow(Area_psa,0.5)/Area_psa + pow(Area,0.5)/Area);
      Acc_e[k] = pow(((Area-Area_psa)*Area_psa)/pow(Area,3.),0.5);
    	
      //-------- print the results of PSD ----------  
      
      TString result_psa = Form("SF_peak = ( %4.2f +/- %4.2f ) %%",Acc[k]*100,Acc_e[k]*100);
      
      TPaveLabel *label1 = new TPaveLabel(0.7,0.5,0.98,0.58,result_psa,"NDC");
      c_4->cd();
      label1->Draw();
      c_4->Update();
      
      //-------- save plots --------
      
      c_4->Print(Form("%s/chn%d_PSA_peak%g.pdf",resdir,chn,peak[k]));
      
    }
    
    //////////////////////////////////////////////////////
    //----------- survival fraction background -----------
    //////////////////////////////////////////////////////
    
    double qbb = 2039.;
    cout << "Calculation of survival fraction background at " << qbb << endl;
    
    min_AoE_corr = AoE_lowCut + b.at(1)*(qbb-1592.5);
    
    //------- qbb region ----------
    min = (qbb - 35) / calib;
    max = (qbb + 35) / calib;
    bin = ( max - min ) * 2 ;

    cut_e = Form("GEMDEnergyGauss_energy[%d] > %g && GEMDEnergyGauss_energy[%d] < %g",chn, min, chn, max);

    //----------- draw the A/E -----------
    TH1D *h_bkg = new TH1D("h_bkg","", 100, 0.7, 1.1);
    
    c1->cd();
    tier2->Draw(Form("((GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g)>>h_bkg", chn, chn, mean_aoe), cut_e, "same");
    h_bkg->Scale(1./h_bkg->GetBinContent(h_bkg->GetMaximumBin()));
    h_bkg->SetLineColor(1);
    h_bkg->SetLineWidth(2);
    h_bkg->Draw("same");
    c1->Update();
    
    //----------- draw the spectrum in the peak range -----------
    
    h_ene = new TH1D("h_ene","Energy; Uncalibrated Energy [chn]; counts",bin,min,max);
    h_ene->GetXaxis()->CenterTitle();
    h_ene->GetYaxis()->CenterTitle();
    c_4->cd();
    tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>h_ene",chn));
    c_4->Update();
    
    //------------ calculation of Area -----------
    
    sumCount = 0;
          
    //----- peak Area ------
    
    for (Int_t j = min; j <= max; j++){
    	count = h_ene->GetBinContent(j);
    	sumCount += count;}
    double b2039 = sumCount;//bkg
    
    //---------------- PSA --------------------
    
    hene_psa = new TH1D("hene_psa","Energy Spectrum; energy (FADC chn); counts",bin,min,max);
    
    cut_aoe = Form("(GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g > %g && (GEMDCurrentPSA_A[%d]/GEMDEnergyGauss_energy[%d])/%g < 1.07", chn, chn, mean_aoe, min_AoE_corr, chn, chn, mean_aoe);
    c_4->cd();
    tier2->Draw(Form("GEMDEnergyGauss_energy[%d]>>hene_psa",chn),cut_aoe,"sames");
    
    hene_psa->GetXaxis()->CenterTitle();
    hene_psa->GetYaxis()->CenterTitle();
    c_4->Update();
    
    //------------ calculation of Area -----------
    //----- area -----
    sumCount = 0;
    for (Int_t j = min; j < max; j++){
    	count = hene_psa->GetBinContent(j);
    	sumCount += count;
    }
    double b2039_psa = sumCount;//total area
    //cout << "Area_tot_psa = " << bkg_psa << endl;
    
    //------- Accentance calculation --------
    
    double SF_bkg = 1.*b2039_psa/b2039;
    double SF_bkg_err = pow(((b2039 - b2039_psa)*b2039_psa)/pow(b2039, 3.), 0.5);
    
    //-------- print the results of PSD ----------  
    
    TString result_psa = Form("SF_bkg = ( %4.2f +/- %4.2f ) %%",SF_bkg*100, SF_bkg_err*100);
    c_4->cd();
    TPaveLabel *label_b = new TPaveLabel(0.7,0.5,0.98,0.58,result_psa,"NDC");
    label_b->Draw();
    c_4->Update();
    
    //-------- save plots --------
    c1->Print(Form("%s/chn%d_AoE_DEP.pdf",resdir,chn));
    c2->Print(Form("%s/chn%d_PSA_DEP-FEP.pdf",resdir,chn));
    //c_4->Print(Form("%s/chn%d_PSA_bkg%g.pdf",resdir,chn,peak));
    
    //-------- Save results in file --------
    
    ofstream file;
    char *out_file = Form("%s/PSA_results.txt",resdir);
    char *out_results = Form("%d %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n",chn, FWHM_AoE*100, FWHM_AoE_err*100, Acc_d*100, Acc_d_e*100, Acc_f*100, Acc_f_e*100,Acc[0]*100, Acc_e[0]*100,Acc[1]*100, Acc_e[1]*100, SF_bkg*100, SF_bkg_err*100);
    
    file.open(out_file,std::ios_base::app|std::ios_base::out);
    file << out_results;
    file.close();
  }
  //myapp->Run();
  return 0;
}
