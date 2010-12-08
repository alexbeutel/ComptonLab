
#include<utility>

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>

void CTplotData0(){
  gStyle->SetOptFit(1111);

  // read in data
  vector<Double_t> vx;
  vector<Double_t> vy;
  Double_t xdat,ydat;

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextFont(1);
  text->SetTextColor(kRed);

  char tex[100];

  TArrow * arrow = new TArrow();
  arrow->SetLineColor(2);
  arrow->SetFillColor(2);

  // Read in Data
  fstream infile;
//  infile.open("deg20.TKA", ios_base::in);
  infile.open("200sec-0deg-without.TKA", ios_base::in);
  while (infile>>ydat){
    //cout <<ydat<<"\n"; 
    vy.push_back(ydat) ;
  }
  infile.close();

  Int_t vsize = vy.size();
  //cout <<" Here I am "<<vsize<<"\n";



  // Read in Background
  fstream infile;
  infile.open("200sec-0deg-without.TKA", ios_base::in);
  while (infile>>xdat){
    //cout <<xdat<<"\n"; 
    vx.push_back(xdat) ;
  }
  infile.close();
  Int_t vsize = vy.size();
  //cout <<" Here I am "<<vsize<<"\n";


 
  // histogram parameters
  Int_t ntbin = 16384;       // number of total bins in the data
  Int_t nbin = 200;       // number of  bins in the final plot

  // book histogram
  TH1F *hist1 = new TH1F("hist1","Signal",          nbin,0.0,16384.0);
  TH1F *hist2 = new TH1F("hist2","Background",      nbin,0.0,16384.0);
  TH1F *hist3 = new TH1F("hist3","  Sig - Background",nbin,0.0,16384.0);

  Double_t BackgroundTimeRatio = 200.0/60.0;

  //plot histogram
  for (Int_t i=0; i!=ntbin; i++)
    { 
  Double_t xx=i+0.5;
  hist1->Fill(xx,vy[i]);
  hist2->Fill(xx,vx[i]*BackgroundTimeRatio);
  //hist2->Fill(xx,0);
    }

  hist1->Sumw2();
  hist2->Sumw2();
  hist3->Sumw2();

  hist1->SetStats(kFALSE);
  hist1->GetXaxis()->SetTitle("X axis title");
  hist1->GetYaxis()->SetTitle("Y axis title");
  TCanvas *myc1 =new TCanvas("myc1","Sig");
  hist1->Draw("e");

  hist2->SetStats(kFALSE);
  hist2->GetXaxis()->SetTitle("X axis title");
  hist2->GetYaxis()->SetTitle("Y axis title");
  TCanvas *myc2 =new TCanvas("myc2","Back");
  hist2->Draw("e");

  //Get the signal minus background
  hist3 ->Add(hist1,hist2,1.0,-1.0);
  TCanvas *myc3 =new TCanvas("myc3","Final");
  hist3->Draw();

  //fitting
  TCanvas *myc4 =new TCanvas("myc4","Sig Fit");

  Double_t par[7];  //par[0]= peak, par[1]=mean, par[2]=width, par[3]->par[6] for poly

  //define fit
//  TF1 *fit1 = new TF1("fit1",GaussPoly, 4000.,10000.,7); // range & number of parameters 
  TF1 *fit1 = new TF1("fit1",GaussPoly, 4000.,10000.,7); // range & number of parameters 
  TF1 *fit2 = new TF1("fit2",GaussPoly, 4000.,10000.,7); // range & number of parameters 
  TF1 *fit3 = new TF1("fit3",GaussPoly, 4000.,10000.,7); // range & number of parameters 
  fit1 ->SetParameters(500., 8500., 300.0, 100.0, 0.0, 0.0, 0.0); 
  fit1->SetLineWidth(1);
  fit1->SetLineColor(kGreen);	
  fit2 ->SetParameters(500., 8500., 300.0, 100.0, 0.0, 0.0, 0.0); 
  fit2->SetLineWidth(1);
  fit2->SetLineColor(kGreen);	
  fit3 ->SetParameters(500., 8500., 300.0, 100.0, 0.0, 0.0, 0.0); 
  fit3->SetLineWidth(1);
  fit3->SetLineColor(kGreen);	
  Double_t xl1=5000.;
  Double_t xh1=9500.;
//  hist3  ->Fit("fit1","R"," ",xl1,xh1);
//  hist3->Draw();
  hist1->Fit("fit1","R"," ",xl1,xh1);
  hist1->Draw();
  fit1->GetParameters(&par[0]);

  Int_t LowerWindowBin = 90-15;
  Int_t UpperWindowBin = 90;

  Double_t MeasuredSignal = hist1->Integral(LowerWindowBin,UpperWindowBin);

  windowLowerBound = hist1->GetBinLowEdge(LowerWindowBin);
  windowUpperBound = hist1->GetBinLowEdge(UpperWindowBin)
    + hist1->GetBinWidth(UpperWindowBin);
  
  arrow->DrawArrow(windowLowerBound, 0.0, windowLowerBound, (8000),0.005,"<");
  arrow->DrawArrow(windowUpperBound, 0.0, windowUpperBound, (8000),0.005,"<");
  arrow->Draw("same");

  Double_t peakHeight=par[0];
  Double_t peakChannel=par[1];
  Double_t peakWidth=par[2]; 

  text->SetTextSize(0.045);
  sprintf(tex,"peak channel: %1d \n",peakChannel);
  text->DrawLatex(0.15, 0.7, tex);
  sprintf(tex,"peak width: %1d \n",peakWidth);
  text->DrawLatex(0.15, 0.75, tex);
  sprintf(tex,"Signal in Peak: %1d \n", MeasuredSignal);
  text->DrawLatex(0.15, 0.8, tex);

  TCanvas *myc5 =new TCanvas("myc5","Bk Fit");

  hist2->Fit("fit2","R"," ",xl1,xh1);
  hist2->Draw();
  fit2->GetParameters(&par[0]);

  Double_t MeasuredSignal = hist2->Integral(LowerWindowBin,UpperWindowBin);
  arrow->DrawArrow(windowLowerBound, 0.0, windowLowerBound, (8000),0.005,"<");
  arrow->DrawArrow(windowUpperBound, 0.0, windowUpperBound, (8000),0.005,"<");
  arrow->Draw("same");

  Double_t peakHeight=par[0];
  Double_t peakChannel=par[1];
  Double_t peakWidth=par[2]; 

  text->SetTextSize(0.045);
  sprintf(tex,"peak channel: %1d \n",peakChannel);
  text->DrawLatex(0.15, 0.7, tex);
  sprintf(tex,"peak width: %1d \n",peakWidth);
  text->DrawLatex(0.15, 0.75, tex);
  sprintf(tex,"Signal in Peak: %1d \n", MeasuredSignal);
  text->DrawLatex(0.15, 0.8, tex);

  TCanvas *myc6 =new TCanvas("myc6","Sig-Bk Fit");

  hist3->Fit("fit3","R"," ",xl1,xh1);
//  hist3->SetMaximum(5000);
  hist3->Draw();

  LowerWindowBin = 90;
  UpperWindowBin = 111;

  MeasuredSignal = hist3->Integral(LowerWindowBin,UpperWindowBin);
  
  arrow->DrawArrow(windowLowerBound, 0.0, windowLowerBound, (8000),0.005,"<");
  arrow->DrawArrow(windowUpperBound, 0.0, windowUpperBound, (8000),0.005,"<");
  arrow->Draw("same");

  fit3->GetParameters(&par[0]);

  Double_t peakHeight=par[0];
  Double_t peakChannel=par[1];
  Double_t peakWidth=par[2]; 

  text->SetTextSize(0.045);
  sprintf(tex,"peak channel: %1d \n",peakChannel);
  text->DrawLatex(0.15, 0.7, tex);
  sprintf(tex,"peak width: %1d \n",peakWidth);
  text->DrawLatex(0.15, 0.75, tex);
  sprintf(tex,"Signal in Peak: %1d \n", MeasuredSignal);
  text->DrawLatex(0.15, 0.8, tex);
}

Double_t Pol2(Double_t *x, Double_t *par)
//Polynomial 
{ 
  Double_t arg1= 0; 
  arg1=par[0]+ par[1]*x[0]+ par[2]*x[0]*x[0] ; 
  Double_t fitval = arg1 ; 
  /* 
  cout <<par[0]<<" "<<par[1]<<" "<<par[2]<<"\n"; 
  cout <<par[3]<<" "<<par[4]<<" "<<par[5]<<"\n"; 
  cout <<x[0]<<" "<<arg1<<" "<<arg2<<" "<<arg3<<" "<<fitval<<" "<<"\n"; 
  */ 
  return fitval; 
}

Double_t DGauss(Double_t *x, Double_t *par)  
//Derivative of Gauss 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*(-xnew)/par[2]/par[2];
    fitval=fitval*exp(-xnew*xnew/2.0/par[2]/par[2]);
    return fitval;
}

Double_t DLorentz(Double_t *x, Double_t *par) 
//Derivative of Lorentz
{
  Double_t amp=par[0];
  Double_t avg=par[1];
  Double_t halfs=(par[2]/2.)*(par[2]/2.);
  Double_t fitval=2.0*amp*halfs/( pow( (x[0]-avg)*(x[0]-avg)+halfs, 2.0));
  fitval=(x[0]-avg)*fitval;
  return fitval;
}

Double_t DGauss(Double_t *x, Double_t *par)  
//Gauss fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    return fitval;
}

Double_t GaussPoly(Double_t *x, Double_t *par)  
//Gauss+Pol3 fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    fitval=fitval+par[3]+par[4]*xnew+par[5]*xnew*xnew+par[6]*xnew*xnew*xnew;
    return fitval;
}
