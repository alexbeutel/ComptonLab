
#include<utility>

#include <iostream>
#include <fstream>
//#include <vector>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>

void CalibrationFit(){
  gStyle->SetOptFit(0000);

  Int_t nbins = 7;
/*
  Double_t energy[7] = {88.0336, 122.1,  356.0,  511.0,  662.0,  835.0,   1173.0};
  Double_t energyError[7] = {0.0,     0.0,    0.0,    0.0,    0.0,    0.0,     0.0};
  Double_t channel[7] = {1206.0,  1690.0, 4778.0, 6691.0, 8514.0, 10616.0, 14798.0};
  Double_t channelError[7] = {79.0,    66.0,   176.0,  214.0,  241.0,  266.0,   319.0}; 
*/

  Double_t energy[7] = {0.384, 0.356,1.275,0.511,0.662,1.333,1.173};
  Double_t energyError[7] = {0.0,     0.0,    0.0,    0.0,    0.0,    0.0,     0.0};
  Double_t channel[7] = {4209.0, 3608.0, 14310.0,5956.0, 7644.0, 14880.0, 13220.0};
  Double_t channelError[7] = {1.7,9.7,5.387,1.6,1.5,12.94,12.79};

  // histogram parameters
  Int_t ntbin = 16384;       // number of total bins in the data
  Int_t nbin = 16384;       // number of  bins in the final plot

  TCanvas *myc4 =new TCanvas("myc4","Fit");

  // book histogram
  TGraphErrors *graphCalib = new TGraphErrors(nbins, energy, channel, energyError, channelError);

	/*
  graphCalib->SetXTitle("Energy");
  graphCalib->SetYTitle("Channel Number");
	*/

  //fitting
  Double_t par[1];  //par[0]= peak, par[1]=mean, par[2]=width, par[3]->par[6] for poly

  //define fit
//  TF1 *fit1 = new TF1("fit1",GaussPoly, 4000.,10000.,7); // range & number of parameters 
  TF1 *fit1 = new TF1("fit1",Pol1, 2000.,8000.,7); // range & number of parameters 
  fit1 ->SetParameters(0.0, .1); 
  fit1->SetLineWidth(1);
  fit1->SetLineColor(kGreen);	
  Double_t xl1=0.;
  Double_t xh1=1300;
//  hist3  ->Fit("fit1","R"," ",xl1,xh1);
//  hist3->Draw();
  graphCalib->Fit(fit1,"R"," ",xl1,xh1);
  graphCalib->SetMarkerColor(4);
  graphCalib->SetLineColor(4);
  graphCalib->SetMarkerStyle(20);
  graphCalib->SetMarkerSize(0.5);

  Double_t Intercept = fit1->GetParameter(0);
  Double_t Slope     = fit1->GetParameter(1);

  graphCalib->Draw("AP");

  /*
  graphCalib->GetParameters(&par[0]);

  Double_t Intercept=par[0];
  Double_t Slope=par[1];
  */

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextFont(1);
  text->SetTextColor(kRed);

  char tex[100];
  
  text->SetTextSize(0.045);
  sprintf(tex,"Intercept: %2.1f \n",Intercept);
  text->DrawLatex(0.15, 0.7, tex);
  sprintf(tex,"Slope: %3.1f \n",Slope);
  text->DrawLatex(0.15, 0.8, tex);
 

  }

Double_t Pol1(Double_t *x, Double_t *par)
//Polynomial 
{ 
  Double_t arg1= 0; 
  arg1=par[0]+ par[1]*x[0] ; 
  Double_t fitval = arg1 ; 
  /* 
  cout <<par[0]<<" "<<par[1]<<" "<<par[2]<<"\n"; 
  cout <<par[3]<<" "<<par[4]<<" "<<par[5]<<"\n"; 
  cout <<x[0]<<" "<<arg1<<" "<<arg2<<" "<<arg3<<" "<<fitval<<" "<<"\n"; 
  */ 
  return fitval; 
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
