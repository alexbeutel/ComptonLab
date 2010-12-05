#include <iostream>
#include <fstream>
#include <vector>
#include "TF1.h"

void CTfit1(){
  //fitting
  Double_t par[6];  //par[0]= peak, par[1]=mean, par[2]=width, par[3]->par[6] for poly

  //define fit
  TF1 *fit1 = new TF1("fit1",GaussPoly2, 1000.,16000.,7); // range & number of parameters 
  fit1 ->SetParameters(500., 8000., 300.0, 100.0, 0.0, 0.0, 0.0); 
  Double_t xl1=3000.;
  Double_t xh1=12000.;
  
  hist3  ->Fit("fit1","R"," ",xl1,xh1);
  gStyle->SetOptFit(1111);
  hist3->Draw();
  Double_t par[7];
  fit1->GetParameters(&par[0]);
  //TF1 *myfit = hist3->GetFunction("fit1");
  //Double_t p0 = myfit->GetParamter(0);
  //polynom3 = new TF1("polynom3","p3+p4*x+p5*x*x+p6*x*x*x",xl1,xh1);
  //polynom3->Draw();  
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

Double_t Gauss(Double_t *x, Double_t *par)  
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

Double_t GaussPoly2(Double_t *x, Double_t *par)  
//Gauss+Pol2 fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    fitval=fitval+par[3]+par[4]*xnew+par[5]*xnew*xnew;
    return fitval;
}
