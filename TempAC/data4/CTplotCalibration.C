
#include <iostream>
#include <fstream>
#include <vector>

void CTplotCalibration(const char * mainFile = "", const char * bgFile = "", 
		Double_t a = 0, Double_t b = 0, Double_t c = 0, Double_t d = 0, Double_t e = 0){
	//return;
	// read in data
	vector<Double_t> vx;
	vector<Double_t> vy;
	Double_t xdat,ydat;

	// Read in Data
	fstream infile;
	//  infile.open("deg20.TKA", ios_base::in);
	//  infile.open("10-26-calibration-cs137.TKA", ios_base::in); 
	infile.open(mainFile, ios_base::in);
	while (infile>>ydat){
		//cout <<ydat<<"\n"; 
		vy.push_back(ydat) ;
	}
	infile.close();

	Int_t vsize = vy.size();
	//cout <<" Here I am "<<vsize<<"\n";



	// Read in Background
	fstream infile;
	//  infile.open("10-26-calibratration-bg.TKA", ios_base::in);
	infile.open(bgFile, ios_base::in);
	while (infile>>xdat){
		//cout <<xdat<<"\n"; 
		vx.push_back(xdat) ;
	}
	infile.close();
	Int_t vsize = vy.size();



	// histogram parameters
	Int_t ntbin = 16384;       // number of total bins in the data
	Int_t nbin = 400;       // number of  bins in the final plot

	// book histogram
	TH1F *hist1 = new TH1F("hist1","Signal",          nbin,0.0,16384.0);
	TH1F *hist2 = new TH1F("hist2","Background",      nbin,0.0,16384.0);
	TH1F *hist3 = new TH1F("hist3","  Sig - Background",nbin,0.0,16384.0);
	hist3 ->Sumw2();  //take care of error properly

	//plot histogram
	for (Int_t i=0; i!=ntbin; i++)
	{ 
		Double_t xx=i+0.5;
		hist1->Fill(xx,vy[i]);
		hist2->Fill(xx,vx[i]);
		//hist2->Fill(xx,0);
	}
	hist1->SetStats(kFALSE);
	hist1->GetXaxis()->SetTitle("X axis title");
	hist1->GetYaxis()->SetTitle("Y axis title");
	//amb79
	TCanvas *myc1 =new TCanvas("myc1","Sig");
	hist1->Draw("e");

	hist2->SetStats(kFALSE);
	hist2->GetXaxis()->SetTitle("X axis title");
	hist2->GetYaxis()->SetTitle("Y axis title");
	TCanvas *myc2 =new TCanvas("myc2","Back");
	hist2->Draw("e");

	//Get the signal minus background
	hist3 ->Add(hist1,hist2,1.0,-1.0);
	//amb79
	TCanvas *myc3 =new TCanvas("myc3","Final");
	hist3->Draw();


	//fitting
	Double_t par[7];  //par[0]= peak, par[1]=mean, par[2]=width, par[3]->par[6] for poly
	//define fit
	TF1 *fit1 = new TF1("fit1",GaussPoly, a,b,7); // range & number of parameters 
	//  TF1 *fit1 = new TF1("fit1",GaussPoly, 1000.,16000.,7); // range & number of parameters 
	//  fit1 ->SetParameters(500., 6000., 300.0, 100.0, 0.0, 0.0, 0.0);  
	fit1 ->SetParameters(500., c, 300.0, 100.0, 0.0, 0.0, 0.0); 

	Double_t xl1=d; //3000.;
	Double_t xh1=e; //9000.;

	hist3->Fit("fit1","R"," ",xl1,xh1);
	hist3->Draw();

	//hist1->Fit("fit1","R"," ",xl1,xh1);
	//hist1->Draw();

	//hist2->Fit("fit1","R"," ",xl1,xh1);
	//hist2->Draw();
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
