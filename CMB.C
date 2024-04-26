#include <fstream>
#include <iostream>
#include <array>
#include <vector>
#include <type_traits>
#include "RooGenericPdf.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
using namespace std;
using namespace RooFit;

using Row = array<double, 5>;

double fitFunc(double *x, double *par) {
   return 1.05457e-34 * pow(x[0],3) / (pow(TMath::Pi(),2)*pow(3e8,3) * (exp(x[0]*1.05457e-34/(1.38e-23*par[0])) - 1));
}

void CMB(){
  vector < Row > v;

  ifstream ifs("cmbFile.txt");

  if (ifs)
  {
    string line;
    getline(ifs, line); // read and discard the header
  }
  for (;;)
  {
    Row row;
    for (auto & d : row)
    {
      ifs >> d; // read all the items of the row
    }
    if ( !ifs) break;
    v.push_back(row); // append the row to the vector
  }

  vector<double> col1; 
  vector<double> col2;
  vector<double> col3;
  vector<double> col4;
  vector<double> col5;
  
  int rowNum = v.size();
  double col1ar[rowNum], col2ar[rowNum], col4ar[rowNum], omega_data[rowNum], u_w_data[rowNum];

//Now I will define my important constants:
        double c = 3e8; // m/s
        double hbar = 1.05457e-34; // #J/s
        double k_B = 1.38e-23; // #J/K
        double c1 = 2*TMath::Pi()*c*100;
        double c2 = (2/c)*10e-21;

        //double omega_data = transform(col1.begin(), col1.end(), col1.begin(), _1 * const1); 

  for(int i = 0; i < v.size(); i++){
	//cout << v[0][i] << endl;
	col1.push_back(c1*v[i][0]);
        col2.push_back(c2*v[i][1]);
        col3.push_back(v[i][2]);
        col4.push_back(v[i][3]);
        col5.push_back(v[i][4]);

	col1ar[i] = c1*v[i][0];
	col2ar[i] = c2*v[i][1];
	col4ar[i] = c2*v[i][3];
  }
	//omega_data = c1*col1ar;
	//u_w_data = c2*col2ar;
	//cout << col1ar[0] << endl;
	
	//PLOTTING

	 TF1 *fitFcn = new TF1("fitFcn",fitFunc,300e9,4.1e12,1);
   	 fitFcn->SetNpx(500);
   	 fitFcn->SetLineWidth(3);
   	 fitFcn->SetLineColor(kRed);
 

	TCanvas *can = new TCanvas("c","Black Body Spectrum Fitted Plot");
	can->SetGrid();
	
	double total = 1.0;
	double p1size = 0.35;
	double p2size = total - p1size;

	TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
        pad->SetTopMargin(0.9); //.7
        pad->SetBottomMargin(p1size);
	pad->Draw();
        pad->SetFillStyle(0);
	pad->SetGrid();
        pad->cd();

	TGraphErrors *g=new TGraphErrors(rowNum, col1ar, col2ar, nullptr, col4ar);
	
	fitFcn->SetParameters(1);
   	g->Fit("fitFcn","0");
	
	g->SetTitle("Black Body Spectrum;#omega; U(#omega)");
	g->SetMarkerSize(3);
	g->Draw("ap"); // ap good and using "AC*" looks really nice hehe
	fitFcn->Draw("same");
	
	TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., 1.);
        pad2->SetTopMargin(p2size);
        pad2->Draw();
	pad2->SetGrid();
        pad2->SetFillStyle(0);
        pad2->cd();

	TGraphErrors* residual = new TGraphErrors;
  for(Int_t i = 1; i <= rowNum; i++){
    Double_t x = col1ar[i-1]; //i*100;
    Double_t y0 = fitFcn->Eval(x);
    Double_t y = (g->GetY()[i - 1] - y0)/TMath::Sqrt(y0);
    Double_t ey = g->GetErrorY(i - 1)/TMath::Sqrt(y0);
    residual->SetPoint(i - 1, x, y);
    residual->SetPointError(i - 1, 0, ey);
  } // i
 
  residual->GetHistogram()->SetMaximum(.15e-11);
  residual->GetHistogram()->SetMinimum(-.15e-11);
  residual->SetLineWidth(1);
  residual->GetHistogram()->GetYaxis()->SetNdivisions(4);
  //residual->GetHistogram()->GetXaxis()->SetNdivisions(0);
  residual->SetMarkerSize(3);
  residual->Draw("ap");
  residual->SetTitle("Black Body Spectrum;#omega; Residuals");
	//can2->Draw("ap");

//LEGEND:
	TLegend *legend = new TLegend(0.9,0.8,0.8,0.9);
	//legend->SetHeader("Legend","C");
	legend->AddEntry(g,"Data","l");
	legend->AddEntry(fitFcn,"Predicted","l");
	legend->Draw();

  //return 0;
}

