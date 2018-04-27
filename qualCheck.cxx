#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TVectorT.h"
#include "TH2.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCut.h"

using namespace std;

int main( int argc, char ** argv){
	if ( argc < 1 ){
		cerr << "Wrong number of arguments. Instead use\n"
		     << "\t./qualCheck \n";
		return -1;
	}
	
	// Files to load in order to compare last and uptolast shift data
	TFile * inFile_3He_uptolast = new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/3He_slow_uptolast.root");
	TFile * inFile_3H_uptolast = new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/3H_slow_uptolast.root");

	TFile * inFile_3He_last = new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/3He_slow_last.root");
	TFile * inFile_3H_last = new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/3H_slow_last.root");

	// Now with the files, load the TTrees that we will access in order
	// to make 3 histograms for each tree
	TTree * inTree_3He_uptolast = (TTree *) inFile_3He_uptolast->Get("sk");
	TTree * inTree_3H_uptolast = (TTree *) inFile_3H_uptolast->Get("sk");
	TTree * inTree_3He_last = (TTree *) inFile_3He_last->Get("sk");
	TTree * inTree_3H_last = (TTree *) inFile_3H_last->Get("sk");

	// Now for all of these trees, ready the histograms that we will fill in
	// order to make canvases to compare
	TH2D * vert_3He_uptolast = new TH2D("vert_3He_uptolast","Vertex Correlation: Up to Last Shift",30,-0.15,0.15,30,-0.15,0.15);
	TH2D * vert_3H_uptolast = new TH2D("vert_3H_uptolast","Vertex Correlation: Up to Last Shift",30,-0.15,0.15,30,-0.15,0.15);
	TH2D * vert_3He_last = new TH2D("vert_3He_last","Vertex Correlation: Last Shift",30,-0.15,0.15,30,-0.15,0.15);
	TH2D * vert_3H_last = new TH2D("vert_3H_last","Vertex Correlation: Last Shift",30,-0.15,0.15,30,-0.15,0.15);
	
	double Lmom, Rmom;
	TVectorT <double> *kin_3H = (TVectorT <double> *) inFile_3H_uptolast->Get("kinematics");
	TVectorT <double> *kin_3He = (TVectorT <double> *) inFile_3He_uptolast->Get("kinematics");
	Lmom = floorf( (*kin_3H)[1] * 100 ) / 100;
	Rmom = floorf( (*kin_3H)[3] * 100 ) / 100;
	TH2D * mom_3H_uptolast = new TH2D("mom_3H_uptolast","Momentum Correlation: Up to Last Shift",40,Rmom-0.2,Rmom+0.2,40,Lmom-0.2,Lmom+0.2);
	TH2D * mom_3H_last = new TH2D("mom_3H_last","Momentum Correlation: Last Shift",40,Rmom-0.2,Rmom+0.2,40,Lmom-0.2,Lmom+0.2);
	Lmom = floorf( (*kin_3He)[1] * 100 ) / 100;
	Rmom = floorf( (*kin_3He)[3] * 100 ) / 100;
	TH2D * mom_3He_uptolast = new TH2D("mom_3He_uptolast","Momentum Correlation: Up to Last Shift",40,Rmom-0.2,Rmom+0.2,40,Lmom-0.2,Lmom+0.2);
	TH2D * mom_3He_last = new TH2D("mom_3He_last","Momentum Correlation: Last Shift",40,Rmom-0.2,Rmom+0.2,40,Lmom-0.2,Lmom+0.2);
	

	// Cuts we will apply
	TCut zDif = "TMath::Abs(L_z-R_z)<0.02";
	TCut Ldelta = "TMath::Abs(L_delta)<.045";
	TCut Rdelta = "TMath::Abs(R_delta)<.045";
	TCut Rcoll  = "((R_ytar+1.5*R_yptar)/0.08)**2 + ((1.5*R_xptar)/0.08)**2 <= 1";
	TCut Lcoll  = "((L_ytar+1.5*L_yptar)/0.08)**2 + ((1.5*L_xptar)/0.08)**2 <= 1";
	TCut z_tar  = "TMath::Abs(L_z)<0.09";
	TCut ePID = "L_EoverP>0.5";
	TCut thrq = "thrq*TMath::RadToDeg()<40.";
	TCut all_cuts = Ldelta + Rdelta + Rcoll + Lcoll + ePID;

	// Canvas for vertex comparison
	TCanvas * myCanv = new TCanvas;
		// H
	myCanv->Divide(2,1);
	myCanv->cd(1);
	inTree_3H_uptolast->Draw("L_z:R_z>>vert_3H_uptolast",all_cuts,"COLZ");
	vert_3H_uptolast->GetXaxis()->SetTitle("R_{z}");
	vert_3H_uptolast->GetYaxis()->SetTitle("L_{z}");
	vert_3H_uptolast->SetStats(kFALSE);
	vert_3H_uptolast->GetYaxis()->SetTitleOffset(1.8);
	myCanv->Update();
	myCanv->cd(2);
	inTree_3H_last->Draw("L_z:R_z>>vert_3H_last",all_cuts,"COLZ");
	vert_3H_last->GetXaxis()->SetTitle("R_{z}");
	vert_3H_last->GetYaxis()->SetTitle("L_{z}");
	vert_3H_last->SetStats(kFALSE);
	vert_3H_last->GetYaxis()->SetTitleOffset(1.8);
	myCanv->Update();
	myCanv->Print("/adaqfs/home/a-onl/tritium_work/Rey/report/plots/3H_vertCompare.pdf");
	myCanv->Clear();
		// He
	myCanv->Divide(2,1);
	myCanv->cd(1);
	inTree_3He_uptolast->Draw("L_z:R_z>>vert_3He_uptolast",all_cuts,"COLZ");
	vert_3He_uptolast->GetXaxis()->SetTitle("R_{z}");
	vert_3He_uptolast->GetYaxis()->SetTitle("L_{z}");
	vert_3He_uptolast->SetStats(kFALSE);
	vert_3He_uptolast->GetYaxis()->SetTitleOffset(1.8);
	myCanv->Update();
	myCanv->cd(2);
	inTree_3He_last->Draw("L_z:R_z>>vert_3He_last",all_cuts,"COLZ");
	vert_3He_last->GetXaxis()->SetTitle("R_{z}");
	vert_3He_last->GetYaxis()->SetTitle("L_{z}");
	vert_3He_last->SetStats(kFALSE);
	vert_3He_last->GetYaxis()->SetTitleOffset(1.8);
	myCanv->Update();
	myCanv->Print("/adaqfs/home/a-onl/tritium_work/Rey/report/plots/3He_vertCompare.pdf");
	myCanv->Clear();
	
	// Momentum plot
		// For H
	myCanv->Divide(2,1);
	myCanv->cd(1);
	inTree_3H_uptolast->Draw("L_mom:R_mom>>mom_3H_uptolast",all_cuts,"COLZ");
        mom_3H_uptolast->GetXaxis()->SetTitle("RHRS Mommentum [GeV/c]");
        mom_3H_uptolast->GetYaxis()->SetTitle("LHRS Mommentum [GeV/c]");
        mom_3H_uptolast->SetStats(kFALSE);
        mom_3H_uptolast->GetYaxis()->SetTitleOffset(1.8);
        myCanv->Update();
	myCanv->cd(2);
	inTree_3H_last->Draw("L_mom:R_mom>>mom_3H_last",all_cuts,"COLZ");
        mom_3H_last->GetXaxis()->SetTitle("RHRS Mommentum [GeV/c]");
        mom_3H_last->GetYaxis()->SetTitle("LHRS Mommentum [GeV/c]");
        mom_3H_last->SetStats(kFALSE);
        mom_3H_last->GetYaxis()->SetTitleOffset(1.8);
	myCanv->Update();
	//myCanv->Print("/adaqfs/home/a-onl/tritium_work/segarrae/3H_momCompare.pdf");
	myCanv->Print("/adaqfs/home/a-onl/tritium_work/Rey/report/plots/3H_momCompare.pdf");
	myCanv->Clear();
		// For He
	myCanv->Divide(2,1);
	myCanv->cd(1);
	inTree_3He_uptolast->Draw("L_mom:R_mom>>mom_3He_uptolast",all_cuts,"COLZ");
        mom_3He_uptolast->GetXaxis()->SetTitle("RHRS Mommentum [GeV/c]");
        mom_3He_uptolast->GetYaxis()->SetTitle("LHRS Mommentum [GeV/c]");
        mom_3He_uptolast->SetStats(kFALSE);
        mom_3He_uptolast->GetYaxis()->SetTitleOffset(1.8);
        myCanv->Update();
	myCanv->cd(2);
	inTree_3He_last->Draw("L_mom:R_mom>>mom_3He_last",all_cuts,"COLZ");
        mom_3He_last->GetXaxis()->SetTitle("RHRS Mommentum [GeV/c]");
        mom_3He_last->GetYaxis()->SetTitle("LHRS Mommentum [GeV/c]");
        mom_3He_last->SetStats(kFALSE);
        mom_3He_last->GetYaxis()->SetTitleOffset(1.8);
	myCanv->Update();
	myCanv->Print("/adaqfs/home/a-onl/tritium_work/Rey/report/plots/3He_momCompare.pdf");
	myCanv->Clear();


	delete myCanv;
	
	return 0;
}

