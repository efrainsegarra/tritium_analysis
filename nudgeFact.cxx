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
#include "TLine.h"
#include "TF1.h"

using namespace std;
const double mP = 0.938272;
const double me = 0.000511;

int main( int argc, char ** argv){
	if ( argc < 1 ){
		cerr << "Wrong number of arguments. Instead use\n"
		     << "\t./nudgeFact\n";
		return -1;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////
	// Set up all the histograms that we will write out
	//	first just optimized kinematic distributions
	TH2D * Pp_Pe  = new TH2D("Pp_Pe",	"p_{e'} vs p_{p}",		200,3.4,3.6,150,1.4,1.55);
	TH2D * DelPhi = new TH2D("DelPhi",	"#phi_{e'} vs #phi_{p}",	200,-10,10,100,-5,5);
	TH2D * Tp_Te  = new TH2D("Tp_Te",	"#theta_{e'} vs #theta_{p}",	300,16.5,19.5,350,46.5,50);
	TH2D * Pp_Tp  = new TH2D("Pp_Tp",	"p_{p} vs #theta_{p}",  	350,46.5,50,150,1.4,1.55);
	TH2D * Pe_Te  = new TH2D("Pe_Te",	"p_{e'} vs #theta_{e'}", 	300,16.5,19.5,200,3.4,3.6 );
	//	now distributions to compare expected and measured
	TH1D * hP1 	= new TH1D("hP1",	"",					200,-0.1,0.1);
	TH1D * hP1_pre	= new TH1D("hP1_pre",	"",					200,-0.1,0.1);
	TH1D * hP3 	= new TH1D("hP3",	"p_{e',expected} - p_{e',measured}",	300,-0.01,0.02);
	TH1D * hP3_pre 	= new TH1D("hP3_pre",	"p_{e',expected} - p_{e',measured}",	300,-0.01,0.02);
	TH1D * hP4 	= new TH1D("hP4",	"p_{p,expected} - p_{p,mesaured}",	200,-0.01,0.01);
	TH1D * hP4_pre 	= new TH1D("hP4_pre",	"p_{p,expected} - p_{p,mesaured}",	200,-0.01,0.01);
	TH1D * hW	= new TH1D("hW ",	"W_{expected} - W_{measured}",		300,-0.02,0.01);
	TH1D * hW_pre	= new TH1D("hW_pre ",	"W_{expected} - W_{measured}",		300,-0.02,0.01);
	//	output histograms for pmiss
	TH1D * hPmx 	= new TH1D("hPmx",	"p_{miss,x}",	250,-0.01,0.015);
	TH1D * hPmx_pre = new TH1D("hPmx_pre",	"p_{miss,x}",	250,-0.01,0.015);
	TH1D * hPmy 	= new TH1D("hPmy",	"p_{miss,y}",	325,-0.025,0.04);
	TH1D * hPmy_pre = new TH1D("hPmy_pre",	"p_{miss,y}",	325,-0.025,0.04);
	TH1D * hPmz 	= new TH1D("hPmz",	"p_{miss,z}",	200,-0.01,0.01);
	TH1D * hPmz_pre = new TH1D("hPmz_pre",	"p_{miss,z}",	200,-0.01,0.01);
	TH1D * hEm  	= new TH1D("hEm",	"E_{miss}",	250,-0.010,0.015);
	TH1D * hEm_pre  = new TH1D("hEm_pre",	"E_{miss}",	250,-0.010,0.015);
	//	Making histograms nice
	DelPhi	->SetStats(kFALSE);
	DelPhi	->GetXaxis()->SetTitle("#phi_{e'} [Degrees]");	
	DelPhi	->GetYaxis()->SetTitle("#phi_{p} [Degrees]");	
	Pp_Pe	->SetStats(kFALSE);
	Pp_Pe	->GetXaxis()->SetTitle("p_{e'} [GeV]");	
	Pp_Pe	->GetYaxis()->SetTitle("p_{p} [GeV]");	
	Pp_Tp	->SetStats(kFALSE);
	Pp_Tp	->GetXaxis()->SetTitle("#theta_{p} [Degrees]");	
	Pp_Tp	->GetYaxis()->SetTitle("p_{p} [GeV]");	
	Pe_Te	->SetStats(kFALSE);
	Pe_Te	->GetXaxis()->SetTitle("#theta_{e'} [Degrees]");	
	Pe_Te	->GetYaxis()->SetTitle("p_{e'} [GeV]");	
	Tp_Te	->SetStats(kFALSE);
	Tp_Te	->GetXaxis()->SetTitle("#theta_{e'} [Degrees]");	
	Tp_Te	->GetYaxis()->SetTitle("#theta_{p} [Degrees]");	
	hP1	->SetStats(kFALSE);
	hP1	->GetXaxis()->SetTitle("Deviation [GeV]");
	hP1	->GetYaxis()->SetTitle("Counts");
	hP3	->SetStats(kFALSE);
	hP3	->GetXaxis()->SetTitle("Deviation [GeV]");
	hP3	->GetYaxis()->SetTitle("Counts");
	hP4	->SetStats(kFALSE);
	hP4	->GetXaxis()->SetTitle("Deviation [GeV]");
	hP4	->GetYaxis()->SetTitle("Counts");
	hW	->SetStats(kFALSE);
	hW	->GetXaxis()->SetTitle("Deviation [GeV]");
	hW	->GetYaxis()->SetTitle("Counts");
	hPmx	->SetStats(kFALSE);
	hPmx	->GetXaxis()->SetTitle("Measured [GeV]");
	hPmx	->GetYaxis()->SetTitle("Counts");
	hPmy	->SetStats(kFALSE);
	hPmy	->GetXaxis()->SetTitle("Measured [GeV]");
	hPmy	->GetYaxis()->SetTitle("Counts");
	hPmz	->SetStats(kFALSE);
	hPmz	->GetXaxis()->SetTitle("Measured [GeV]");
	hPmz	->GetYaxis()->SetTitle("Counts");
	hEm	->SetStats(kFALSE);
	hEm	->GetXaxis()->SetTitle("Measured [GeV]");
	hEm	->GetYaxis()->SetTitle("Counts");
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load file, tree, and input parameters
	TFile * inFile		= new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/1H_mid.root");
	TFile * inFile_pre	= new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/BACKUPS-beforeMemsetFix/mid-optimized/1H_mid.root");
	TTree * inTree		= (TTree *) inFile->Get("sk");
	TTree * inTree_pre 	= (TTree *) inFile_pre->Get("sk");
	//	Load the kinematics. These are what were shifted/optimized!!!
	TVectorT <double> *kin		= (TVectorT <double> *) inFile->Get("kinematics");
	TVectorT <double> *kin_pre	= (TVectorT <double> *) inFile_pre->Get("kinematics");
	double kinThetaL 	= (*kin)[0]; 
	double kinMomL	 	= (*kin)[1];
	double kinThetaR 	= (*kin)[2]; 
	double kinMomR	 	= (*kin)[3];
	double kinEbeam		= (*kin)[4];
	double kinThetaL_pre 	= (*kin_pre)[0]; 
	double kinMomL_pre	= (*kin_pre)[1];
	double kinThetaR_pre	= (*kin_pre)[2]; 
	double kinMomR_pre	= (*kin_pre)[3];
	double kinEbeam_pre	= (*kin_pre)[4];
	// Setup reading branches event by event
	double QSq, Nu;
	double L_EoverP, R_ytar, R_yptar, R_xptar, L_z, R_z;
	double L_ytar, L_yptar, L_xptar;
	double L_mom, L_theta, L_phi;
	double R_mom, R_theta, R_phi;
	double Pm, Em, Pmx, Pmy, Pmz;
	double thrq, L_delta, R_delta;
	// and same for the pre-optimized variables
	double QSq_pre,		Nu_pre;
	double L_EoverP_pre,	R_ytar_pre,	R_yptar_pre,	R_xptar_pre,	L_z_pre,	R_z_pre;
	double L_ytar_pre,	L_yptar_pre,	L_xptar_pre;
	double L_mom_pre,	L_theta_pre,	L_phi_pre;
	double R_mom_pre, 	R_theta_pre, 	R_phi_pre;
	double Pm_pre,		Em_pre,		Pmx_pre,	Pmy_pre,	Pmz_pre;
	double thrq_pre,	L_delta_pre,	R_delta_pre;
	// Set branch addresses
	inTree->SetBranchAddress("Q2",		&QSq		);
	inTree->SetBranchAddress("Nu",		&Nu		);
	inTree->SetBranchAddress("L_EoverP",	&L_EoverP	);
	inTree->SetBranchAddress("R_ytar",	&R_ytar		);
	inTree->SetBranchAddress("R_yptar",	&R_yptar	);
	inTree->SetBranchAddress("R_xptar",	&R_xptar	);
	inTree->SetBranchAddress("R_delta",	&R_delta	);
	inTree->SetBranchAddress("R_z",		&R_z		);
	inTree->SetBranchAddress("L_ytar",	&L_ytar		);
	inTree->SetBranchAddress("L_yptar",	&L_yptar	);
	inTree->SetBranchAddress("L_xptar",	&L_xptar	);
	inTree->SetBranchAddress("L_delta",	&L_delta	);
	inTree->SetBranchAddress("L_z",		&L_z		);
	inTree->SetBranchAddress("L_mom",	&L_mom		);
	inTree->SetBranchAddress("L_theta",	&L_theta	);
	inTree->SetBranchAddress("L_phi",	&L_phi          );
        inTree->SetBranchAddress("R_z", 	&R_z		);
	inTree->SetBranchAddress("R_mom",	&R_mom		);
	inTree->SetBranchAddress("R_theta",	&R_theta	);
	inTree->SetBranchAddress("R_phi",	&R_phi		);
	inTree->SetBranchAddress("Pm",		&Pm		);
	inTree->SetBranchAddress("Em",		&Em		);
	inTree->SetBranchAddress("Pmx",		&Pmx		);
	inTree->SetBranchAddress("Pmy",		&Pmy		);
	inTree->SetBranchAddress("Pmz",		&Pmz		);
	inTree->SetBranchAddress("thrq",	&thrq		);
	inTree_pre->SetBranchAddress("Q2",		&QSq_pre		);
	inTree_pre->SetBranchAddress("Nu",		&Nu_pre			);
	inTree_pre->SetBranchAddress("L_EoverP",	&L_EoverP_pre		);
	inTree_pre->SetBranchAddress("R_ytar",		&R_ytar_pre		);
	inTree_pre->SetBranchAddress("R_yptar",		&R_yptar_pre		);
	inTree_pre->SetBranchAddress("R_xptar",		&R_xptar_pre		);
	inTree_pre->SetBranchAddress("R_delta",		&R_delta_pre		);
	inTree_pre->SetBranchAddress("R_z",		&R_z_pre		);
	inTree_pre->SetBranchAddress("L_ytar",		&L_ytar_pre		);
	inTree_pre->SetBranchAddress("L_yptar",		&L_yptar_pre		);
	inTree_pre->SetBranchAddress("L_xptar",		&L_xptar_pre		);
	inTree_pre->SetBranchAddress("L_delta",		&L_delta_pre		);
	inTree_pre->SetBranchAddress("L_z",		&L_z_pre		);
	inTree_pre->SetBranchAddress("L_mom",		&L_mom_pre		);
	inTree_pre->SetBranchAddress("L_theta",		&L_theta_pre		);
	inTree_pre->SetBranchAddress("L_phi",		&L_phi_pre          	);
        inTree_pre->SetBranchAddress("R_z", 		&R_z_pre		);
	inTree_pre->SetBranchAddress("R_mom",		&R_mom_pre		);
	inTree_pre->SetBranchAddress("R_theta",		&R_theta_pre		);
	inTree_pre->SetBranchAddress("R_phi",		&R_phi_pre		);
	inTree_pre->SetBranchAddress("Pm",		&Pm_pre			);
	inTree_pre->SetBranchAddress("Em",		&Em_pre			);
	inTree_pre->SetBranchAddress("Pmx",		&Pmx_pre		);
	inTree_pre->SetBranchAddress("Pmy",		&Pmy_pre		);
	inTree_pre->SetBranchAddress("Pmz",		&Pmz_pre		);
	inTree_pre->SetBranchAddress("thrq",	&thrq_pre		);

	//////////////////////////////////////////////////////////////////////////////////
	// Ready output file
	TFile * outFile = new TFile("H_analysis.root","RECREATE");

	///////////////////////////////////////////////////////////////////////////////////
	// Ready variables used in looping over both files
	double del_p1,		del_p3,		del_p4,		del_W;
	double del_p1_pre,	del_p3_pre,	del_p4_pre,	del_W_pre;

	////////////////////////////////////////////////////////////////////////////////////
	// Loop over events:
	const int nEvents	= inTree->GetEntries();
	const int nEvents_pre	= inTree_pre->GetEntries();
	
	for (int i = 0 ; i < nEvents ; i++){
		// Get entries from BOTH trees
		inTree->GetEntry(i);		
		inTree_pre->GetEntry(i);
			
		// Apply event cleaning cuts for both trees so we don't mix any weird events
		if( L_EoverP < 0.5								) continue;
		if( TMath::Abs( L_z - R_z - 0.0037 ) > 0.018 					) continue;
		if( TMath::Abs( L_delta ) > 0.045 						) continue;
		if( TMath::Abs( R_delta ) > 0.045 						) continue;
		if( TMath::Abs( L_z ) > 0.09 							) continue;
		if( pow(((R_ytar+1.5*R_yptar)/0.08),2) + pow(((1.5*R_xptar)/0.08),2) > 1 	) continue;
		if( pow(((L_ytar+1.5*L_yptar)/0.08),2) + pow(((1.5*L_xptar)/0.08),2) > 1 	) continue;
		if( L_EoverP_pre < 0.5									) continue;
		if( TMath::Abs( L_z_pre - R_z_pre - 0.0037 ) > 0.018 					) continue;
		if( TMath::Abs( L_delta_pre ) > 0.045 							) continue;
		if( TMath::Abs( R_delta_pre ) > 0.045 							) continue;
		if( TMath::Abs( L_z_pre ) > 0.09 							) continue;
		if( pow(((R_ytar_pre+1.5*R_yptar_pre)/0.08),2) + pow(((1.5*R_xptar_pre)/0.08),2) > 1 	) continue;
		if( pow(((L_ytar_pre+1.5*L_yptar_pre)/0.08),2) + pow(((1.5*L_xptar_pre)/0.08),2) > 1 	) continue;
		
		// Fill the kinematical variables
		DelPhi->Fill( L_phi * 180./M_PI , 180 - R_phi * 180./M_PI ); 
		Pp_Pe ->Fill( L_mom , R_mom );
		Pp_Tp ->Fill( R_theta * 180./M_PI , R_mom );
		Pe_Te ->Fill( L_theta * 180./M_PI , L_mom );
		Tp_Te ->Fill( L_theta * 180./M_PI , R_theta * 180./M_PI );

		// Now calculate EXPECTED quantites based on angles of reconstructed scattering angles
		double E1   	= mP * (  1./(tan(L_theta/2.)) * 1./(tan(R_theta)) - 1. );	
		double p3_exp   = (kinEbeam*mP) / ( mP + kinEbeam*(1-cos(L_theta))  );
		double p4_exp   = pow( pow((kinEbeam + mP - L_mom),2) - pow(mP,2) , 0.5);
		del_p1 		= E1 - kinEbeam;
		del_p3 		= p3_exp - L_mom;
		del_p4 		= p4_exp - R_mom;		
		
		double W2 = 2*pow(me,2) + pow(mP,2) + 2*mP*E1 - 2*(E1*p3_exp - E1*p3_exp*cos(L_theta)) - 2*mP*p3_exp;
		double W = pow(W2,0.5); // this is exactly equal to my mP that I input, which is good cross-check...
		double W_measured = pow( -QSq + 2*Nu*mP + pow(mP,2) , 0.5);
		del_W = mP - W_measured;
		
		double gammaP = sqrt( 1 + pow(p4_exp/mP,2)  );
		double kinP = (gammaP -1)*mP;
		double Emiss = E1 - p3_exp - kinP;   // Emiss is just 0, again good cross-check
		
		// Fill histograms
		hP1->Fill( del_p1 );
		hP3->Fill( del_p3 );
		hP4->Fill( del_p4 );
		hW ->Fill( del_W  );
		hPmx->Fill( Pmx	  );
		hPmy->Fill( Pmy	  );
		hPmz->Fill( Pmz	  );
		hEm ->Fill( Em    );
	
		// Repeat EXPECTED quantites for pre-optimized values
		double E1_pre		= mP * (  1./(tan(L_theta_pre/2.)) * 1./(tan(R_theta_pre)) - 1. );	
		double p3_exp_pre	= (kinEbeam_pre*mP) / ( mP + kinEbeam_pre*(1-cos(L_theta_pre))  );
		double p4_exp_pre	= pow( pow((kinEbeam_pre + mP - L_mom_pre),2) - pow(mP,2) , 0.5);
		del_p1_pre		= E1_pre - kinEbeam_pre;
		del_p3_pre		= p3_exp_pre - L_mom_pre;
		del_p4_pre		= p4_exp_pre - R_mom_pre;		
		
		double W2_pre = 2*pow(me,2) + pow(mP,2) + 2*mP*E1_pre - 2*(E1_pre*p3_exp_pre - E1_pre*p3_exp_pre*cos(L_theta_pre)) - 2*mP*p3_exp_pre;
		double W_pre = pow(W2,0.5); // this is exactly equal to my mP that I input, which is good cross-check...
		double W_measured_pre = pow( -QSq_pre + 2*Nu_pre*mP + pow(mP,2) , 0.5);
		del_W_pre = mP - W_measured_pre;
		
		// Fill histograms
		hP1_pre	->Fill( del_p1_pre );
		hP3_pre	->Fill( del_p3_pre );
		hP4_pre	->Fill( del_p4_pre );
		hW_pre	->Fill( del_W_pre  );
		hPmx_pre->Fill( Pmx_pre	   );
		hPmy_pre->Fill( Pmy_pre	   );
		hPmz_pre->Fill( Pmz_pre	   );
		hEm_pre ->Fill( Em_pre     );
		
	}
	
	TString P1title;
	P1title.Form("E_{e,expected} - %f GeV",floor(kinEbeam*100)/100.);
	hP1->SetTitle(P1title);
	
	// Set up pad for canvas saving:
	TCanvas * myCanv = new TCanvas;
	TImage * img = TImage::Create();
	TString save;
	TLine *zero = NULL;

	// Saving all histograms
		// P1
	hP1->Draw();
	zero = new TLine(0,-1,0,hP1->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	hP1_pre->SetLineColor(kRed);
	hP1_pre->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delE1.png";
	img->WriteImage(save);
		// P3
	hP3->Draw();
	zero = new TLine(0,-1,0,hP3->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	hP3_pre->SetLineColor(kRed);
	hP3_pre->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delPe.png";
	img->WriteImage(save);
		// P4
	hP4->Draw();
	zero = new TLine(0,-1,0,hP4->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	hP4_pre->SetLineColor(kRed);
	hP4_pre->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delPp.png";
	img->WriteImage(save);
		// W
	hW ->Draw();
	zero = new TLine(0,-1,0,hW->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	hW_pre->SetLineColor(kRed);
	hW_pre->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delW.png";
	img->WriteImage(save);
		// Pmiss
	hPmx->Draw();
	hPmx_pre->SetLineColor(kRed);
	hPmx_pre->Draw("same");
	zero = new TLine(0,-1,0,hPmx->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pmx.png";
	img->WriteImage(save);
		// Pmy
	hPmy->Draw();
	hPmy_pre->SetLineColor(kRed);
	hPmy_pre->Draw("same");
	zero = new TLine(0,-1,0,hPmy->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pmy.png";
	img->WriteImage(save);
		// Pmz
	hPmz->Draw();
	hPmz_pre->SetLineColor(kRed);
	hPmz_pre->Draw("same");
	zero = new TLine(0,-1,0,hPmz->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pmz.png";
	img->WriteImage(save);
		// Em
	hEm->Draw();
	hEm_pre->SetLineColor(kRed);
	hEm_pre->Draw("same");
	zero = new TLine(0,-1,0,hEm->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Em.png";
	img->WriteImage(save);
		// Phi
	DelPhi->Draw("colz");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/DelPhi.png";
	img->WriteImage(save);
		// Pp_Pe
	Pp_Pe->Draw("colz");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pp_Pe.png";
	img->WriteImage(save);
		// Pp_Tp
	Pp_Tp->Draw("colz");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pp_Tp.png";
	img->WriteImage(save);
		// Pe_Te
	Pe_Te->Draw("colz");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pe_Te.png";
	img->WriteImage(save);
		// Tp_Te
	Tp_Te->Draw("colz");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Tp_Te.png";
	img->WriteImage(save);

	outFile->cd();
	hP1->Write();
	hP3->Write();
	hP4->Write();
	hW ->Write();
	hPmx->Write();
	hPmy->Write();
	hPmz->Write();
	hEm->Write();

	DelPhi->Write();
	Pp_Pe->Write();
	Pp_Tp->Write();
	Pe_Te->Write();
	Tp_Te->Write();
	
	outFile->Close();	
	
	inFile->Close();
	inFile_pre->Close();
	//delete myCanv;	
	return 0;
}
