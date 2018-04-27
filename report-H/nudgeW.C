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
#include "TMath.h"
#include "TLorentzVector.h"

using namespace std;
const double mP = 0.938272;
const double me = 0.000511;
const double mTarg = mP;
const double mD = 1.876;


int main( int argc, char ** argv){
	if ( argc < 1 ){
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\t./nudgeFact\n";
		return -1;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// Set up all the histograms that we will write out
	
	TH1D * hW_id  	= new TH1D("hW_id",	"W",	500,-0.02,0.03);
	TH1D * hW_ra	= new TH1D("hW_ra",	"W",	500,-0.02,0.03);
	TH1D * hW_exRa	= new TH1D("hW_exRa",	"W",	500,-0.02,0.03);

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load file, tree, and input parameters
	//TFile * inFile		= new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/1H_mid.root");
	TFile * inFile		= new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/1H_mid.root");
	TTree * thisTree	= (TTree *) inFile->Get("sk");
	
	// Setup reading branches event by event
	double 	sk_e_cer		,
		sk_e_prl1		,	sk_e_prl2		,
		sk_e_ntrk		,	sk_p_ntrk		,
		sk_e_ytar		,	sk_p_ytar		,
		sk_e_yptar		,	sk_p_yptar		,
		sk_gold_e_yptar		,	sk_gold_p_yptar		,
		sk_e_xptar		,	sk_p_xptar		,
		sk_gold_e_xptar		,	sk_gold_p_xptar		,
		sk_e_delta		,	sk_p_delta		,
		sk_gold_e_delta		,	sk_gold_p_delta		,
		sk_e_mom		,	sk_p_mom		,
		sk_gold_e_mom		,	sk_gold_p_mom		,
		sk_e_beta		,	sk_p_beta		,
		sk_gold_e_beta		,	sk_gold_p_beta		,
		sk_e_x			,	sk_p_x			,
		sk_e_y			,	sk_p_y			,
		sk_e_z			,	sk_p_z			,
		sk_e_px			,	sk_e_py			,
		sk_e_pz			,	sk_p_px			,
		sk_p_py			,	sk_p_pz			,
		sk_gold_e_px		,	sk_gold_p_px		,
		sk_gold_e_py		,	sk_gold_p_py		,
		sk_gold_e_pz		,	sk_gold_p_pz		,
		sk_t1			,	sk_t4			,
		sk_tcoinc		,
		sk_Pmiss		,	sk_Emiss		,
		sk_Pmiss_x		,	sk_Pmiss_y		,
		sk_Pmiss_z		,	sk_p_thWe		,
		sk_th_bq		,	sk_ph_bq		,
		sk_th_xq		,	sk_ph_xq		,
		sk_rast_Pmiss		,	sk_rast_Emiss		,
		sk_rast_Pmiss_x		,	sk_rast_Pmiss_y		,
		sk_rast_Pmiss_z		,	sk_p_rast_thWe		,
		sk_rast_th_bq		,	sk_rast_ph_bq		,
		sk_rast_th_xq		,	sk_rast_ph_xq		,
		sk_exRa_Pmiss		,	sk_exRa_Emiss		,
		sk_exRa_Pmiss_x		,	sk_exRa_Pmiss_y		,
		sk_exRa_Pmiss_z		,	sk_p_exRa_thWe		,
		sk_exRa_th_bq		,	sk_exRa_ph_bq		,
		sk_exRa_th_xq		,	sk_exRa_ph_xq		,
		sk_Q2			,	sk_W2			,
		sk_Nu			,	sk_ph_q			,
		sk_th_q			,	sk_xB			,
		sk_q3m			,	sk_q_x			,
		sk_q_y			,	sk_q_z			,
		sk_e_th			,
		sk_rast_Q2		,	sk_rast_W2		,
		sk_rast_Nu		,	sk_rast_ph_q		,
		sk_rast_th_q		,	sk_rast_xB		,
		sk_rast_q3m		,	sk_rast_q_x		,
		sk_rast_q_y		,	sk_rast_q_z		,		
		sk_e_rast_th		,
		sk_exRa_Q2		,	sk_exRa_W2		,
		sk_exRa_Nu		,	sk_exRa_ph_q		,
		sk_exRa_th_q		,	sk_exRa_xB		,
		sk_exRa_q3m		,	sk_exRa_q_x		,
		sk_exRa_q_y		,	sk_exRa_q_z		,
		sk_e_exRa_th		,
		sk_e_ext_deltaDp	,	sk_e_ext_deltaP		,
		sk_e_ext_deltaTh	,	sk_e_ext_delta		,
		sk_e_ext_mom		,	sk_e_ext_yptar		,
		sk_e_ext_xptar		,	sk_e_ext_px		,
		sk_e_ext_py		,	sk_e_ext_pz		,		
		sk_p_ext_deltaDp	,	sk_p_ext_deltaP		,
		sk_p_ext_deltaTh	,	sk_p_ext_delta		,
		sk_p_ext_mom		,	sk_p_ext_yptar		,
		sk_p_ext_xptar		,	sk_p_ext_px		,
		sk_p_ext_py		,	sk_p_ext_pz		,		
		sk_kin_eThetaC		,	sk_kin_pThetaC		,
		sk_kin_eMomC		,	sk_kin_pMomC		,
		sk_kin_Ebeam		,	sk_kin_Q		,
		sk_BCMcharge		,	sk_BCMcurr		,
		sk_BCMrenew 		;

	// Set branch addresses
	thisTree -> SetBranchAddress("L_prl1"		,&sk_e_prl1		);
	thisTree -> SetBranchAddress("L_prl2"		,&sk_e_prl2		);

	thisTree -> SetBranchAddress("L_cer"		,&sk_e_cer		);

	thisTree -> SetBranchAddress("L_ntrk"		,&sk_e_ntrk		);
	thisTree -> SetBranchAddress("R_ntrk"		,&sk_p_ntrk		);

	thisTree -> SetBranchAddress("L_ytar"		,&sk_e_ytar		);
	thisTree -> SetBranchAddress("R_ytar"		,&sk_p_ytar		);

	thisTree -> SetBranchAddress("L_yptar"		,&sk_e_yptar		);
	thisTree -> SetBranchAddress("R_yptar"		,&sk_p_yptar		);
	thisTree -> SetBranchAddress("L_gold_yptar"	,&sk_gold_e_yptar	);
	thisTree -> SetBranchAddress("R_gold_yptar"	,&sk_gold_p_yptar	);

	thisTree -> SetBranchAddress("L_xptar"		,&sk_e_xptar		);
	thisTree -> SetBranchAddress("R_xptar"		,&sk_p_xptar		);
	thisTree -> SetBranchAddress("L_gold_xptar"	,&sk_gold_e_xptar	);
	thisTree -> SetBranchAddress("R_gold_xptar"	,&sk_gold_p_xptar	);
	
	thisTree -> SetBranchAddress("L_dp"		,&sk_e_delta		);
	thisTree -> SetBranchAddress("R_dp"		,&sk_p_delta		);
	thisTree -> SetBranchAddress("L_gold_dp"	,&sk_gold_e_delta	);
	thisTree -> SetBranchAddress("R_gold_dp"	,&sk_gold_p_delta	);
	
	thisTree -> SetBranchAddress("L_mom"		,&sk_e_mom		);
	thisTree -> SetBranchAddress("R_mom"		,&sk_p_mom		);
	thisTree -> SetBranchAddress("L_gold_mom"	,&sk_gold_e_mom		);
	thisTree -> SetBranchAddress("R_gold_mom"	,&sk_gold_p_mom		);

	thisTree -> SetBranchAddress("L_beta"		,&sk_e_beta		);
	thisTree -> SetBranchAddress("R_beta"		,&sk_p_beta		);
	thisTree -> SetBranchAddress("L_gold_beta"	,&sk_gold_e_beta	);
	thisTree -> SetBranchAddress("R_gold_beta"	,&sk_gold_p_beta	);

	thisTree -> SetBranchAddress("L_vx"		,&sk_e_x		);
	thisTree -> SetBranchAddress("L_vy"		,&sk_e_y		);
	thisTree -> SetBranchAddress("L_vz"		,&sk_e_z		);
	thisTree -> SetBranchAddress("R_vx"		,&sk_p_x		);
	thisTree -> SetBranchAddress("R_vy"		,&sk_p_y		);
	thisTree -> SetBranchAddress("R_vz"		,&sk_p_z		);
	
	thisTree -> SetBranchAddress("L_px"		,&sk_e_px		);
	thisTree -> SetBranchAddress("L_py"		,&sk_e_py		);
	thisTree -> SetBranchAddress("L_pz"		,&sk_e_pz		);
	thisTree -> SetBranchAddress("R_px"		,&sk_p_px		);
	thisTree -> SetBranchAddress("R_py"		,&sk_p_py		);
	thisTree -> SetBranchAddress("R_pz"		,&sk_p_pz		);
	
	thisTree -> SetBranchAddress("L_gold_px"	,&sk_gold_e_px		);
	thisTree -> SetBranchAddress("L_gold_py"	,&sk_gold_e_py		);
	thisTree -> SetBranchAddress("L_gold_pz"	,&sk_gold_e_pz		);
	thisTree -> SetBranchAddress("R_gold_px"	,&sk_gold_p_px		);
	thisTree -> SetBranchAddress("R_gold_py"	,&sk_gold_p_py		);
	thisTree -> SetBranchAddress("R_gold_pz"	,&sk_gold_p_pz		);

	thisTree -> SetBranchAddress("t1"		,&sk_t1			);
	thisTree -> SetBranchAddress("t4"		,&sk_t4			);
	thisTree -> SetBranchAddress("tcoinc"		,&sk_tcoinc		);
	
	thisTree -> SetBranchAddress("Pmiss"		,&sk_Pmiss		);
	thisTree -> SetBranchAddress("Emiss"		,&sk_Emiss		);
	thisTree -> SetBranchAddress("Pmiss_x"		,&sk_Pmiss_x		);
	thisTree -> SetBranchAddress("Pmiss_y"		,&sk_Pmiss_y		);
	thisTree -> SetBranchAddress("Pmiss_z"		,&sk_Pmiss_z		);
	thisTree -> SetBranchAddress("R_thetaWe"	,&sk_p_thWe		);
	thisTree -> SetBranchAddress("ph_rq"		,&sk_ph_bq		);
	thisTree -> SetBranchAddress("th_rq"		,&sk_th_bq		);
	thisTree -> SetBranchAddress("ph_xq"		,&sk_ph_xq		);
	thisTree -> SetBranchAddress("th_xq"		,&sk_th_xq		);

	thisTree -> SetBranchAddress("rast_Pmiss"	,&sk_rast_Pmiss		);
	thisTree -> SetBranchAddress("rast_Emiss"	,&sk_rast_Emiss		);
	thisTree -> SetBranchAddress("rast_Pmiss_x"	,&sk_rast_Pmiss_x	);
	thisTree -> SetBranchAddress("rast_Pmiss_y"	,&sk_rast_Pmiss_y	);
	thisTree -> SetBranchAddress("rast_Pmiss_z"	,&sk_rast_Pmiss_z	);
	thisTree -> SetBranchAddress("rast_R_thetaWe"	,&sk_p_rast_thWe	);
	thisTree -> SetBranchAddress("rast_ph_rq"	,&sk_rast_ph_bq		);
	thisTree -> SetBranchAddress("rast_th_rq"	,&sk_rast_th_bq		);
	thisTree -> SetBranchAddress("rast_ph_xq"	,&sk_rast_ph_xq		);
	thisTree -> SetBranchAddress("rast_th_xq"	,&sk_rast_th_xq		);

	thisTree -> SetBranchAddress("exRa_Pmiss"	,&sk_exRa_Pmiss		);
	thisTree -> SetBranchAddress("exRa_Emiss"	,&sk_exRa_Emiss		);
	thisTree -> SetBranchAddress("exRa_Pmiss_x"	,&sk_exRa_Pmiss_x	);
	thisTree -> SetBranchAddress("exRa_Pmiss_y"	,&sk_exRa_Pmiss_y	);
	thisTree -> SetBranchAddress("exRa_Pmiss_z"	,&sk_exRa_Pmiss_z	);
	thisTree -> SetBranchAddress("exRa_R_thetaWe"	,&sk_p_exRa_thWe	);
	thisTree -> SetBranchAddress("exRa_ph_rq"	,&sk_exRa_ph_bq		);
	thisTree -> SetBranchAddress("exRa_th_rq"	,&sk_exRa_th_bq		);
	thisTree -> SetBranchAddress("exRa_ph_xq"	,&sk_exRa_ph_xq		);
	thisTree -> SetBranchAddress("exRa_th_xq"	,&sk_exRa_th_xq		);
	
	thisTree -> SetBranchAddress("Q2"		,&sk_Q2			);
	thisTree -> SetBranchAddress("W2"		,&sk_W2			);
	thisTree -> SetBranchAddress("Nu"		,&sk_Nu			);
	thisTree -> SetBranchAddress("ph_q"		,&sk_ph_q		);
	thisTree -> SetBranchAddress("th_q"		,&sk_th_q		);
	thisTree -> SetBranchAddress("xB"		,&sk_xB			);
	thisTree -> SetBranchAddress("q3m"		,&sk_q3m		);
	thisTree -> SetBranchAddress("q_x"		,&sk_q_x		);
	thisTree -> SetBranchAddress("q_y"		,&sk_q_y		);
	thisTree -> SetBranchAddress("q_z"		,&sk_q_z		);
	thisTree -> SetBranchAddress("L_theta"		,&sk_e_th		);

	thisTree -> SetBranchAddress("rast_Q2"		,&sk_rast_Q2		);
	thisTree -> SetBranchAddress("rast_W2"		,&sk_rast_W2		);
	thisTree -> SetBranchAddress("rast_Nu"		,&sk_rast_Nu		);
	thisTree -> SetBranchAddress("rast_ph_q"	,&sk_rast_ph_q		);
	thisTree -> SetBranchAddress("rast_th_q"	,&sk_rast_th_q		);
	thisTree -> SetBranchAddress("rast_xB"		,&sk_rast_xB		);
	thisTree -> SetBranchAddress("rast_q3m"		,&sk_rast_q3m		);
	thisTree -> SetBranchAddress("rast_q_x"		,&sk_rast_q_x		);
	thisTree -> SetBranchAddress("rast_q_y"		,&sk_rast_q_y		);
	thisTree -> SetBranchAddress("rast_q_z"		,&sk_rast_q_z		);
	thisTree -> SetBranchAddress("rast_L_theta"	,&sk_e_rast_th		);

	thisTree -> SetBranchAddress("exRa_Q2"		,&sk_exRa_Q2		);
	thisTree -> SetBranchAddress("exRa_W2"		,&sk_exRa_W2		);
	thisTree -> SetBranchAddress("exRa_Nu"		,&sk_exRa_Nu		);
	thisTree -> SetBranchAddress("exRa_ph_q"	,&sk_exRa_ph_q		);
	thisTree -> SetBranchAddress("exRa_th_q"	,&sk_exRa_th_q		);
	thisTree -> SetBranchAddress("exRa_xB"		,&sk_exRa_xB		);
	thisTree -> SetBranchAddress("exRa_q3m"		,&sk_exRa_q3m		);
	thisTree -> SetBranchAddress("exRa_q_x"		,&sk_exRa_q_x		);
	thisTree -> SetBranchAddress("exRa_q_y"		,&sk_exRa_q_y		);
	thisTree -> SetBranchAddress("exRa_q_z"		,&sk_exRa_q_z		);
	thisTree -> SetBranchAddress("exRa_L_theta"	,&sk_e_exRa_th		);

	thisTree -> SetBranchAddress("L_ext_delta_dp"	,&sk_e_ext_deltaDp	);
	thisTree -> SetBranchAddress("L_ext_delta_p"	,&sk_e_ext_deltaP	);
	thisTree -> SetBranchAddress("L_ext_delta_yptar",&sk_e_ext_deltaTh	);
	thisTree -> SetBranchAddress("L_ext_dp"		,&sk_e_ext_delta	);
	thisTree -> SetBranchAddress("L_ext_mom"	,&sk_e_ext_mom		);
	thisTree -> SetBranchAddress("L_ext_yptar"	,&sk_e_ext_yptar	);
	thisTree -> SetBranchAddress("L_ext_xptar"	,&sk_e_ext_xptar	);
	thisTree -> SetBranchAddress("L_ext_px"		,&sk_e_ext_px		);
	thisTree -> SetBranchAddress("L_ext_py"		,&sk_e_ext_py		);
	thisTree -> SetBranchAddress("L_ext_pz"		,&sk_e_ext_pz		);	

	thisTree -> SetBranchAddress("R_ext_delta_dp"	,&sk_p_ext_deltaDp	);
	thisTree -> SetBranchAddress("R_ext_delta_p"	,&sk_p_ext_deltaP	);
	thisTree -> SetBranchAddress("R_ext_delta_yptar",&sk_p_ext_deltaTh	);
	thisTree -> SetBranchAddress("R_ext_dp"		,&sk_p_ext_delta	);
	thisTree -> SetBranchAddress("R_ext_mom"	,&sk_p_ext_mom		);
	thisTree -> SetBranchAddress("R_ext_yptar"	,&sk_p_ext_yptar	);
	thisTree -> SetBranchAddress("R_ext_xptar"	,&sk_p_ext_xptar	);
	thisTree -> SetBranchAddress("R_ext_px"		,&sk_p_ext_px		);
	thisTree -> SetBranchAddress("R_ext_py"		,&sk_p_ext_py		);
	thisTree -> SetBranchAddress("R_ext_pz"		,&sk_p_ext_pz		);	

	thisTree -> SetBranchAddress("Kin_L_thetaC"	,&sk_kin_eThetaC	);
	thisTree -> SetBranchAddress("Kin_R_thetaC"	,&sk_kin_pThetaC	);
	thisTree -> SetBranchAddress("Kin_L_momC"	,&sk_kin_eMomC		);
	thisTree -> SetBranchAddress("Kin_R_momC"	,&sk_kin_pMomC		);
	thisTree -> SetBranchAddress("Kin_Ebeam"	,&sk_kin_Ebeam		);
	thisTree -> SetBranchAddress("Kin_Q"		,&sk_kin_Q		);
	
	thisTree -> SetBranchAddress("BCM_curr"		,&sk_BCMcurr		);
	thisTree -> SetBranchAddress("BCM_charge"	,&sk_BCMcharge		);
	thisTree -> SetBranchAddress("BCM_isrenew"	,&sk_BCMrenew		);

	//////////////////////////////////////////////////////////////////////////////////
	// Ready output file
	TFile * outFile = new TFile("H_analysis.root","RECREATE");

	///////////////////////////////////////////////////////////////////////////////////
	// Ready variables used in looping over both files
	double del_p1,		del_p3,		del_p4,		del_W;
	double del_p1_un,	del_p3_un,	del_p4_un,	del_W_un;

	////////////////////////////////////////////////////////////////////////////////////
	// Loop over events:
	const int nEvents	= thisTree->GetEntries();

	for (int i = 0 ; i < nEvents ; i++){
		// Get entries from BOTH trees
		thisTree->GetEntry(i);		

		// Apply event cleaning cuts for both trees so we don't mix any weird events
		if( pow(((sk_p_ytar+1.5*sk_p_ext_yptar)/0.08),2) + pow(((1.5*sk_p_ext_xptar)/0.08),2) <= 1 ){
		if( pow(((sk_e_ytar+1.5*sk_e_ext_yptar)/0.08),2) + pow(((1.5*sk_e_ext_xptar)/0.08),2) <= 1 ){
		if ( abs(sk_e_ext_delta)<0.045 &&
		     abs(sk_p_ext_delta)<0.045 ){
		if ( ((sk_e_prl1 + sk_e_prl2)/1000.)/(sk_e_ext_mom) > 0.8 && 
		     ((sk_e_prl1 + sk_e_prl2)/1000.)/(sk_e_ext_mom) < 1.3 ){
		if( sk_e_z > -0.1 && sk_e_z < 0.11 && 
		    sk_p_z > -0.1 && sk_p_z < 0.11 &&
		    TMath::Abs( sk_e_z - sk_p_z - 0.005) < 0.02 ){
		if( sk_tcoinc > 10){
			double sk_p_exRa_th = sk_p_exRa_thWe - sk_e_exRa_th;
			double sk_p_rast_th = sk_p_rast_thWe - sk_e_rast_th;
			double sk_p_th = sk_p_thWe - sk_e_th;
			double sk_e_phi = atan2( sk_e_ext_py , sk_e_ext_px  ) * 180./M_PI;
			double sk_p_phi = atan2( sk_p_ext_py , sk_p_ext_px  ) * 180./M_PI;
			if( sk_p_phi < 0.) sk_p_phi += 360.;
			
		
			double e_central1 = 17.8018 * M_PI/180.;	
			double p_central1 = -48.82  * M_PI/180.;	
			double calc_e_th1 = acos((cos(e_central1) - sk_e_yptar*sin(e_central1)) / sqrt( 1 + pow(sk_e_yptar,2) + pow(sk_e_xptar,2)  )); 
			double calc_p_th1 = acos((cos(p_central1) - sk_p_yptar*sin(p_central1)) / sqrt( 1 + pow(sk_p_yptar,2) + pow(sk_p_xptar,2)  )); 
			heThdiff_exRa->Fill( 180./M_PI*(calc_e_th1 - sk_e_exRa_th) );
			heThdiff_id  ->Fill( 180./M_PI*(calc_e_th1- sk_e_th)       );
			hpThdiff_exRa->Fill( 180./M_PI*(calc_p_th1 - sk_p_exRa_th) );
			hpThdiff_id  ->Fill( 180./M_PI*(calc_p_th1- sk_p_th)       );
			//cout << "p: "  << 180./M_PI*(calc_p_th1 - sk_p_exRa_th) << " " << 180./M_PI*(calc_p_th1- sk_p_th) << "\n";
			
			// So we can just use ideal beam to compare 'us' vs extended target / rastered
			// corrections on beam energy reconstruction
			
			// And now let's look at how much of an effect all these calibrations are..
			// by looking at Emiss/Pmiss first
			hPmx_id->Fill	( sk_Pmiss_x	);
			hPmy_id->Fill	( sk_Pmiss_y	);
			hPmz_id->Fill	( sk_Pmiss_z	);
			hEm_id ->Fill	( sk_Emiss	);
			hPmx_ra->Fill	( sk_rast_Pmiss_x	);
			hPmy_ra->Fill	( sk_rast_Pmiss_y	);
			hPmz_ra->Fill	( sk_rast_Pmiss_z	);
			hEm_ra ->Fill	( sk_rast_Emiss		);
			hPmx_exRa->Fill	( sk_exRa_Pmiss_x	);
			hPmy_exRa->Fill	( sk_exRa_Pmiss_y	);
			hPmz_exRa->Fill	( sk_exRa_Pmiss_z	);
			hEm_exRa ->Fill	( sk_exRa_Emiss		);
			hW_id 	  ->Fill( mP - sqrt(sk_W2)  );
			hW_ra 	  ->Fill( mP - sqrt(sk_rast_W2)  );
			hW_exRa   ->Fill( mP - sqrt(sk_exRa_W2)  );

			// It looks like Emiss is all messed up, so let's calculate it ourselves
			// and compare it
			double Em_calc_exRa = sk_exRa_Nu - (sqrt(pow(mP,2) + pow(sk_p_ext_mom,2)) - mP);
			double Em_calc_id   = sk_Nu - (sqrt(pow(mP,2) + pow(sk_p_mom,2)) - mP);
			hEm_calc_exRa -> Fill( Em_calc_exRa );
			hEm_calc_id   -> Fill( Em_calc_id  );
			
			// and then do a W cut on WHEREVER the mass of proton shows up
			// to see less of a tail, but of course the spectrum are all still shifted

			// Now we can reconstruct beam energy based on angle and then based on momenta
			double E1_an_exRa = mP * (  1./(tan(sk_e_exRa_th/2.)) * 1./(tan(sk_p_exRa_th)) - 1. );
			double E1_an_ra   = mP * (  1./(tan(sk_e_rast_th/2.)) * 1./(tan(sk_p_rast_th)) - 1. );
			double E1_an_id   = mP * (  1./(tan(sk_e_th/2.     )) * 1./(tan(sk_p_th     )) - 1. );
			hE1_an_exRa -> Fill( E1_an_exRa - sk_kin_Ebeam);
			hE1_an_ra   -> Fill( E1_an_ra   - sk_kin_Ebeam);
			hE1_an_id   -> Fill( E1_an_id   - sk_kin_Ebeam);
						
			double E1_p_exRa = sk_e_ext_mom + sqrt( pow(sk_p_ext_mom,2) + pow(mP,2) ) - mP;
			double E1_p_id   = sk_e_mom     + sqrt( pow(sk_p_mom    ,2) + pow(mP,2) ) - mP;
			hE1_p_exRa -> Fill( E1_p_exRa -sk_kin_Ebeam );
			hE1_p_id   -> Fill( E1_p_id   -sk_kin_Ebeam );
			
			// Let's just sudo add in the shifts I calculated last time, which
			// should be ROUGHLY okay, and see where this falls
			// -- of course I can only compare with ideal beam case, since
			// the equation is only equivalent for ideal beam case
			double e_central2 = 17.8018 * M_PI/180. + 0.00067;	
			double p_central2 = -48.82  * M_PI/180. - 0.000445557;
			double calc_e_th2 = acos((cos(e_central2) - sk_e_yptar*sin(e_central2)) / sqrt( 1 + pow(sk_e_yptar,2) + pow(sk_e_xptar,2)  )); 
			double calc_p_th2 = acos((cos(p_central2) - sk_p_yptar*sin(p_central2)) / sqrt( 1 + pow(sk_p_yptar,2) + pow(sk_p_xptar,2)  )); 
			double E1_an_id_nudge=mP*(  1./(tan(calc_e_th2/2.     )) * 1./(tan(calc_p_th2    )) - 1. );
			hE1_an_id_nudge-> Fill( E1_an_id_nudge   - sk_kin_Ebeam);
			
			// So with these 'fixes' to central angles, we can look at
			// the expected proton/electron momenta
			
			double p3_exp   = (E1_an_id_nudge*mP) / ( mP + E1_an_id_nudge*(1-cos(calc_e_th2))  );
			double p4_exp   = pow( pow((E1_an_id_nudge + mP - p3_exp),2) - pow(mP,2) , 0.5);
			hP3->Fill( p3_exp - sk_e_ext_mom );
			hP4->Fill( p4_exp - sk_p_ext_mom );

			// now with fixed angles, the point is we want to see how Emiss might change 
			// and pmiss might change but I need to REANALYZE those files to get the 
			// corrected angles in there...
		}
		}
		}
		}
		}
		}
	}
	
	// Set up pad for canvas saving:
	TCanvas * myCanv = new TCanvas;
	TImage * img = TImage::Create();
	TString save;
	TLine *zero = NULL;
	
	// Saving all histograms
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
	
	heThdiff_id->Scale(1./50);
	heThdiff_id->SetStats(kFALSE);
	heThdiff_id->GetXaxis()->SetTitle("Deviation [degrees]");
	heThdiff_id->Draw();
	heThdiff_exRa->SetLineColor(kRed);
	heThdiff_exRa->Draw("same");
	zero = new TLine(0,-1,0,heThdiff_id->GetMaximum()*1.05);
        zero->SetLineStyle(9);
        zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/thetaDiffe.png";
        img->WriteImage(save);


	hpThdiff_id->Scale(1./5);
	hpThdiff_id->SetStats(kFALSE);
	hpThdiff_id->GetXaxis()->SetTitle("Deviation [degrees]");
	hpThdiff_id->Draw();
	hpThdiff_exRa->SetLineColor(kRed);
	hpThdiff_exRa->Draw("same");
	zero = new TLine(0,-1,0,heThdiff_id->GetMaximum()*1.05);
        zero->SetLineStyle(9);
        zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/thetaDiffp.png";
        img->WriteImage(save);
	
	hW_id->Draw();
	hW_ra->SetLineColor(kBlack);
	hW_ra->Draw("same");
	hW_exRa->SetLineColor(kRed);
	hW_exRa->Draw("same");
	zero = new TLine(0,-1,0,hW_exRa->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delW.png";
	img->WriteImage(save);


	hPmx_id->Draw();
	hPmx_ra->SetLineColor(kBlack);
	hPmx_ra->Draw("same");
	hPmx_exRa->SetLineColor(kRed);
	hPmx_exRa->Draw("same");
	zero = new TLine(0,-1,0,hPmx_exRa->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pmx.png";
	img->WriteImage(save);

	hPmy_id->Draw();
	hPmy_ra->SetLineColor(kBlack);
	hPmy_ra->Draw("same");
	hPmy_exRa->SetLineColor(kRed);
	hPmy_exRa->Draw("same");
	zero = new TLine(0,-1,0,hPmy_exRa->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pmy.png";
	img->WriteImage(save);

	hPmz_id->Draw();
	hPmz_ra->SetLineColor(kBlack);
	hPmz_ra->Draw("same");
	hPmz_exRa->SetLineColor(kRed);
	hPmz_exRa->Draw("same");
	zero = new TLine(0,-1,0,hPmz_exRa->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Pmz.png";
	img->WriteImage(save);
	
	hEm_calc_exRa->SetLineColor(kRed);
	hEm_calc_exRa->Draw();
	hEm_exRa->Draw("same");
	zero = new TLine(0,-1,0,hEm_calc_exRa->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/Em.png";
	img->WriteImage(save);

	hE1_p_exRa->SetLineColor(kRed);
	hE1_p_exRa->Draw();
	hE1_an_exRa->SetLineColor(kBlue);
	hE1_an_exRa->Draw("same");
	hE1_an_id_nudge->SetLineColor(kBlack);
	hE1_an_id_nudge->Draw("same");
	zero = new TLine(0,-1,0,hE1_p_exRa->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delEnergy.png";
	img->WriteImage(save);
	
	// P3
	hP3->Draw();
	zero = new TLine(0,-1,0,hP3->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delPe.png";
	img->WriteImage(save);
	// P4
	hP4->Draw();
	zero = new TLine(0,-1,0,hP4->GetMaximum()*1.05);
	zero->SetLineStyle(9);
	zero->Draw("same");
	myCanv->Update();
	img->FromPad( myCanv );
	save = "/adaqfs/home/a-onl/tritium_work/segarrae/report-H/delPp.png";
	img->WriteImage(save);
	


	outFile->cd();
	
	heThdiff_exRa->Write();
	heThdiff_id  ->Write();
	hpThdiff_exRa->Write();
	hpThdiff_id  ->Write();

	hPmx_id-> Write();
	hPmy_id-> Write();
	hPmz_id-> Write();
	hEm_id -> Write();
	hPmx_ra-> Write();
	hPmy_ra-> Write();
	hPmz_ra-> Write();
	hEm_ra -> Write();
	hPmx_exRa-> Write();
	hPmy_exRa-> Write();
	hPmz_exRa-> Write();
	hEm_exRa -> Write();
	hEm_calc_exRa ->Write();
	hEm_calc_id   ->Write();
	hW_id ->Write();
	hW_ra ->Write();
	hW_exRa ->Write();

	hE1_an_exRa -> Write();
	hE1_an_ra   -> Write();
	hE1_an_id   -> Write();
	hE1_an_id_nudge-> Write();       
	hE1_p_exRa -> Write();
	hE1_p_id   -> Write();

	hP3->Write();
	hP4->Write();

	DelPhi->Write();
	Pp_Pe->Write();
	Pp_Tp->Write();
	Pe_Te->Write();
	Tp_Te->Write();

	outFile->Close();	

	inFile->Close();
	//delete myCanv;	
	return 0;
}
