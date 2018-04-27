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

using namespace std;
const double mP = 0.938272;
const double me = 0.000511;

double CalcEbeam( double L_theta , double R_theta , double e_yptar, double e_xptar, double p_yptar, double p_xptar );

int main( int argc, char ** argv){
	if ( argc < 1 ){
		cerr << "Wrong number of arguments. Instead use\n"
		     << "\t./nudgeFact\n";
		return -1;
	}
	
	// Load file, tree, and input parameters
	TFile * inFile = new TFile("/adaqfs/home/a-onl/tritium_work/Rey/Event_Counter/skimmer/combiner_out/1H_mid.root");
	TTree * inTree = (TTree *) inFile->Get("sk");
	TVectorT <double> *kin = (TVectorT <double> *) inFile->Get("kinematics");
	double kinThetaL = (*kin)[0]; 
	double kinMomL	 = (*kin)[1];
	double kinThetaR = (*kin)[2]; 
	double kinMomR	 = (*kin)[3];
	double kinEbeam  = (*kin)[4];
		// Setup reading branches event by event
	double QSq, Nu;
	double L_EoverP, R_ytar, R_yptar, R_xptar, L_z, R_z;
	double L_ytar, L_yptar, L_xptar;
	double L_mom, L_theta, L_phi;
	double R_mom, R_theta, R_phi;
	double Pm, Em, Pmx, Pmy, Pmz;
	double thrq, L_delta, R_delta;
	inTree->SetBranchAddress("Q2",		&QSq		);
	inTree->SetBranchAddress("Nu",		&Nu		);
		// Branches used for cuts to clean event selection
	inTree->SetBranchAddress("L_EoverP",	&L_EoverP	);
	inTree->SetBranchAddress("R_ytar",	&R_ytar		);
	inTree->SetBranchAddress("R_yptar",	&R_yptar	);
	inTree->SetBranchAddress("R_xptar",	&R_xptar	);
	inTree->SetBranchAddress("L_ytar",	&L_ytar		);
	inTree->SetBranchAddress("L_yptar",	&L_yptar	);
	inTree->SetBranchAddress("L_xptar",	&L_xptar	);
	inTree->SetBranchAddress("L_z",		&L_z		);
	inTree->SetBranchAddress("R_z",		&R_z		);
		// Left arm
	inTree->SetBranchAddress("L_mom",	&L_mom		);
	//inTree->SetBranchAddress("L_theta",	&L_theta	);
	inTree->SetBranchAddress("L_phi",	&L_phi          );
        inTree->SetBranchAddress("R_z", 	&R_z		);
		// Right arm
	inTree->SetBranchAddress("R_mom",	&R_mom		);
	//inTree->SetBranchAddress("R_theta",	&R_theta	);
	inTree->SetBranchAddress("R_phi",	&R_phi		);
		// Missing kinematics
	inTree->SetBranchAddress("Pm",		&Pm		);
	inTree->SetBranchAddress("Em",		&Em		);
	inTree->SetBranchAddress("Pmx",		&Pmx		);
	inTree->SetBranchAddress("Pmy",		&Pmy		);
	inTree->SetBranchAddress("Pmz",		&Pmz		);
		// Other kinematics used for cuts
	inTree->SetBranchAddress("thrq",	&thrq		);
	inTree->SetBranchAddress("L_delta",	&L_delta	);
	inTree->SetBranchAddress("R_delta",	&R_delta	);


	TString P1title;
	P1title.Form("|E_{e,expected} - %f GeV|",floor(kinEbeam*100)/100.);
	TH1D * hP1_a = new TH1D("hP1_a",P1title,100000,-50,50);
	TH1D * hP1_b = new TH1D("hP1_b",P1title,100000,-50,50);
	TH1D * hP1_c = new TH1D("hP1_c",P1title,100000,-50,50);
	TFile * outFile = new TFile("movingE1.root","RECREATE");	


	TF1 * f1_a = new TF1("f1_a","gaus");
	TF1 * f1_b = new TF1("f1_b","gaus");
	TF1 * f1_c = new TF1("f1_c","gaus");
	
	// And now we want to repeat the above until we minimize the mean
	const int nEvents = inTree->GetEntries();
	//double nudgeL[9] = { -5./1000. , -4./1000. ,-3./1000. , -2./1000. , 0./1000. ,2./1000. , 3./1000. ,4./1000. , 5./1000. , };
	//		did the above, and saw that the absolute smallest nudge is around 0
	//double nudgeL[9] = { -0.5/1000. , -0.4/1000. ,-0.3/1000. , -0.2/1000. , 0./1000. ,0.2/1000. , 0.3/1000. ,0.4/1000. , 0.5/1000. , };
	//		did the above and saw that the absolute smallest nudge is around 0.5
	//double nudgeL[11] = { 0./1000. , 0.1/1000. ,0.2/1000. , 0.3/1000. , 0.4/1000. ,0.5/1000. , 0.6/1000. ,0.7/1000. , 0.8/1000. , 0.9/1000., 1./1000.};
	//		did the above and saw that the absolute smallest nudge is around 0.7
	//double nudgeL[11] = { 0.60/1000. , 0.62/1000. ,0.64/1000. , 0.66/1000. , 0.68/1000. ,0.7/1000. , 0.72/1000. ,0.74/1000. , 0.76/1000. , 0.78/1000., 0.8/1000.};
	// 		did above and around 0.68
	//double nudgeL[11] = { 0.655/1000. , 0.66/1000. , 0.665/1000. ,0.67/1000. , 0.675/1000. , 0.68/1000. ,0.685/1000. , 0.69/1000. ,0.695/1000. , 0.7/1000. , 0.705/1000.};
	//	did the above and new center is around 0.67 mrad, but we are reaching sensitivity now
	double nudgeL[9] = { 0.665/1000. , 0.667/1000. ,0.668/1000. , 0.669/1000. , 0.67/1000. ,0.671/1000. , 0.672/1000. ,0.673/1000. , 0.675/1000.};
	//	REACHED PRECISION
	//		nudge factors are: LEFT-0.00067 RIGHT-0.000445557
	double a,b,c;
	double fa,fb,fc;
	double rad = 1e5;
	double fin_nudgeL, fin_nudgeR;
	for ( int j = 0 ; j < 9; j++){
		//cout << j << "\n";
		//cout << nudgeL[j] << "\n";
		a = -50./1000.;
		b = 50./1000.;
		// Now for both a and b and nudgeL[j], we want to recreate the distribution, and 
		// find out what the new means are
		hP1_a->Reset();
		hP1_b->Reset();
		f1_a->SetParameters(0,0,0);
		f1_b->SetParameters(0,0,0);
		for (int i = 0 ; i < nEvents ; i++){
			inTree->GetEntry(i);		
			
			// Apply event cleaning cuts:
			if( L_EoverP < 0.5 ) continue;
			if( TMath::Abs( L_z - R_z - 0.0037 ) > 0.018 ) continue;
			if( TMath::Abs( L_delta ) > 0.045 ) continue;
			if( TMath::Abs( R_delta ) > 0.045 ) continue;
			if( TMath::Abs( L_z ) > 0.09 ) continue;
			//if( thrq* 180./M_PI > 40. ) continue;
			if( pow(((R_ytar+1.5*R_yptar)/0.08),2) + pow(((1.5*R_xptar)/0.08),2) > 1 ) continue;
			if( pow(((L_ytar+1.5*L_yptar)/0.08),2) + pow(((1.5*L_xptar)/0.08),2) > 1 ) continue;


			// Now calculate measurement momentum
			// versus assumed momentum to understand
			// nudge factors we'll need to apply
				// Unshifted quantities
			double E1_a = CalcEbeam( kinThetaL + nudgeL[j] , kinThetaR + a , L_yptar , L_xptar , R_yptar , R_xptar );
			double E1_b = CalcEbeam( kinThetaL + nudgeL[j] , kinThetaR + b , L_yptar , L_xptar , R_yptar , R_xptar );
			double del_p1_a = E1_a - kinEbeam;
			double del_p1_b = E1_b - kinEbeam;
			hP1_a->Fill( del_p1_a );
			hP1_b->Fill( del_p1_b );

		}
		// Now we can get the mean of hP1
		f1_a->SetParameter(0,hP1_a->GetMaximum());
		f1_a->SetParameter(1,hP1_a->GetMean());
		f1_a->SetParameter(2,hP1_a->GetRMS());
		f1_b->SetParameter(0,hP1_b->GetMaximum());
		f1_b->SetParameter(1,hP1_b->GetMean());
		f1_b->SetParameter(2,hP1_b->GetRMS());

		hP1_a->Fit("f1_a","QES");
		hP1_b->Fit("f1_b","QES");
		

		// And now we want to repeat the above until we minimize the mean
		fa = f1_a->GetParameter(1);
		fb = f1_b->GetParameter(1);
		int nCnt = 0;
		if( ( (fa < 0) && (fb > 0) ) || ( (fa > 0) && (fb < 0)   ) ){
			while( TMath::Abs(a-b) > 0.01/1000){
				hP1_c->Reset();
				c = (a+b)/2.;
				// Now reloop again to find the mean based on this nudge
				for (int i = 0 ; i < nEvents ; i++){
					inTree->GetEntry(i);		
					// Apply event cleaning cuts:
					if( L_EoverP < 0.5 ) continue;
					if( TMath::Abs( L_z - R_z - 0.0037 ) > 0.018 ) continue;
					if( TMath::Abs( L_delta ) > 0.045 ) continue;
					if( TMath::Abs( R_delta ) > 0.045 ) continue;
					if( TMath::Abs( L_z ) > 0.09 ) continue;
					//if( thrq* 180./M_PI > 40. ) continue;
					if( pow(((R_ytar+1.5*R_yptar)/0.08),2) + pow(((1.5*R_xptar)/0.08),2) > 1 ) continue;
					if( pow(((L_ytar+1.5*L_yptar)/0.08),2) + pow(((1.5*L_xptar)/0.08),2) > 1 ) continue;

					// Now calculate measurement momentum
					// versus assumed momentum to understand
					// nudge factors we'll need to apply
						// Unshifted quantities
					double E1_c = CalcEbeam( kinThetaL + nudgeL[j] , kinThetaR + c , L_yptar , L_xptar , R_yptar , R_xptar );
					double del_p1_c = E1_c - kinEbeam;
					hP1_c->Fill( del_p1_c );
				}
				// Now we can get the mean of hP1
				
				//if( nCnt == 0) cout << "Max: " << hP1_c->GetMaximum() << " Mean: " << hP1_c->GetMean() << " RMS: " << hP1_c->GetRMS() << "\n";
				nCnt++;
				f1_c->SetParameter(0,hP1_c->GetMaximum());
				f1_c->SetParameter(1,hP1_c->GetMean());
				f1_c->SetParameter(2,hP1_c->GetRMS());
				hP1_c->Fit("f1_c","QES");
				fc = f1_c->GetParameter(1);
				double errC = f1_c->GetParError(1);
				// Now find where the root is
				if( (fc < 0) && (fa > 0) ){
					b = c;
					fb = fc;
				}
				else if( (fc > 0) && (fb < 0) ){
					a = c;
					fa = fc;
				}
				else{
					cout << "WHAT\n";
					return 0;
				}
			}
		}
		else{
			cout << "Will not converge, choose better initial a,b nudge\n";
		}
		cout << "Nudge factors: " << nudgeL[j] << " " << c << "\n";

		double testRad = pow( pow(nudgeL[j],2) + pow(c,2) , 0.5);
		if( testRad < rad){
			rad = testRad;
			fin_nudgeL = nudgeL[j];
			fin_nudgeR = c;
		}
	}
	cout << fin_nudgeL << " " << fin_nudgeR << "\n";
	hP1_a->Write();
	hP1_b->Write();
	hP1_c->Write();
	outFile->Close();
	inFile->Close();

	return 0;
}
double CalcEbeam( double L_theta , double R_theta , double e_yptar, double e_xptar, double p_yptar, double p_xptar ){
	double e_scat = TMath::ACos( (TMath::Cos(L_theta) - e_yptar*TMath::Sin(L_theta)) / TMath::Sqrt(1. + e_yptar*e_yptar + e_xptar*e_xptar));
	double p_scat = TMath::ACos( (TMath::Cos(R_theta) + p_yptar*TMath::Sin(R_theta)) / TMath::Sqrt(1. + p_yptar*p_yptar + p_xptar*p_xptar));
	return mP * (  1./(tan(e_scat/2.)) * 1./(tan(p_scat)) - 1. );	
}
