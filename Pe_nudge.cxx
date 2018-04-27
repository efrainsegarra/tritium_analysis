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
	inTree->SetBranchAddress("L_theta",	&L_theta	);
	inTree->SetBranchAddress("L_phi",	&L_phi          );
        inTree->SetBranchAddress("R_z", 	&R_z		);
		// Right arm
	inTree->SetBranchAddress("R_mom",	&R_mom		);
	inTree->SetBranchAddress("R_theta",	&R_theta	);
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


	TH1D * hP3_a = new TH1D("hP3_a","p_{e',expected} - p_{e',measured}",20000,-5,5);	
	TH1D * hP3_b = new TH1D("hP3_b","p_{e',expected} - p_{e',measured}",20000,-5,5);	
	TH1D * hP3_c = new TH1D("hP3_c","p_{e',expected} - p_{e',measured}",20000,-5,5);	
	TFile * outFile = new TFile("movingE1.root","RECREATE");	

	TF1 * f1_a = new TF1("f1_a","gaus");
	TF1 * f1_b = new TF1("f1_b","gaus");
	TF1 * f1_c = new TF1("f1_c","gaus");
	
	const int nEvents = inTree->GetEntries();
	double a,b,c;
	double fa,fb,fc;
	double del_p3_a, del_p3_b, del_p3_c;
	a = -100./1000.;
	b =  100./1000.;
	// Now for both a and b and nudgeL[j], we want to recreate the distribution, and 
	// find out what the new means are
	hP3_a->Reset();
	hP3_b->Reset();
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

		// Now we don't need to calculate beam energy from angles
		// anymore, but we do need to calculate p3_expected based on
		// energy and angles to compare it to p3_measured
		
		double p3_exp   = (kinEbeam*mP) / ( mP + kinEbeam*(1-cos(L_theta))  );
		del_p3_a = p3_exp  - ( kinMomL + a )*( 1 + L_delta ) ;
		del_p3_b = p3_exp  - ( kinMomL + b )*( 1 + L_delta ) ;
		
		hP3_a->Fill( del_p3_a );
		hP3_b->Fill( del_p3_b );
		
	}
	
	// Now we can get the mean of hP1
	f1_a->SetParameter(0,hP3_a->GetMaximum());
	f1_a->SetParameter(1,hP3_a->GetMean());
	f1_a->SetParameter(2,hP3_a->GetRMS());
	f1_a->SetRange( hP3_a->GetMean() - 3*hP3_a->GetRMS()  , hP3_a->GetMean() + 0.01*hP3_a->GetRMS()  );
	f1_b->SetParameter(0,hP3_b->GetMaximum());
	f1_b->SetParameter(1,hP3_b->GetMean());
	f1_b->SetParameter(2,hP3_b->GetRMS());
	f1_b->SetRange( hP3_b->GetMean() - 3*hP3_b->GetRMS()  , hP3_b->GetMean() + 0.01*hP3_b->GetRMS()  );
	hP3_c->Fit("f1_c","QERS");
	hP3_a->Fit("f1_a","QERS");
	hP3_b->Fit("f1_b","QERS");
	

	// And now we want to repeat the above until we minimize the mean
	fa = f1_a->GetParameter(1);
	fb = f1_b->GetParameter(1);
	cout << "Before loop: " << fa << " " << fb << "\n";
	int nCnt = 0;
	if( ( (fa < 0) && (fb > 0) ) || ( (fa > 0) && (fb < 0)   ) ){
		while( TMath::Abs(a-b) > 0.00001/1000){
			hP3_c->Reset();
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
				double p3_exp_c   = (kinEbeam*mP) / ( mP + kinEbeam*(1-cos(L_theta))  );
				del_p3_c = p3_exp_c  - ( kinMomL + c )*( 1 + L_delta ) ;
				
				hP3_c->Fill( del_p3_c );
			}
			// Now we can get the mean of hP1
			
			//if( nCnt == 0) cout << "Max: " << hP1_c->GetMaximum() << " Mean: " << hP1_c->GetMean() << " RMS: " << hP1_c->GetRMS() << "\n";
			nCnt++;
			f1_c->SetParameter(0,hP3_c->GetMaximum());
			f1_c->SetParameter(1,hP3_c->GetMean());
			f1_c->SetParameter(2,hP3_c->GetRMS());
			f1_c->SetRange( hP3_c->GetMean() - 3*hP3_c->GetRMS()  , hP3_c->GetMean() + 0.01*hP3_c->GetRMS()  );
			hP3_c->Fit("f1_c","QERS");
			fc = f1_c->GetParameter(1);
			cout << "After loop: midpoint is " << fc <<  "\n";
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
			cout << "After sets: " << fa << " " << fb <<  "\n\n";
		}
	}
	else{
		cout << "Will not converge, choose better initial a,b nudge\n";
	}
	cout << "Final nudge factor for Pe correction: " << c << "\n";
	hP3_a->Write();
	hP3_b->Write();
	hP3_c->Write();
	outFile->Close();
	inFile->Close();

	return 0;
}
double CalcEbeam( double L_theta , double R_theta , double e_yptar, double e_xptar, double p_yptar, double p_xptar ){
	double e_scat = TMath::ACos( (TMath::Cos(L_theta) - e_yptar*TMath::Sin(L_theta)) / TMath::Sqrt(1. + e_yptar*e_yptar + e_xptar*e_xptar));
	double p_scat = TMath::ACos( (TMath::Cos(R_theta) + p_yptar*TMath::Sin(R_theta)) / TMath::Sqrt(1. + p_yptar*p_yptar + p_xptar*p_xptar));
	return mP * (  1./(tan(e_scat/2.)) * 1./(tan(p_scat)) - 1. );	
}
