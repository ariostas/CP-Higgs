// Include ROOT and C++ libraries
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>

#include <TLorentzVector.h>
#include <TVector3.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

// Declare functions
void histogram(TH1D*, TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TH1D*, TCanvas*, const char*, const char*, const char*);
void histogram(TH1D*, TCanvas*, const char*, const char*, const char*);
void saveResults();
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );
void analyze(TString, Int_t);

// Initialize histograms
TH1D *hThetaS1 = new TH1D("hThetaS1", "hThetaS1", 25, -3.14, 3.14);
TH1D *hThetaS2 = new TH1D("hThetaS2", "hThetaS2", 25, -3.14, 3.14);
TH1D *hThetaB = new TH1D("hThetaB", "hThetaB", 25, -3.14, 3.14);

using namespace std;

/*
 * MAIN FUNCTION
 */

 void cphiggs(){

 	analyze("CP0",1);
 	analyze("CP1p57",2);
 	analyze("ee_dy_tauhtauh_MG5Pythia8",0);
 	analyze("ee_ww_MG5Pythia8",0);
 	analyze("ee_zz_MG5Pythia8",0);
 	analyze("ee_zzee_MG5Pythia8",0);

 	saveResults();

 }

void analyze(TString inputfile, Int_t samp){	
	
	const TString inputFileTemp = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs/" + inputfile + ".root";
	
	inputfile = "Reading " + inputfile + " events... ";
	
	inputfile.Resize(60);
	
	cout << inputfile;

	// Set up storage variables
	Float_t genTau1_pt, genTau2_pt;
    Float_t genTau1_eta, genTau2_eta;
    Float_t genTau1_phi, genTau2_phi;
    Float_t genTau1_mass, genTau2_mass;

	// Charged pions coming from taus
    Float_t pions1_pt, pions2_pt;
    Float_t pions1_eta, pions2_eta;
    Float_t pions1_phi, pions2_phi;
    Float_t pions1_mass, pions2_mass;

    // Neutral pions coming from taus
    Float_t neutpions1_pt, neutpions2_pt;
    Float_t neutpions1_eta, neutpions2_eta;
    Float_t neutpions1_phi, neutpions2_phi;
    Float_t neutpions1_mass, neutpions2_mass;

    Float_t eventWeight;

	TLorentzVector vtau1, vtau2, vcpion1, vcpion2, vnpion1, vnpion2, vHiggs, vrho1, vrho2, vq1, vq2;
	TVector3 v3Higgs, vE1, vE2, vB1, vB2, v3tau1, v3tau2;
	UInt_t nProngTau1=0, nProngTau2=0, nLeptons=0;

	Float_t y1, y2;
	Float_t r=0.14;
	Float_t theta;

	TFile* infile = new TFile(inputFileTemp); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("nLeptons",		&nLeptons);
	intree->SetBranchAddress("eventWeight",		&eventWeight);
	intree->SetBranchAddress("nProngTau1",		&nProngTau1);
	intree->SetBranchAddress("nProngTau2",		&nProngTau2);
	intree->SetBranchAddress("genTau1_pt",		&genTau1_pt);
	intree->SetBranchAddress("genTau1_eta",		&genTau1_eta);
	intree->SetBranchAddress("genTau1_phi",		&genTau1_phi);
	intree->SetBranchAddress("genTau1_mass",	&genTau1_mass);
	intree->SetBranchAddress("genTau2_pt",		&genTau2_pt);
	intree->SetBranchAddress("genTau2_eta",		&genTau2_eta);
	intree->SetBranchAddress("genTau2_phi",		&genTau2_phi);
	intree->SetBranchAddress("genTau2_mass",	&genTau2_mass);
	intree->SetBranchAddress("pions1_pt",		&pions1_pt);
	intree->SetBranchAddress("pions1_eta",		&pions1_eta);
	intree->SetBranchAddress("pions1_phi",		&pions1_phi);
	intree->SetBranchAddress("pions1_mass",		&pions1_mass);
	intree->SetBranchAddress("pions2_pt",		&pions2_pt);
	intree->SetBranchAddress("pions2_eta",		&pions2_eta);
	intree->SetBranchAddress("pions2_phi",		&pions2_phi);
	intree->SetBranchAddress("pions2_mass",		&pions2_mass);
	intree->SetBranchAddress("neutpions1_pt",	&neutpions1_pt);
	intree->SetBranchAddress("neutpions1_eta",	&neutpions1_eta);
	intree->SetBranchAddress("neutpions1_phi",	&neutpions1_phi);
	intree->SetBranchAddress("neutpions1_mass",	&neutpions1_mass);
	intree->SetBranchAddress("neutpions2_pt",	&neutpions2_pt);
	intree->SetBranchAddress("neutpions2_eta",	&neutpions2_eta);
	intree->SetBranchAddress("neutpions2_phi",	&neutpions2_phi);
	intree->SetBranchAddress("neutpions2_mass",	&neutpions2_mass);

	Int_t nevents=0;

	for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) { // Event loop
		intree->GetEntry(iEntry);

		if(pions1_pt==0 || pions2_pt==0 || neutpions1_pt==0 || neutpions2_pt==0 || nProngTau1!=1 || nProngTau2!=1) continue;

		// Taus 4-vectors
		vtau1.SetPtEtaPhiM(genTau1_pt, genTau1_eta, genTau1_phi, genTau1_mass);
		vtau2.SetPtEtaPhiM(genTau2_pt, genTau2_eta, genTau2_phi, genTau2_mass);

		// Charged pions 4-vectors
		vcpion1.SetPtEtaPhiM(pions1_pt, pions1_eta, pions1_phi, pions1_mass);
		vcpion2.SetPtEtaPhiM(pions2_pt, pions2_eta, pions2_phi, pions2_mass);

		// Neutral pions 4-vectors
		vnpion1.SetPtEtaPhiM(neutpions1_pt, neutpions1_eta, neutpions1_phi, neutpions1_mass);
		vnpion2.SetPtEtaPhiM(neutpions2_pt, neutpions2_eta, neutpions2_phi, neutpions2_mass);

		// Rho 4-vectors
		vrho1=vcpion1+vnpion1;
		vrho2=vcpion2+vnpion2;

		// Higgs 4- and 3-vector
		vHiggs=vtau1+vtau2;
		v3Higgs.SetXYZ(-vHiggs.Px()/vHiggs.E(), -vHiggs.Py()/vHiggs.E(), -vHiggs.Pz()/vHiggs.E());

		vq1=vcpion1-vnpion1;
		vq2=vcpion2-vnpion2;

		y1=(vq1.Dot(vtau1))/(vrho1.Dot(vtau1));
		y2=(vq2.Dot(vtau2))/(vrho2.Dot(vtau2));

		// Boost to Higgs frame
		vcpion1.Boost(v3Higgs);
		vcpion2.Boost(v3Higgs);
		vnpion1.Boost(v3Higgs);
		vnpion2.Boost(v3Higgs);
		vtau1.Boost(v3Higgs);
		vtau2.Boost(v3Higgs);

		// Set up 3-vectors in Higgs frame
		TVector3 v3cpion1, v3cpion2, v3npion1, v3npion2;

		v3tau1.SetXYZ(vtau1.Px()/vtau1.E(), vtau1.Py()/vtau1.E(), vtau1.Pz()/vtau1.E());
		v3tau2.SetXYZ(vtau2.Px()/vtau2.E(), vtau2.Py()/vtau2.E(), vtau2.Pz()/vtau2.E());

		v3cpion1.SetXYZ(vcpion1.Px()/vcpion1.E(), vcpion1.Py()/vcpion1.E(), vcpion1.Pz()/vcpion1.E());
		v3cpion2.SetXYZ(vcpion2.Px()/vcpion2.E(), vcpion2.Py()/vcpion2.E(), vcpion2.Pz()/vcpion2.E());
		v3npion1.SetXYZ(vnpion1.Px()/vnpion1.E(), vnpion1.Py()/vnpion1.E(), vnpion1.Pz()/vnpion1.E());
		v3npion2.SetXYZ(vnpion2.Px()/vnpion2.E(), vnpion2.Py()/vnpion2.E(), vnpion2.Pz()/vnpion2.E());

		TVector3 tempE1, tempE2;

		tempE1=(y1-r)*v3cpion1-(y1+r)*v3npion1;
		tempE2=(y2-r)*v3cpion2-(y2+r)*v3npion2;

		vE1=tempE1-(tempE1.Dot(v3tau1))*v3tau1;
		vE2=tempE2-(tempE2.Dot(v3tau2))*v3tau2;

		vB1=v3tau1.Cross(vE1);
		vB2=v3tau2.Cross(vE2);

		// Compute theta variable
		theta=TMath::Sign(Double_t(1),v3tau1.Dot(vE2.Cross(vE1)))*TMath::ACos(vE1.Dot(vE2)/(vE1.Mag()*vE2.Mag()));

		// Fill histogram for theta variable
		(samp==0 ? hThetaB : (samp==1 ? hThetaS1 : hThetaS2))->Fill(theta);
		
		nevents++;
		
	} // end event loop
	
	TString out = "";
	out += nevents;
	out.Resize(8);
		
	cout << (true ? "\033[1;32m" : (out == "0       " ?  "\033[1;34m": "\033[1;31m")) << out << "\033[0m" << " events passed all cuts" << endl;
	
	infile->Close();
}

/*
 * FUNCTION FOR PRINTING AND SAVING THE RESULTS
 */

void saveResults(){

	TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);
	
	gStyle->SetOptStat(kFALSE);
	
	histogram(hThetaS1, hThetaS2, hThetaB, c1, "#Theta variable", "Fraction", "Theta.jpg");
	
	TFile f("histos.root","new");
	hThetaS1->Write();
	hThetaS2->Write();
	hThetaB->Write();
	f.Close();
	
	//cout << "\n\n\nProcess finished. Printing results...\n\n" << "\033[1;34mResults\033[0m\n\n";

	cout << "done\nExiting...\n\n\n";

}

/*
 * FUNCTION FOR SAVING THREE HISTOGRAMS
 */
 
void histogram(TH1D *histoS1, TH1D *histoS2, TH1D *histoB, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t nS1=1, nS2=1, nB=1;
	
	nS1/=histoS1->Integral();
	nS2/=histoS2->Integral();
	nB/=histoB->Integral();
	histoS1->Scale(nS1);
	histoS2->Scale(nS2);
	histoB->Scale(nB);
	
	Double_t max=histoB->GetMaximum();
	if(histoS1->GetMaximum() > max) max=histoS1->GetMaximum();
	if(histoS2->GetMaximum() > max) max=histoS2->GetMaximum();
	max*=1.1;
	
	histoS1->SetMaximum(max);
	histoS2->SetMaximum(max);
	histoB->SetMaximum(max);
	
	histoS1->SetLineWidth(3);
	histoS2->SetLineWidth(3);
	histoB->SetLineWidth(3);
	
	histoS1->Draw();
	// add axis labels
	histoS1->GetXaxis()->SetTitle(xTitle);
	histoS1->GetYaxis()->SetTitle(yTitle);
	histoS1->SetTitle(""); // title on top
	
	histoS2->SetLineColor(kGreen);
	histoS2->Draw("same");

	histoB->SetLineColor(kRed);
	histoB->Draw("same");
	
	TLegend *leg = new TLegend(0.605,0.675,0.885,0.875);
	leg->SetTextFont(72);
	leg->SetTextSize(0.04);
	leg->AddEntry(histoS1,"#Delta=0","l");
	leg->AddEntry(histoS2, "#Delta=#pi/2","l");
	leg->AddEntry(histoB, "Backgrounds","l");
	leg->Draw();

	can->SaveAs(name);
}

/*
 * FUNCTION FOR SAVING TWO HISTOGRAMS
 */
 
void histogram(TH1D *histoS, TH1D *histoB, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t nS=1, nB=1;
	
	nS/=histoS->Integral();
	nB/=histoB->Integral();
	histoS->Scale(nS);
	histoB->Scale(nB);
	
	Double_t max;
	if((histoS->GetMaximum()) > (histoB->GetMaximum())) max=1.1*(histoS->GetMaximum());
	else max=1.1*(histoB->GetMaximum());
	
	histoS->SetMaximum(max);
	histoB->SetMaximum(max);
	
	histoS->SetLineWidth(3);
	histoB->SetLineWidth(3);
	
	histoS->Draw();
	// add axis labels
	histoS->GetXaxis()->SetTitle(xTitle);
	histoS->GetYaxis()->SetTitle(yTitle);
	histoS->SetTitle(""); // title on top
	
	histoB->SetLineColor(kRed);
	histoB->Draw("same");
	
	TLegend *leg = new TLegend(0.6,0.65,0.88,0.85);
	leg->SetTextFont(72);
	leg->SetTextSize(0.04);
	leg->AddEntry(histoS,"#Delta=0","l");
	leg->AddEntry(histoB, "#Delta=#pi/2","l");
	leg->Draw();

	can->SaveAs(name);
}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, TCanvas *can, const char* xTitle, const char* yTitle, const char* name){
	Double_t norm=1;
	norm/=histo->Integral();
	histo->Scale(norm);
	histo->SetLineWidth(3);
	//histo->Fit("gaus","V");
	histo->Draw();
	// add axis labels
	histo->GetXaxis()->SetTitle(xTitle);
	histo->GetYaxis()->SetTitle(yTitle);
	histo->SetTitle(""); // title on top

	can->SaveAs(name);
}

/*
 * FUNCTION FOR dR calculation
 */
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 ) {

  const Float_t pi = 3.14159265358979;

  Float_t etaDiff = (eta1-eta2);
  Float_t phiDiff = fabs(phi1-phi2);
  while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

  Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}
