#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TMath.h>
#include <TChain.h>
#include <TH1.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Math/LorentzVector.h"
#endif

void cleanUpMergedFiles(TString infilename="/filepath/test.root", TString outfilename="test.root") {

	// Set up input/output variables and files
	UInt_t n=0;
	Float_t eventWeight;

    UInt_t nProngTau1=0, nProngTau2=0, nLeptons=0, zToLep=0;

    // Jets matched to gen taus
    Float_t jetTau1_pt, jetTau2_pt;
    Float_t jetTau1_eta, jetTau2_eta;
    Float_t jetTau1_phi, jetTau2_phi;
    Float_t jetTau1_mass, jetTau2_mass;

    // Gen taus
    Float_t genTau1_pt, genTau2_pt;
    Float_t genTau1_eta, genTau2_eta;
    Float_t genTau1_phi, genTau2_phi;
    Float_t genTau1_mass, genTau2_mass;

    // Visible taus (taus minus neutrino)
    Float_t visTau1_pt, visTau2_pt;
    Float_t visTau1_eta, visTau2_eta;
    Float_t visTau1_phi, visTau2_phi;
    Float_t visTau1_mass, visTau2_mass;

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

    // Particles from Z
    Float_t z1_pt, z2_pt;
    Float_t z1_eta, z2_eta;
    Float_t z1_phi, z2_phi;
    Float_t z1_mass, z2_mass;

    // Jets matched to particles from Z
    Float_t jetz1_pt, jetz2_pt;
    Float_t jetz1_eta, jetz2_eta;
    Float_t jetz1_phi, jetz2_phi;
    Float_t jetz1_mass, jetz2_mass;

    // Z
    Float_t z_pt, z_eta, z_phi, z_mass;

	TFile* infile = new TFile(infilename); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("eventWeight",		&eventWeight);
	intree->SetBranchAddress("nProngTau1",		&nProngTau1);
	intree->SetBranchAddress("nProngTau2",		&nProngTau2);
	intree->SetBranchAddress("jetTau1_pt",		&jetTau1_pt);
	intree->SetBranchAddress("jetTau1_eta",		&jetTau1_eta);
	intree->SetBranchAddress("jetTau1_phi",		&jetTau1_phi);
	intree->SetBranchAddress("jetTau1_mass",	&jetTau1_mass);
	intree->SetBranchAddress("jetTau2_pt",		&jetTau2_pt);
	intree->SetBranchAddress("jetTau2_eta",		&jetTau2_eta);
	intree->SetBranchAddress("jetTau2_phi",		&jetTau2_phi);
	intree->SetBranchAddress("jetTau2_mass",	&jetTau2_mass);
	intree->SetBranchAddress("genTau1_pt",		&genTau1_pt);
	intree->SetBranchAddress("genTau1_eta",		&genTau1_eta);
	intree->SetBranchAddress("genTau1_phi",		&genTau1_phi);
	intree->SetBranchAddress("genTau1_mass",	&genTau1_mass);
	intree->SetBranchAddress("genTau2_pt",		&genTau2_pt);
	intree->SetBranchAddress("genTau2_eta",		&genTau2_eta);
	intree->SetBranchAddress("genTau2_phi",		&genTau2_phi);
	intree->SetBranchAddress("genTau2_mass",	&genTau2_mass);
	intree->SetBranchAddress("visTau1_pt",		&visTau1_pt);
	intree->SetBranchAddress("visTau1_eta",		&visTau1_eta);
	intree->SetBranchAddress("visTau1_phi",		&visTau1_phi);
	intree->SetBranchAddress("visTau1_mass",	&visTau1_mass);
	intree->SetBranchAddress("visTau2_pt",		&visTau2_pt);
	intree->SetBranchAddress("visTau2_eta",		&visTau2_eta);
	intree->SetBranchAddress("visTau2_phi",		&visTau2_phi);
	intree->SetBranchAddress("visTau2_mass",	&visTau2_mass);
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
	intree->SetBranchAddress("z1_pt",			&z1_pt);
	intree->SetBranchAddress("z1_eta",			&z1_eta);
	intree->SetBranchAddress("z1_phi",			&z1_phi);
	intree->SetBranchAddress("z1_mass",			&z1_mass);
	intree->SetBranchAddress("z2_pt",			&z2_pt);
	intree->SetBranchAddress("z2_eta",			&z2_eta);
	intree->SetBranchAddress("z2_phi",			&z2_phi);
	intree->SetBranchAddress("z2_mass",			&z2_mass);
	intree->SetBranchAddress("jetz1_pt",		&jetz1_pt);
	intree->SetBranchAddress("jetz1_eta",		&jetz1_eta);
	intree->SetBranchAddress("jetz1_phi",		&jetz1_phi);
	intree->SetBranchAddress("jetz1_mass",		&jetz1_mass);
	intree->SetBranchAddress("jetz2_pt",		&jetz2_pt);
	intree->SetBranchAddress("jetz2_eta",		&jetz2_eta);
	intree->SetBranchAddress("jetz2_phi",		&jetz2_phi);
	intree->SetBranchAddress("jetz2_mass",		&jetz2_mass);
    intree->SetBranchAddress("z_pt",            &z_pt);
    intree->SetBranchAddress("z_eta",           &z_eta);
    intree->SetBranchAddress("z_phi",           &z_phi);
    intree->SetBranchAddress("z_mass",          &z_mass);

	TTree* infotree = (TTree*) infile->Get("Count"); assert(infotree);
	infotree->SetBranchAddress("n",      &n);
	

	UInt_t totalnEvents=0;
	
	for (UInt_t iEntry=0; iEntry<infotree->GetEntries(); iEntry++) {
		infotree->GetEntry(iEntry);
		totalnEvents+=n;
	}

	TFile *outFile = new TFile(outfilename, "RECREATE");
	
	// tree to hold the number of events
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("nEvents",		&n,		"nEvents/i");

	// tree to hold information about selected events
	TTree *outTree = new TTree("Events", "Events");
	outTree->Branch("eventWeight",      &eventWeight,       "eventWeight/f");   // event weight from cross-section and Event->Weight

    outTree->Branch("nLeptons",         &nLeptons,          "nLeptons/i");      // number of leptons
    outTree->Branch("nProngTau1",       &nProngTau1,        "nProngTau1/i");
    outTree->Branch("nProngTau2",       &nProngTau2,        "nProngTau2/i");

    outTree->Branch("jetTau1_pt",       &jetTau1_pt,        "jetTau1_pt/f");    // pt(Tau1)
    outTree->Branch("jetTau1_eta",      &jetTau1_eta,       "jetTau1_eta/f");   // eta(Tau1)
    outTree->Branch("jetTau1_phi",      &jetTau1_phi,       "jetTau1_phi/f");   // phi(Tau1)
    outTree->Branch("jetTau1_mass",     &jetTau1_mass,      "jetTau1_mass/f");  // m(Tau1)

    outTree->Branch("jetTau2_pt",       &jetTau2_pt,        "jetTau2_pt/f");    // pt(Tau2)
    outTree->Branch("jetTau2_eta",      &jetTau2_eta,       "jetTau2_eta/f");   // eta(Tau2)
    outTree->Branch("jetTau2_phi",      &jetTau2_phi,       "jetTau2_phi/f");   // phi(Tau2)
    outTree->Branch("jetTau2_mass",     &jetTau2_mass,      "jetTau2_mass/f");  // m(Tau2)

    outTree->Branch("genTau1_pt",       &genTau1_pt,        "genTau1_pt/f");    // pt(Tau1)
    outTree->Branch("genTau1_eta",      &genTau1_eta,       "genTau1_eta/f");   // eta(Tau1)
    outTree->Branch("genTau1_phi",      &genTau1_phi,       "genTau1_phi/f");   // phi(Tau1)
    outTree->Branch("genTau1_mass",     &genTau1_mass,      "genTau1_mass/f");  // m(Tau1)

    outTree->Branch("genTau2_pt",       &genTau2_pt,        "genTau2_pt/f");    // pt(Tau2)
    outTree->Branch("genTau2_eta",      &genTau2_eta,       "genTau2_eta/f");   // eta(Tau2)
    outTree->Branch("genTau2_phi",      &genTau2_phi,       "genTau2_phi/f");   // phi(Tau2)
    outTree->Branch("genTau2_mass",     &genTau2_mass,      "genTau2_mass/f");  // m(Tau2)

    outTree->Branch("visTau1_pt",       &visTau1_pt,        "visTau1_pt/f");    // pt(Tau1)
    outTree->Branch("visTau1_eta",      &visTau1_eta,       "visTau1_eta/f");   // eta(Tau1)
    outTree->Branch("visTau1_phi",      &visTau1_phi,       "visTau1_phi/f");   // phi(Tau1)
    outTree->Branch("visTau1_mass",     &visTau1_mass,      "visTau1_mass/f");  // m(Tau1)

    outTree->Branch("visTau2_pt",       &visTau2_pt,        "visTau2_pt/f");    // pt(Tau2)
    outTree->Branch("visTau2_eta",      &visTau2_eta,       "visTau2_eta/f");   // eta(Tau2)
    outTree->Branch("visTau2_phi",      &visTau2_phi,       "visTau2_phi/f");   // phi(Tau2)
    outTree->Branch("visTau2_mass",     &visTau2_mass,      "visTau2_mass/f");  // m(Tau2)

    outTree->Branch("pions1_pt",        &pions1_pt,          "pions1_pt/f");
    outTree->Branch("pions1_eta",       &pions1_eta,         "pions1_eta/f");
    outTree->Branch("pions1_phi",       &pions1_phi,         "pions1_phi/f");
    outTree->Branch("pions1_mass",      &pions1_mass,        "pions1_mass/f");

    outTree->Branch("pions2_pt",        &pions2_pt,          "pions2_pt/f");
    outTree->Branch("pions2_eta",       &pions2_eta,         "pions2_eta/f");
    outTree->Branch("pions2_phi",       &pions2_phi,         "pions2_phi/f");
    outTree->Branch("pions2_mass",      &pions2_mass,        "pions2_mass/f");

    outTree->Branch("neutpions1_pt",    &neutpions1_pt,      "neutpions1_pt/f");
    outTree->Branch("neutpions1_eta",   &neutpions1_eta,     "neutpions1_eta/f");
    outTree->Branch("neutpions1_phi",   &neutpions1_phi,     "neutpions1_phi/f");
    outTree->Branch("neutpions1_mass",  &neutpions1_mass,    "neutpions1_mass/f");

    outTree->Branch("neutpions2_pt",    &neutpions2_pt,      "neutpions2_pt/f");
    outTree->Branch("neutpions2_eta",   &neutpions2_eta,     "neutpions2_eta/f");
    outTree->Branch("neutpions2_phi",   &neutpions2_phi,     "neutpions2_phi/f");
    outTree->Branch("neutpions2_mass",  &neutpions2_mass,    "neutpions2_mass/f");

    outTree->Branch("z1_pt",    		&z1_pt,      		 "z1_pt/f");
    outTree->Branch("z1_eta",   		&z1_eta,     		 "z1_eta/f");
    outTree->Branch("z1_phi",   		&z1_phi,     		 "z1_phi/f");
    outTree->Branch("z1_mass",  		&z1_mass,    		 "z1_mass/f");

    outTree->Branch("z2_pt",    		&z2_pt,              "z2_pt/f");
    outTree->Branch("z2_eta",   		&z2_eta,             "z2_eta/f");
    outTree->Branch("z2_phi",   		&z2_phi,             "z2_phi/f");
    outTree->Branch("z2_mass",  		&z2_mass,            "z2_mass/f");

 	outTree->Branch("jetz1_pt",         &jetz1_pt,           "jetz1_pt/f");
    outTree->Branch("jetz1_eta",        &jetz1_eta,          "jetz1_eta/f");
    outTree->Branch("jetz1_phi",        &jetz1_phi,          "jetz1_phi/f");
    outTree->Branch("jetz1_mass",       &jetz1_mass,         "jetz1_mass/f");

    outTree->Branch("jetz2_pt",         &jetz2_pt,           "jetz2_pt/f");
    outTree->Branch("jetz2_eta",        &jetz2_eta,          "jetz2_eta/f");
    outTree->Branch("jetz2_phi",        &jetz2_phi,          "jetz2_phi/f");
    outTree->Branch("jetz2_mass",       &jetz2_mass,         "jetz2_mass/f");

    outTree->Branch("z_pt",            &z_pt,              "z_pt/f");
    outTree->Branch("z_eta",           &z_eta,             "z_eta/f");
    outTree->Branch("z_phi",           &z_phi,             "z_phi/f");
    outTree->Branch("z_mass",          &z_mass,            "z_mass/f");
	
	sampTree->Fill();

	for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) {
		intree->GetEntry(iEntry);
		eventWeight/=Double_t(totalnEvents);

		outTree->Fill();

	}
	
	outFile->Write();
	outFile->Close();
  
	cout << "Finished " << outfilename << " with " << totalnEvents << " total events" << endl;

}
