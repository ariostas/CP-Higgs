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
	Double_t eventWeight;

    UInt_t ncpions1=0, ncpions2=0, nnpions1=0, nnpions2=0, nLeptons=0;
    Int_t zToLep=0;

    // Jets matched to gen taus
    Double_t jetTau1_pt, jetTau2_pt;
    Double_t jetTau1_eta, jetTau2_eta;
    Double_t jetTau1_phi, jetTau2_phi;
    Double_t jetTau1_mass, jetTau2_mass;

    // Gen taus
    Double_t genTau1_pt, genTau2_pt;
    Double_t genTau1_eta, genTau2_eta;
    Double_t genTau1_phi, genTau2_phi;
    Double_t genTau1_mass, genTau2_mass;

    // Visible taus (taus minus neutrino)
    Double_t visTau1_pt, visTau2_pt;
    Double_t visTau1_eta, visTau2_eta;
    Double_t visTau1_phi, visTau2_phi;
    Double_t visTau1_mass, visTau2_mass;

    // Charged pions coming from taus
    Double_t cpions1_pt, cpions2_pt;
    Double_t cpions1_eta, cpions2_eta;
    Double_t cpions1_phi, cpions2_phi;
    Double_t cpions1_mass, cpions2_mass;

    // Neutral pions coming from taus
    Double_t npions1_pt, npions2_pt;
    Double_t npions1_eta, npions2_eta;
    Double_t npions1_phi, npions2_phi;
    Double_t npions1_mass, npions2_mass;

    // Particles/jets from Z
    Double_t z1_pt, z2_pt;
    Double_t z1_eta, z2_eta;
    Double_t z1_phi, z2_phi;
    Double_t z1_mass, z2_mass;

    // Z
    Double_t z_pt, z_eta, z_phi, z_mass;

	TFile* infile = new TFile(infilename); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("eventWeight",	   &eventWeight);
	intree->SetBranchAddress("ncpions1",       &ncpions1);
	intree->SetBranchAddress("ncpions2",       &ncpions2);
    intree->SetBranchAddress("nnpions1",       &nnpions1);
    intree->SetBranchAddress("nnpions2",       &nnpions2);
    intree->SetBranchAddress("zToLep",         &zToLep);
	intree->SetBranchAddress("jetTau1_pt",	   &jetTau1_pt);
	intree->SetBranchAddress("jetTau1_eta",	   &jetTau1_eta);
	intree->SetBranchAddress("jetTau1_phi",	   &jetTau1_phi);
	intree->SetBranchAddress("jetTau1_mass",   &jetTau1_mass);
	intree->SetBranchAddress("jetTau2_pt",	   &jetTau2_pt);
	intree->SetBranchAddress("jetTau2_eta",	   &jetTau2_eta);
	intree->SetBranchAddress("jetTau2_phi",	   &jetTau2_phi);
	intree->SetBranchAddress("jetTau2_mass",   &jetTau2_mass);
	intree->SetBranchAddress("genTau1_pt",	   &genTau1_pt);
	intree->SetBranchAddress("genTau1_eta",	   &genTau1_eta);
	intree->SetBranchAddress("genTau1_phi",	   &genTau1_phi);
	intree->SetBranchAddress("genTau1_mass",   &genTau1_mass);
	intree->SetBranchAddress("genTau2_pt",	   &genTau2_pt);
	intree->SetBranchAddress("genTau2_eta",	   &genTau2_eta);
	intree->SetBranchAddress("genTau2_phi",	   &genTau2_phi);
	intree->SetBranchAddress("genTau2_mass",   &genTau2_mass);
	intree->SetBranchAddress("visTau1_pt",	   &visTau1_pt);
	intree->SetBranchAddress("visTau1_eta",	   &visTau1_eta);
	intree->SetBranchAddress("visTau1_phi",	   &visTau1_phi);
	intree->SetBranchAddress("visTau1_mass",   &visTau1_mass);
	intree->SetBranchAddress("visTau2_pt",	   &visTau2_pt);
	intree->SetBranchAddress("visTau2_eta",	   &visTau2_eta);
	intree->SetBranchAddress("visTau2_phi",	   &visTau2_phi);
	intree->SetBranchAddress("visTau2_mass",   &visTau2_mass);
	intree->SetBranchAddress("cpions1_pt",	   &cpions1_pt);
	intree->SetBranchAddress("cpions1_eta",	   &cpions1_eta);
	intree->SetBranchAddress("cpions1_phi",	   &cpions1_phi);
	intree->SetBranchAddress("cpions1_mass",   &cpions1_mass);
	intree->SetBranchAddress("cpions2_pt",	   &cpions2_pt);
	intree->SetBranchAddress("cpions2_eta",	   &cpions2_eta);
	intree->SetBranchAddress("cpions2_phi",	   &cpions2_phi);
	intree->SetBranchAddress("cpions2_mass",   &cpions2_mass);
	intree->SetBranchAddress("npions1_pt",	   &npions1_pt);
	intree->SetBranchAddress("npions1_eta",	   &npions1_eta);
	intree->SetBranchAddress("npions1_phi",	   &npions1_phi);
	intree->SetBranchAddress("npions1_mass",   &npions1_mass);
	intree->SetBranchAddress("npions2_pt",	   &npions2_pt);
	intree->SetBranchAddress("npions2_eta",	   &npions2_eta);
	intree->SetBranchAddress("npions2_phi",	   &npions2_phi);
	intree->SetBranchAddress("npions2_mass",   &npions2_mass);
	intree->SetBranchAddress("z1_pt",		   &z1_pt);
	intree->SetBranchAddress("z1_eta",		   &z1_eta);
	intree->SetBranchAddress("z1_phi",		   &z1_phi);
	intree->SetBranchAddress("z1_mass",		   &z1_mass);
	intree->SetBranchAddress("z2_pt",		   &z2_pt);
	intree->SetBranchAddress("z2_eta",		   &z2_eta);
	intree->SetBranchAddress("z2_phi",		   &z2_phi);
	intree->SetBranchAddress("z2_mass",		   &z2_mass);
    intree->SetBranchAddress("z_pt",           &z_pt);
    intree->SetBranchAddress("z_eta",          &z_eta);
    intree->SetBranchAddress("z_phi",          &z_phi);
    intree->SetBranchAddress("z_mass",         &z_mass);

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
	TTree *outtree = new TTree("Events", "Events");
	outtree->Branch("eventWeight",      &eventWeight,       "eventWeight/D");   // event weight from cross-section and Event->Weight

    outtree->Branch("nLeptons",         &nLeptons,          "nLeptons/i");      // number of leptons
    outtree->Branch("ncpions1",         &ncpions1,          "ncpions1/i");
    outtree->Branch("ncpions2",         &ncpions2,          "ncpions2/i");
    outtree->Branch("nnpions1",         &nnpions1,          "nnpions1/i");
    outtree->Branch("nnpions2",         &nnpions2,          "nnpions2/i");

    outtree->Branch("zToLep",           &zToLep,            "zToLep/I");

    outtree->Branch("jetTau1_pt",       &jetTau1_pt,        "jetTau1_pt/D");    // pt(Tau1)
    outtree->Branch("jetTau1_eta",      &jetTau1_eta,       "jetTau1_eta/D");   // eta(Tau1)
    outtree->Branch("jetTau1_phi",      &jetTau1_phi,       "jetTau1_phi/D");   // phi(Tau1)
    outtree->Branch("jetTau1_mass",     &jetTau1_mass,      "jetTau1_mass/D");  // m(Tau1)

    outtree->Branch("jetTau2_pt",       &jetTau2_pt,        "jetTau2_pt/D");    // pt(Tau2)
    outtree->Branch("jetTau2_eta",      &jetTau2_eta,       "jetTau2_eta/D");   // eta(Tau2)
    outtree->Branch("jetTau2_phi",      &jetTau2_phi,       "jetTau2_phi/D");   // phi(Tau2)
    outtree->Branch("jetTau2_mass",     &jetTau2_mass,      "jetTau2_mass/D");  // m(Tau2)

    outtree->Branch("genTau1_pt",       &genTau1_pt,        "genTau1_pt/D");    // pt(Tau1)
    outtree->Branch("genTau1_eta",      &genTau1_eta,       "genTau1_eta/D");   // eta(Tau1)
    outtree->Branch("genTau1_phi",      &genTau1_phi,       "genTau1_phi/D");   // phi(Tau1)
    outtree->Branch("genTau1_mass",     &genTau1_mass,      "genTau1_mass/D");  // m(Tau1)

    outtree->Branch("genTau2_pt",       &genTau2_pt,        "genTau2_pt/D");    // pt(Tau2)
    outtree->Branch("genTau2_eta",      &genTau2_eta,       "genTau2_eta/D");   // eta(Tau2)
    outtree->Branch("genTau2_phi",      &genTau2_phi,       "genTau2_phi/D");   // phi(Tau2)
    outtree->Branch("genTau2_mass",     &genTau2_mass,      "genTau2_mass/D");  // m(Tau2)

    outtree->Branch("visTau1_pt",       &visTau1_pt,        "visTau1_pt/D");    // pt(Tau1)
    outtree->Branch("visTau1_eta",      &visTau1_eta,       "visTau1_eta/D");   // eta(Tau1)
    outtree->Branch("visTau1_phi",      &visTau1_phi,       "visTau1_phi/D");   // phi(Tau1)
    outtree->Branch("visTau1_mass",     &visTau1_mass,      "visTau1_mass/D");  // m(Tau1)

    outtree->Branch("visTau2_pt",       &visTau2_pt,        "visTau2_pt/D");    // pt(Tau2)
    outtree->Branch("visTau2_eta",      &visTau2_eta,       "visTau2_eta/D");   // eta(Tau2)
    outtree->Branch("visTau2_phi",      &visTau2_phi,       "visTau2_phi/D");   // phi(Tau2)
    outtree->Branch("visTau2_mass",     &visTau2_mass,      "visTau2_mass/D");  // m(Tau2)

    outtree->Branch("cpions1_pt",       &cpions1_pt,        "cpions1_pt/D");
    outtree->Branch("cpions1_eta",      &cpions1_eta,       "cpions1_eta/D");
    outtree->Branch("cpions1_phi",      &cpions1_phi,       "cpions1_phi/D");
    outtree->Branch("cpions1_mass",     &cpions1_mass,      "cpions1_mass/D");

    outtree->Branch("cpions2_pt",       &cpions2_pt,        "pions2_pt/D");
    outtree->Branch("cpions2_eta",      &cpions2_eta,       "pions2_eta/D");
    outtree->Branch("cpions2_phi",      &cpions2_phi,       "pions2_phi/D");
    outtree->Branch("cpions2_mass",     &cpions2_mass,      "pions2_mass/D");

    outtree->Branch("npions1_pt",       &npions1_pt,        "npions1_pt/D");
    outtree->Branch("npions1_eta",      &npions1_eta,       "npions1_eta/D");
    outtree->Branch("npions1_phi",      &npions1_phi,       "npions1_phi/D");
    outtree->Branch("npions1_mass",     &npions1_mass,      "npions1_mass/D");

    outtree->Branch("npions2_pt",       &npions2_pt,        "npions2_pt/D");
    outtree->Branch("npions2_eta",      &npions2_eta,       "npions2_eta/D");
    outtree->Branch("npions2_phi",      &npions2_phi,       "npions2_phi/D");
    outtree->Branch("npions2_mass",     &npions2_mass,      "npions2_mass/D");

    outtree->Branch("z1_pt",            &z1_pt,             "z1_pt/D");
    outtree->Branch("z1_eta",           &z1_eta,            "z1_eta/D");
    outtree->Branch("z1_phi",           &z1_phi,            "z1_phi/D");
    outtree->Branch("z1_mass",          &z1_mass,           "z1_mass/D");

    outtree->Branch("z2_pt",            &z2_pt,             "z2_pt/D");
    outtree->Branch("z2_eta",           &z2_eta,            "z2_eta/D");
    outtree->Branch("z2_phi",           &z2_phi,            "z2_phi/D");
    outtree->Branch("z2_mass",          &z2_mass,           "z2_mass/D");

    outtree->Branch("z_pt",             &z_pt,              "z_pt/D");
    outtree->Branch("z_eta",            &z_eta,             "z_eta/D");
    outtree->Branch("z_phi",            &z_phi,             "z_phi/D");
    outtree->Branch("z_mass",           &z_mass,            "z_mass/D");
	
	sampTree->Fill();

	for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) {
		intree->GetEntry(iEntry);
		eventWeight/=Double_t(totalnEvents);

		outtree->Fill();

	}
	
	outFile->Write();
	outFile->Close();
  
	cout << "Finished " << outfilename << " with " << totalnEvents << " total events" << endl;

}
