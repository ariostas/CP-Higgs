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

void cleanUpMergedFiles(TString sampName = "test") {

	// Set up input/output variables and files
	UInt_t numberOfEntries;

	//Event info
    Double_t eventWeight;
    Int_t NLeptons, NJets, NTauJets, hasH, sameCharge, NCPions1, NCPions2, NNPions1, NNPions2;

    //Charged pions
    Double_t CPion1_Pt, CPion2_Pt;
    Double_t CPion1_Eta, CPion2_Eta;
    Double_t CPion1_Phi, CPion2_Phi;
    Double_t CPion1_Mass, CPion2_Mass;
    Double_t CPions1_Pt, CPions2_Pt;
    Double_t CPions1_Eta, CPions2_Eta;
    Double_t CPions1_Phi, CPions2_Phi;
    Double_t CPions1_Mass, CPions2_Mass;

    //Neutral pions
    Double_t NPion11_Pt, NPion12_Pt;
    Double_t NPion11_Eta, NPion12_Eta;
    Double_t NPion11_Phi, NPion12_Phi;
    Double_t NPion11_Mass, NPion12_Mass;
    Double_t NPion21_Pt, NPion22_Pt;
    Double_t NPion21_Eta, NPion22_Eta;
    Double_t NPion21_Phi, NPion22_Phi;
    Double_t NPion21_Mass, NPion22_Mass;
    Double_t NPions1_Pt, NPions2_Pt;
    Double_t NPions1_Eta, NPions2_Eta;
    Double_t NPions1_Phi, NPions2_Phi;
    Double_t NPions1_Mass, NPions2_Mass;

    //Jets from Z
    Double_t ZJet1_Pt, ZJet2_Pt;
    Double_t ZJet1_Eta, ZJet2_Eta;
    Double_t ZJet1_Phi, ZJet2_Phi;
    Double_t ZJet1_Mass, ZJet2_Mass;

    //Electrons/Muons from Z
    Double_t ZLepton1_Pt, ZLepton2_Pt;
    Double_t ZLepton1_Eta, ZLepton2_Eta;
    Double_t ZLepton1_Phi, ZLepton2_Phi;
    Double_t ZLepton1_Mass, ZLepton2_Mass;

    //Reco Z
    Double_t ZReco_Pt, ZReco_Eta, ZReco_Phi, ZReco_Mass;

	TFile* infile = new TFile("/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/" + sampName + "_temp.root"); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("eventWeight",     &eventWeight);
    intree->SetBranchAddress("NLeptons",        &NLeptons);
    intree->SetBranchAddress("NJets",           &NJets);
    intree->SetBranchAddress("NTauJets",        &NTauJets);
    intree->SetBranchAddress("hasH",            &hasH);
    intree->SetBranchAddress("sameCharge",      &sameCharge);
    intree->SetBranchAddress("NCPions1",        &NCPions1);
    intree->SetBranchAddress("NCPions2",        &NCPions2);
    intree->SetBranchAddress("NNPions1",        &NNPions1);
    intree->SetBranchAddress("NNPions2",        &NNPions2);
    intree->SetBranchAddress("CPion1_Pt",       &CPion1_Pt);
    intree->SetBranchAddress("CPion1_Eta",      &CPion1_Eta);
    intree->SetBranchAddress("CPion1_Phi",      &CPion1_Phi);
    intree->SetBranchAddress("CPion1_Mass",     &CPion1_Mass);
    intree->SetBranchAddress("CPion2_Pt",       &CPion2_Pt);
    intree->SetBranchAddress("CPion2_Eta",      &CPion2_Eta);
    intree->SetBranchAddress("CPion2_Phi",      &CPion2_Phi);
    intree->SetBranchAddress("CPion2_Mass",     &CPion2_Mass);
    intree->SetBranchAddress("CPions1_Pt",      &CPions1_Pt);
    intree->SetBranchAddress("CPions1_Eta",     &CPions1_Eta);
    intree->SetBranchAddress("CPions1_Phi",     &CPions1_Phi);
    intree->SetBranchAddress("CPions1_Mass",    &CPions1_Mass);
    intree->SetBranchAddress("CPions2_Pt",      &CPions2_Pt);
    intree->SetBranchAddress("CPions2_Eta",     &CPions2_Eta);
    intree->SetBranchAddress("CPions2_Phi",     &CPions2_Phi);
    intree->SetBranchAddress("CPions2_Mass",    &CPions2_Mass);
    intree->SetBranchAddress("NPion11_Pt",      &NPion11_Pt);
    intree->SetBranchAddress("NPion11_Eta",     &NPion11_Eta);
    intree->SetBranchAddress("NPion11_Phi",     &NPion11_Phi);
    intree->SetBranchAddress("NPion11_Mass",    &NPion11_Mass);
    intree->SetBranchAddress("NPion12_Pt",      &NPion12_Pt);
    intree->SetBranchAddress("NPion12_Eta",     &NPion12_Eta);
    intree->SetBranchAddress("NPion12_Phi",     &NPion12_Phi);
    intree->SetBranchAddress("NPion12_Mass",    &NPion12_Mass);
    intree->SetBranchAddress("NPion21_Pt",      &NPion21_Pt);
    intree->SetBranchAddress("NPion21_Eta",     &NPion21_Eta);
    intree->SetBranchAddress("NPion21_Phi",     &NPion21_Phi);
    intree->SetBranchAddress("NPion21_Mass",    &NPion21_Mass);
    intree->SetBranchAddress("NPion22_Pt",      &NPion22_Pt);
    intree->SetBranchAddress("NPion22_Eta",     &NPion22_Eta);
    intree->SetBranchAddress("NPion22_Phi",     &NPion22_Phi);
    intree->SetBranchAddress("NPion22_Mass",    &NPion22_Mass);
    intree->SetBranchAddress("NPions1_Pt",      &NPions1_Pt);
    intree->SetBranchAddress("NPions1_Eta",     &NPions1_Eta);
    intree->SetBranchAddress("NPions1_Phi",     &NPions1_Phi);
    intree->SetBranchAddress("NPions1_Mass",    &NPions1_Mass);
    intree->SetBranchAddress("NPions2_Pt",      &NPions2_Pt);
    intree->SetBranchAddress("NPions2_Eta",     &NPions2_Eta);
    intree->SetBranchAddress("NPions2_Phi",     &NPions2_Phi);
    intree->SetBranchAddress("NPions2_Mass",    &NPions2_Mass);
    intree->SetBranchAddress("ZJet1_Pt",        &ZJet1_Pt);
    intree->SetBranchAddress("ZJet1_Eta",       &ZJet1_Eta);
    intree->SetBranchAddress("ZJet1_Phi",       &ZJet1_Phi);
    intree->SetBranchAddress("ZJet1_Mass",      &ZJet1_Mass);
    intree->SetBranchAddress("ZJet2_Pt",        &ZJet2_Pt);
    intree->SetBranchAddress("ZJet2_Eta",       &ZJet2_Eta);
    intree->SetBranchAddress("ZJet2_Phi",       &ZJet2_Phi);
    intree->SetBranchAddress("ZJet2_Mass",      &ZJet2_Mass);
    intree->SetBranchAddress("ZLepton1_Pt",     &ZLepton1_Pt);
    intree->SetBranchAddress("ZLepton1_Eta",    &ZLepton1_Eta);
    intree->SetBranchAddress("ZLepton1_Phi",    &ZLepton1_Phi);
    intree->SetBranchAddress("ZLepton1_Mass",   &ZLepton1_Mass);
    intree->SetBranchAddress("ZLepton2_Pt",     &ZLepton2_Pt);
    intree->SetBranchAddress("ZLepton2_Eta",    &ZLepton2_Eta);
    intree->SetBranchAddress("ZLepton2_Phi",    &ZLepton2_Phi);
    intree->SetBranchAddress("ZLepton2_Mass",   &ZLepton2_Mass);
    intree->SetBranchAddress("ZReco_Pt",        &ZReco_Pt);
    intree->SetBranchAddress("ZReco_Eta",       &ZReco_Eta);
    intree->SetBranchAddress("ZReco_Phi",       &ZReco_Phi);
    intree->SetBranchAddress("ZReco_Mass",      &ZReco_Mass);

	TTree* infotree = (TTree*) infile->Get("Count"); assert(infotree);
	infotree->SetBranchAddress("numberOfEntries",      &numberOfEntries);

	UInt_t totalnEvents=0;
	
	for (UInt_t iEntry=0; iEntry<infotree->GetEntries(); iEntry++) {
		infotree->GetEntry(iEntry);
		totalnEvents+=numberOfEntries;
	}

	TFile *outFile = new TFile("/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/" + sampName + ".root", "RECREATE");
	
	// tree to hold the number of events
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("TotalEvents",		&totalnEvents,		"TotalEvents/i");
    sampTree->Fill();

	// tree to hold information about selected events
	TTree *outtree = new TTree("Events", "Events");
	
    outtree->Branch("eventWeight",      &eventWeight,       "eventWeight/D");
    outtree->Branch("NLeptons",         &NLeptons,          "NLeptons/I");
    outtree->Branch("NJets",         &NJets,          "NJets/I");
    outtree->Branch("NTauJets",         &NTauJets,          "NTauJets/I");
    outtree->Branch("hasH",         &hasH,          "hasH/I");
    outtree->Branch("sameCharge",         &sameCharge,          "sameCharge/I");

    outtree->Branch("NCPions1",         &NCPions1,          "NCPions1/I");
    outtree->Branch("NCPions2",         &NCPions2,          "NCPions2/I");
    outtree->Branch("NNPions1",         &NNPions1,          "NNPions1/I");
    outtree->Branch("NNPions2",         &NNPions2,          "NNPions2/I");

    outtree->Branch("CPion1_Pt",       &CPion1_Pt,        "CPion1_Pt/D");
    outtree->Branch("CPion1_Eta",      &CPion1_Eta,       "CPion1_Eta/D");
    outtree->Branch("CPion1_Phi",      &CPion1_Phi,       "CPion1_Phi/D");
    outtree->Branch("CPion1_Mass",     &CPion1_Mass,      "CPion1_Mass/D");
    outtree->Branch("CPion2_Pt",       &CPion2_Pt,        "CPion2_Pt/D");
    outtree->Branch("CPion2_Eta",      &CPion2_Eta,       "CPion2_Eta/D");
    outtree->Branch("CPion2_Phi",      &CPion2_Phi,       "CPion2_Phi/D");
    outtree->Branch("CPion2_Mass",     &CPion2_Mass,      "CPion2_Mass/D");
    outtree->Branch("CPions1_Pt",       &CPions1_Pt,        "CPions1_Pt/D");
    outtree->Branch("CPions1_Eta",      &CPions1_Eta,       "CPions1_Eta/D");
    outtree->Branch("CPions1_Phi",      &CPions1_Phi,       "CPions1_Phi/D");
    outtree->Branch("CPions1_Mass",     &CPions1_Mass,      "CPions1_Mass/D");
    outtree->Branch("CPions2_Pt",       &CPions2_Pt,        "CPions2_Pt/D");
    outtree->Branch("CPions2_Eta",      &CPions2_Eta,       "CPions2_Eta/D");
    outtree->Branch("CPions2_Phi",      &CPions2_Phi,       "CPions2_Phi/D");
    outtree->Branch("CPions2_Mass",     &CPions2_Mass,      "CPions2_Mass/D");

    outtree->Branch("NPion11_Pt",       &NPion11_Pt,        "NPion11_Pt/D");
    outtree->Branch("NPion11_Eta",      &NPion11_Eta,       "NPion11_Eta/D");
    outtree->Branch("NPion11_Phi",      &NPion11_Phi,       "NPion11_Phi/D");
    outtree->Branch("NPion11_Mass",     &NPion11_Mass,      "NPion11_Mass/D");
    outtree->Branch("NPion12_Pt",       &NPion12_Pt,        "NPion12_Pt/D");
    outtree->Branch("NPion12_Eta",      &NPion12_Eta,       "NPion12_Eta/D");
    outtree->Branch("NPion12_Phi",      &NPion12_Phi,       "NPion12_Phi/D");
    outtree->Branch("NPion12_Mass",     &NPion12_Mass,      "NPion12_Mass/D");
    outtree->Branch("NPion21_Pt",       &NPion21_Pt,        "NPion21_Pt/D");
    outtree->Branch("NPion21_Eta",      &NPion21_Eta,       "NPion21_Eta/D");
    outtree->Branch("NPion21_Phi",      &NPion21_Phi,       "NPion21_Phi/D");
    outtree->Branch("NPion21_Mass",     &NPion21_Mass,      "NPion21_Mass/D");
    outtree->Branch("NPion22_Pt",       &NPion22_Pt,        "NPion22_Pt/D");
    outtree->Branch("NPion22_Eta",      &NPion22_Eta,       "NPion22_Eta/D");
    outtree->Branch("NPion22_Phi",      &NPion22_Phi,       "NPion22_Phi/D");
    outtree->Branch("NPion22_Mass",     &NPion22_Mass,      "NPion22_Mass/D");
    outtree->Branch("NPions1_Pt",       &NPions1_Pt,        "NPions1_Pt/D");
    outtree->Branch("NPions1_Eta",      &NPions1_Eta,       "NPions1_Eta/D");
    outtree->Branch("NPions1_Phi",      &NPions1_Phi,       "NPions1_Phi/D");
    outtree->Branch("NPions1_Mass",     &NPions1_Mass,      "NPions1_Mass/D");
    outtree->Branch("NPions2_Pt",       &NPions2_Pt,        "NPions2_Pt/D");
    outtree->Branch("NPions2_Eta",      &NPions2_Eta,       "NPions2_Eta/D");
    outtree->Branch("NPions2_Phi",      &NPions2_Phi,       "NPions2_Phi/D");
    outtree->Branch("NPions2_Mass",     &NPions2_Mass,      "NPions2_Mass/D");

    outtree->Branch("ZJet1_Pt",       &ZJet1_Pt,        "ZJet1_Pt/D");
    outtree->Branch("ZJet1_Eta",      &ZJet1_Eta,       "ZJet1_Eta/D");
    outtree->Branch("ZJet1_Phi",      &ZJet1_Phi,       "ZJet1_Phi/D");
    outtree->Branch("ZJet1_Mass",     &ZJet1_Mass,      "ZJet1_Mass/D");
    outtree->Branch("ZJet2_Pt",       &ZJet2_Pt,        "ZJet2_Pt/D");
    outtree->Branch("ZJet2_Eta",      &ZJet2_Eta,       "ZJet2_Eta/D");
    outtree->Branch("ZJet2_Phi",      &ZJet2_Phi,       "ZJet2_Phi/D");
    outtree->Branch("ZJet2_Mass",     &ZJet2_Mass,      "ZJet2_Mass/D");

    outtree->Branch("ZLepton1_Pt",       &ZLepton1_Pt,        "ZLepton1_Pt/D");
    outtree->Branch("ZLepton1_Eta",      &ZLepton1_Eta,       "ZLepton1_Eta/D");
    outtree->Branch("ZLepton1_Phi",      &ZLepton1_Phi,       "ZLepton1_Phi/D");
    outtree->Branch("ZLepton1_Mass",     &ZLepton1_Mass,      "ZLepton1_Mass/D");
    outtree->Branch("ZLepton2_Pt",       &ZLepton2_Pt,        "ZLepton2_Pt/D");
    outtree->Branch("ZLepton2_Eta",      &ZLepton2_Eta,       "ZLepton2_Eta/D");
    outtree->Branch("ZLepton2_Phi",      &ZLepton2_Phi,       "ZLepton2_Phi/D");
    outtree->Branch("ZLepton2_Mass",     &ZLepton2_Mass,      "ZLepton2_Mass/D");

    outtree->Branch("ZReco_Pt",       &ZReco_Pt,        "ZReco_Pt/D");
    outtree->Branch("ZReco_Eta",      &ZReco_Eta,       "ZReco_Eta/D");
    outtree->Branch("ZReco_Phi",      &ZReco_Phi,       "ZReco_Phi/D");
    outtree->Branch("ZReco_Mass",     &ZReco_Mass,      "ZReco_Mass/D");

	for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) {
		intree->GetEntry(iEntry);
		eventWeight/=Double_t(totalnEvents);
		outtree->Fill();
	}
	
	outFile->Write();
	outFile->Close();
  
	cout << "Finished " << sampName << " with " << totalnEvents << " total events" << endl;

}
