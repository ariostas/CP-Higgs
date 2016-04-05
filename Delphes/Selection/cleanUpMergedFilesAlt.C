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

void cleanUpMergedFilesAlt(TString sampName = "test") {

	// Set up input/output variables and files
	UInt_t numberOfEntries;

	//Event info
    Double_t eventWeight;
    Int_t NLeptons, NJets, NTauJets, hasH, ZFromLep, oppositeTrackCharge;

    //Charged pions
    Double_t CPion1_Pt  , CPion2_Pt  ;
    Double_t CPion1_Eta , CPion2_Eta ;
    Double_t CPion1_Phi , CPion2_Phi ;
    Double_t CPion1_Mass, CPion2_Mass;

    //Neutral pions
    Double_t NPion1_Pt  , NPion2_Pt  ;
    Double_t NPion1_Eta , NPion2_Eta ;
    Double_t NPion1_Phi , NPion2_Phi ;
    Double_t NPion1_Mass, NPion2_Mass;

    //Pion jets
    Double_t JetTau1_Pt  , JetTau2_Pt  ;
    Double_t JetTau1_Eta , JetTau2_Eta ;
    Double_t JetTau1_Phi , JetTau2_Phi ;
    Double_t JetTau1_Mass, JetTau2_Mass;

    //Remaining tracks and photons
    Double_t TracksTau1_Pt  , TracksTau2_Pt  ;
    Double_t TracksTau1_Eta , TracksTau2_Eta ;
    Double_t TracksTau1_Phi , TracksTau2_Phi ;
    Double_t TracksTau1_Mass, TracksTau2_Mass;
    Double_t PhotonsTau1_Pt  , PhotonsTau2_Pt  ;
    Double_t PhotonsTau1_Eta , PhotonsTau2_Eta ;
    Double_t PhotonsTau1_Phi , PhotonsTau2_Phi ;
    Double_t PhotonsTau1_Mass, PhotonsTau2_Mass;

    //Jets/leptons from Z
    Double_t ZParticle1_Pt  , ZParticle2_Pt  ;
    Double_t ZParticle1_Eta , ZParticle2_Eta ;
    Double_t ZParticle1_Phi , ZParticle2_Phi ;
    Double_t ZParticle1_Mass, ZParticle2_Mass;

    //Reconstructed Z from all other particles
    Double_t ZReco_Pt, ZReco_Eta, ZReco_Phi, ZReco_Mass;

	TFile* infile = new TFile("/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_smallAlt/" + sampName + "_temp.root"); assert(infile);
	TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

	intree->SetBranchAddress("eventWeight",     &eventWeight);
    intree->SetBranchAddress("NLeptons",        &NLeptons);
    intree->SetBranchAddress("NJets",           &NJets);
    intree->SetBranchAddress("NTauJets",        &NTauJets);
    intree->SetBranchAddress("hasH",            &hasH);
    intree->SetBranchAddress("ZFromLep",        &ZFromLep);
    intree->SetBranchAddress("oppositeTrackCharge",&oppositeTrackCharge);
    intree->SetBranchAddress("CPion1_Pt",       &CPion1_Pt);
    intree->SetBranchAddress("CPion1_Eta",      &CPion1_Eta);
    intree->SetBranchAddress("CPion1_Phi",      &CPion1_Phi);
    intree->SetBranchAddress("CPion1_Mass",     &CPion1_Mass);
    intree->SetBranchAddress("CPion2_Pt",       &CPion2_Pt);
    intree->SetBranchAddress("CPion2_Eta",      &CPion2_Eta);
    intree->SetBranchAddress("CPion2_Phi",      &CPion2_Phi);
    intree->SetBranchAddress("CPion2_Mass",     &CPion2_Mass);
    intree->SetBranchAddress("NPion1_Pt",       &NPion1_Pt);
    intree->SetBranchAddress("NPion1_Eta",      &NPion1_Eta);
    intree->SetBranchAddress("NPion1_Phi",      &NPion1_Phi);
    intree->SetBranchAddress("NPion1_Mass",     &NPion1_Mass);
    intree->SetBranchAddress("NPion2_Pt",       &NPion2_Pt);
    intree->SetBranchAddress("NPion2_Eta",      &NPion2_Eta);
    intree->SetBranchAddress("NPion2_Phi",      &NPion2_Phi);
    intree->SetBranchAddress("NPion2_Mass",     &NPion2_Mass);
    intree->SetBranchAddress("JetTau1_Pt",      &JetTau1_Pt);
    intree->SetBranchAddress("JetTau1_Eta",     &JetTau1_Eta);
    intree->SetBranchAddress("JetTau1_Phi",     &JetTau1_Phi);
    intree->SetBranchAddress("JetTau1_Mass",    &JetTau1_Mass);
    intree->SetBranchAddress("JetTau2_Pt",      &JetTau2_Pt);
    intree->SetBranchAddress("JetTau2_Eta",     &JetTau2_Eta);
    intree->SetBranchAddress("JetTau2_Phi",     &JetTau2_Phi);
    intree->SetBranchAddress("JetTau2_Mass",    &JetTau2_Mass);
    intree->SetBranchAddress("TracksTau1_Pt",   &TracksTau1_Pt);
    intree->SetBranchAddress("TracksTau1_Eta",  &TracksTau1_Eta);
    intree->SetBranchAddress("TracksTau1_Phi",  &TracksTau1_Phi);
    intree->SetBranchAddress("TracksTau1_Mass", &TracksTau1_Mass);
    intree->SetBranchAddress("TracksTau2_Pt",   &TracksTau2_Pt);
    intree->SetBranchAddress("TracksTau2_Eta",  &TracksTau2_Eta);
    intree->SetBranchAddress("TracksTau2_Phi",  &TracksTau2_Phi);
    intree->SetBranchAddress("TracksTau2_Mass", &TracksTau2_Mass);
    intree->SetBranchAddress("PhotonsTau1_Pt",  &PhotonsTau1_Pt);
    intree->SetBranchAddress("PhotonsTau1_Eta", &PhotonsTau1_Eta);
    intree->SetBranchAddress("PhotonsTau1_Phi", &PhotonsTau1_Phi);
    intree->SetBranchAddress("PhotonsTau1_Mass",&PhotonsTau1_Mass);
    intree->SetBranchAddress("PhotonsTau2_Pt",  &PhotonsTau2_Pt);
    intree->SetBranchAddress("PhotonsTau2_Eta", &PhotonsTau2_Eta);
    intree->SetBranchAddress("PhotonsTau2_Phi", &PhotonsTau2_Phi);
    intree->SetBranchAddress("PhotonsTau2_Mass",&PhotonsTau2_Mass);
    intree->SetBranchAddress("ZParticle1_Pt",   &ZParticle1_Pt);
    intree->SetBranchAddress("ZParticle1_Eta",  &ZParticle1_Eta);
    intree->SetBranchAddress("ZParticle1_Phi",  &ZParticle1_Phi);
    intree->SetBranchAddress("ZParticle1_Mass", &ZParticle1_Mass);
    intree->SetBranchAddress("ZParticle2_Pt",   &ZParticle2_Pt);
    intree->SetBranchAddress("ZParticle2_Eta",  &ZParticle2_Eta);
    intree->SetBranchAddress("ZParticle2_Phi",  &ZParticle2_Phi);
    intree->SetBranchAddress("ZParticle2_Mass", &ZParticle2_Mass);
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

	TFile *outFile = new TFile("/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_smallAlt/" + sampName + ".root", "RECREATE");
	
	// tree to hold the number of events
	TTree *sampTree = new TTree("Info", "Info");
	sampTree->Branch("TotalEvents",		&totalnEvents,		"TotalEvents/i");
    sampTree->Fill();

	// tree to hold information about selected events
	TTree *outtree = new TTree("Events", "Events");
	
    outtree->Branch("eventWeight",      &eventWeight,       "eventWeight/D");
    outtree->Branch("NLeptons",         &NLeptons,          "NLeptons/I");
    outtree->Branch("NJets",            &NJets,             "NJets/I");
    outtree->Branch("NTauJets",         &NTauJets,          "NTauJets/I");
    outtree->Branch("hasH",             &hasH,              "hasH/I");
    outtree->Branch("ZFromLep",         &ZFromLep,          "ZFromLep/I");
    outtree->Branch("oppositeTrackCharge",&oppositeTrackCharge,"oppositeTrackCharge/I");

    outtree->Branch("CPion1_Pt",        &CPion1_Pt,         "CPion1_Pt/D");
    outtree->Branch("CPion1_Eta",       &CPion1_Eta,        "CPion1_Eta/D");
    outtree->Branch("CPion1_Phi",       &CPion1_Phi,        "CPion1_Phi/D");
    outtree->Branch("CPion1_Mass",      &CPion1_Mass,       "CPion1_Mass/D");
    outtree->Branch("CPion2_Pt",        &CPion2_Pt,         "CPion2_Pt/D");
    outtree->Branch("CPion2_Eta",       &CPion2_Eta,        "CPion2_Eta/D");
    outtree->Branch("CPion2_Phi",       &CPion2_Phi,        "CPion2_Phi/D");
    outtree->Branch("CPion2_Mass",      &CPion2_Mass,       "CPion2_Mass/D");

    outtree->Branch("NPion1_Pt",        &NPion1_Pt,         "NPion1_Pt/D");
    outtree->Branch("NPion1_Eta",       &NPion1_Eta,        "NPion1_Eta/D");
    outtree->Branch("NPion1_Phi",       &NPion1_Phi,        "NPion1_Phi/D");
    outtree->Branch("NPion1_Mass",      &NPion1_Mass,       "NPion1_Mass/D");
    outtree->Branch("NPion2_Pt",        &NPion2_Pt,         "NPion2_Pt/D");
    outtree->Branch("NPion2_Eta",       &NPion2_Eta,        "NPion2_Eta/D");
    outtree->Branch("NPion2_Phi",       &NPion2_Phi,        "NPion2_Phi/D");
    outtree->Branch("NPion2_Mass",      &NPion2_Mass,       "NPion2_Mass/D");

    outtree->Branch("JetTau1_Pt",       &JetTau1_Pt,        "JetTau1_Pt/D");
    outtree->Branch("JetTau1_Eta",      &JetTau1_Eta,       "JetTau1_Eta/D");
    outtree->Branch("JetTau1_Phi",      &JetTau1_Phi,       "JetTau1_Phi/D");
    outtree->Branch("JetTau1_Mass",     &JetTau1_Mass,      "JetTau1_Mass/D");
    outtree->Branch("JetTau2_Pt",       &JetTau2_Pt,        "JetTau2_Pt/D");
    outtree->Branch("JetTau2_Eta",      &JetTau2_Eta,       "JetTau2_Eta/D");
    outtree->Branch("JetTau2_Phi",      &JetTau2_Phi,       "JetTau2_Phi/D");
    outtree->Branch("JetTau2_Mass",     &JetTau2_Mass,      "JetTau2_Mass/D");

    outtree->Branch("TracksTau1_Pt",    &TracksTau1_Pt,     "TracksTau1_Pt/D");
    outtree->Branch("TracksTau1_Eta",   &TracksTau1_Eta,    "TracksTau1_Eta/D");
    outtree->Branch("TracksTau1_Phi",   &TracksTau1_Phi,    "TracksTau1_Phi/D");
    outtree->Branch("TracksTau1_Mass",  &TracksTau1_Mass,   "TracksTau1_Mass/D");
    outtree->Branch("TracksTau2_Pt",    &TracksTau2_Pt,     "TracksTau2_Pt/D");
    outtree->Branch("TracksTau2_Eta",   &TracksTau2_Eta,    "TracksTau2_Eta/D");
    outtree->Branch("TracksTau2_Phi",   &TracksTau2_Phi,    "TracksTau2_Phi/D");
    outtree->Branch("TracksTau2_Mass",  &TracksTau2_Mass,   "TracksTau2_Mass/D");
    outtree->Branch("PhotonsTau1_Pt",   &PhotonsTau1_Pt,    "PhotonsTau1_Pt/D");
    outtree->Branch("PhotonsTau1_Eta",  &PhotonsTau1_Eta,   "PhotonsTau1_Eta/D");
    outtree->Branch("PhotonsTau1_Phi",  &PhotonsTau1_Phi,   "PhotonsTau1_Phi/D");
    outtree->Branch("PhotonsTau1_Mass", &PhotonsTau1_Mass,  "PhotonsTau1_Mass/D");
    outtree->Branch("PhotonsTau2_Pt",   &PhotonsTau2_Pt,    "PhotonsTau2_Pt/D");
    outtree->Branch("PhotonsTau2_Eta",  &PhotonsTau2_Eta,   "PhotonsTau2_Eta/D");
    outtree->Branch("PhotonsTau2_Phi",  &PhotonsTau2_Phi,   "PhotonsTau2_Phi/D");
    outtree->Branch("PhotonsTau2_Mass", &PhotonsTau2_Mass,  "PhotonsTau2_Mass/D");

    outtree->Branch("ZParticle1_Pt",    &ZParticle1_Pt,     "ZParticle1_Pt/D");
    outtree->Branch("ZParticle1_Eta",   &ZParticle1_Eta,    "ZParticle1_Eta/D");
    outtree->Branch("ZParticle1_Phi",   &ZParticle1_Phi,    "ZParticle1_Phi/D");
    outtree->Branch("ZParticle1_Mass",  &ZParticle1_Mass,   "ZParticle1_Mass/D");
    outtree->Branch("ZParticle2_Pt",    &ZParticle2_Pt,     "ZParticle2_Pt/D");
    outtree->Branch("ZParticle2_Eta",   &ZParticle2_Eta,    "ZParticle2_Eta/D");
    outtree->Branch("ZParticle2_Phi",   &ZParticle2_Phi,    "ZParticle2_Phi/D");
    outtree->Branch("ZParticle2_Mass",  &ZParticle2_Mass,   "ZParticle2_Mass/D");

    outtree->Branch("ZReco_Pt",         &ZReco_Pt,          "ZReco_Pt/D");
    outtree->Branch("ZReco_Eta",        &ZReco_Eta,         "ZReco_Eta/D");
    outtree->Branch("ZReco_Phi",        &ZReco_Phi,         "ZReco_Phi/D");
    outtree->Branch("ZReco_Mass",       &ZReco_Mass,        "ZReco_Mass/D");

	for(UInt_t iEntry=0; iEntry<intree->GetEntries(); iEntry++) {
		intree->GetEntry(iEntry);
		eventWeight/=Double_t(totalnEvents);
		outtree->Fill();
	}
	
	outFile->Write();
	outFile->Close();
  
	cout << "Finished " << sampName << " with " << totalnEvents << " total events" << endl;

}
