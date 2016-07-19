#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <fstream>
#include <TString.h>
#include <TLorentzVector.h>
#include "TMath.h"

// Delphes libraries
#include "modules/Delphes.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"

#endif

using namespace TMath;
using namespace std;

void select(const TString TempInput = "zh_delphes_0pi12_1.root", const Double_t XSec = 1){

    // Read input input file
    TChain chain("Delphes");
    chain.Add("root://eoscms//store/user/arapyan/mc/Delphes/" + TempInput);
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t TotalEntries = treeReader->GetEntries();
    cout << "Reading " << "root://eoscms//store/user/arapyan/mc/Delphes/" << TempInput << " ..." << endl;

    // Set up branches to read in from file
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchChargedHadron = treeReader->UseBranch("Track");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    TClonesArray *branchNeutralHadron = treeReader->UseBranch("Tower");
    
    // Check if the file was opened correctly
    if (!branchMuon || !branchElectron || !branchPhoton || !branchJet || 
        !branchParticle || !branchChargedHadron || !branchNeutralHadron){
        cout << "Error opening file. Exiting..." << endl;
        return;
    }
    
    // Set up storage variables
    Jet *ZJet1 = 0, *ZJet2 = 0;
    Electron *ZElectron1 = 0, *ZElectron2 = 0;
    Muon *ZMuon1 = 0, *ZMuon2 = 0;
    Track *CPion1 = 0, *CPion2 = 0;
    Photon *NPion11 = 0, *NPion12 = 0, *NPion21 = 0, *NPion22 = 0;
    TLorentzVector ZReco;

    // Set up output file and trees
    TString output = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/out_" + TempInput;
    TFile *outfile = new TFile(output, "RECREATE");

    TTree *infoTree = new TTree("Count", "Count");
    infoTree->Branch("TotalEntries", &TotalEntries, "TotalEntries/i");
    infoTree->Fill();

    TTree *outtree = new TTree("Events", "Events");

    // Event info
    Double_t eventWeight;
    Int_t NLeptons, NJets, NTauJets, hasH, sameCharge;

    // Charged pions
    Double_t CPion1_Pt, CPion2_Pt;
    Double_t CPion1_Eta, CPion2_Eta;
    Double_t CPion1_Phi, CPion2_Phi;
    Double_t CPion1_Mass, CPion2_Mass;

    // Neutral pions (2 photons each)
    Double_t NPion11_Pt, NPion12_Pt;
    Double_t NPion11_Eta, NPion12_Eta;
    Double_t NPion11_Phi, NPion12_Phi;
    Double_t NPion11_Mass, NPion12_Mass;
    Double_t NPion21_Pt, NPion22_Pt;
    Double_t NPion21_Eta, NPion22_Eta;
    Double_t NPion21_Phi, NPion22_Phi;
    Double_t NPion21_Mass, NPion22_Mass;

    // Jets from Z
    Double_t ZJet1_Pt, ZJet2_Pt;
    Double_t ZJet1_Eta, ZJet2_Eta;
    Double_t ZJet1_Phi, ZJet2_Phi;
    Double_t ZJet1_Mass, ZJet2_Mass;

    // Electrons/Muons from Z
    Double_t ZLepton1_Pt, ZLepton2_Pt;
    Double_t ZLepton1_Eta, ZLepton2_Eta;
    Double_t ZLepton1_Phi, ZLepton2_Phi;
    Double_t ZLepton1_Mass, ZLepton2_Mass;

    // Reco Z (from adding all particles apart from the tau decay products)
    Double_t ZReco_Pt, ZReco_Eta, ZReco_Phi, ZReco_Mass;

    outtree->Branch("eventWeight",      &eventWeight,       "eventWeight/D");
    outtree->Branch("NLeptons",         &NLeptons,          "NLeptons/I");
    outtree->Branch("NJets",            &NJets,             "NJets/I");
    outtree->Branch("NTauJets",         &NTauJets,          "NTauJets/I");
    outtree->Branch("hasH",             &hasH,              "hasH/I");
    outtree->Branch("sameCharge",       &sameCharge,        "sameCharge/I");

    outtree->Branch("CPion1_Pt",        &CPion1_Pt,         "CPion1_Pt/D");
    outtree->Branch("CPion1_Eta",       &CPion1_Eta,        "CPion1_Eta/D");
    outtree->Branch("CPion1_Phi",       &CPion1_Phi,        "CPion1_Phi/D");
    outtree->Branch("CPion1_Mass",      &CPion1_Mass,       "CPion1_Mass/D");
    outtree->Branch("CPion2_Pt",        &CPion2_Pt,         "CPion2_Pt/D");
    outtree->Branch("CPion2_Eta",       &CPion2_Eta,        "CPion2_Eta/D");
    outtree->Branch("CPion2_Phi",       &CPion2_Phi,        "CPion2_Phi/D");
    outtree->Branch("CPion2_Mass",      &CPion2_Mass,       "CPion2_Mass/D");

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

    outtree->Branch("ZJet1_Pt",         &ZJet1_Pt,          "ZJet1_Pt/D");
    outtree->Branch("ZJet1_Eta",        &ZJet1_Eta,         "ZJet1_Eta/D");
    outtree->Branch("ZJet1_Phi",        &ZJet1_Phi,         "ZJet1_Phi/D");
    outtree->Branch("ZJet1_Mass",       &ZJet1_Mass,        "ZJet1_Mass/D");
    outtree->Branch("ZJet2_Pt",         &ZJet2_Pt,          "ZJet2_Pt/D");
    outtree->Branch("ZJet2_Eta",        &ZJet2_Eta,         "ZJet2_Eta/D");
    outtree->Branch("ZJet2_Phi",        &ZJet2_Phi,         "ZJet2_Phi/D");
    outtree->Branch("ZJet2_Mass",       &ZJet2_Mass,        "ZJet2_Mass/D");

    outtree->Branch("ZLepton1_Pt",      &ZLepton1_Pt,       "ZLepton1_Pt/D");
    outtree->Branch("ZLepton1_Eta",     &ZLepton1_Eta,      "ZLepton1_Eta/D");
    outtree->Branch("ZLepton1_Phi",     &ZLepton1_Phi,      "ZLepton1_Phi/D");
    outtree->Branch("ZLepton1_Mass",    &ZLepton1_Mass,     "ZLepton1_Mass/D");
    outtree->Branch("ZLepton2_Pt",      &ZLepton2_Pt,       "ZLepton2_Pt/D");
    outtree->Branch("ZLepton2_Eta",     &ZLepton2_Eta,      "ZLepton2_Eta/D");
    outtree->Branch("ZLepton2_Phi",     &ZLepton2_Phi,      "ZLepton2_Phi/D");
    outtree->Branch("ZLepton2_Mass",    &ZLepton2_Mass,     "ZLepton2_Mass/D");

    outtree->Branch("ZReco_Pt",         &ZReco_Pt,          "ZReco_Pt/D");
    outtree->Branch("ZReco_Eta",        &ZReco_Eta,         "ZReco_Eta/D");
    outtree->Branch("ZReco_Phi",        &ZReco_Phi,         "ZReco_Phi/D");
    outtree->Branch("ZReco_Mass",       &ZReco_Mass,        "ZReco_Mass/D");

    for (Int_t iEntry=0; iEntry<chain.GetEntries(); iEntry++){
        
        // Read entry
        treeReader->ReadEntry(iEntry);

        // Reset variables
        NLeptons = NJets = NTauJets = hasH = sameCharge = 0;
        ZJet1 = ZJet2 = 0;
        ZElectron1 = ZElectron2 = 0;
        ZMuon1 = ZMuon2 = 0;
        CPion1 = CPion2 = 0;
        NPion11 = NPion12 = NPion21 = NPion22 = 0;
        ZReco.SetPxPyPzE(0,0,0,0);

        // Selection starts here

        // Select leading two electrons
        for (Int_t i = 0; i < branchElectron->GetEntries(); i++) {
            Electron *electron = (Electron*) branchElectron->At(i);
            
            if(electron->P4().P() < 10. || electron->IsolationVar > 0.4 || Abs(electron->Eta) > 2.5) continue;
            
            NLeptons++;

            if(!ZElectron1){
                ZElectron1 = electron;
            }
            else if(electron->P4().P() > ZElectron1->P4().P()){
                ZElectron2 = ZElectron1;
                ZElectron1 = electron;
            }
            else if(!ZElectron2){
                ZElectron2 = electron;
            }
            else if(electron->P4().P() > ZElectron2->P4().P()){
                ZElectron2 = electron;
            }

        }

        // Discard electron if only one is selected
        if(!ZElectron1 || !ZElectron2) ZElectron1 = ZElectron2 = 0;

        // Select leading two muon
        for (Int_t i = 0; i < branchMuon->GetEntries(); i++) {
            Muon *muon = (Muon*) branchMuon->At(i);
            
            if(muon->P4().P() < 10. || muon->IsolationVar > 0.4 || Abs(muon->Eta) > 2.5) continue;
            
            NLeptons++;

            if(!ZMuon1){
                ZMuon1 = muon;
            }
            else if(muon->P4().P() > ZMuon1->P4().P()){
                ZMuon2 = ZMuon1;
                ZMuon1 = muon;
            }
            else if(!ZMuon2){
                ZMuon2 = muon;
            }
            else if(muon->P4().P() > ZMuon2->P4().P()){
                ZMuon2 = muon;
            }

        }

        // Discard muon if only one if selected
        if(!ZMuon1 || !ZMuon2) ZMuon1 = ZMuon2 = 0;

        // Order charged hardons (far from an isolated lepton) with respect to momentum
        vector<Track*> orderedChargedHadrons;
        for(Int_t n = 0; n < branchChargedHadron->GetEntries(); n++) {
            
            Track *track = 0;

            for(Int_t m = 0; m < branchChargedHadron->GetEntries(); m++) {

                Track *chargedHadron = (Track*) branchChargedHadron->At(m);

                if(chargedHadron->P4().P() < 0.5 || Abs(chargedHadron->Eta) > 2.5) continue;

                bool isDifferent = true;
                for(UInt_t x = 0; x < orderedChargedHadrons.size(); x++){
                    if(chargedHadron == orderedChargedHadrons.at(x)){isDifferent = false; break;}
                }
                if(!isDifferent) continue;

                bool closeToLep = false;
                for (Int_t i = 0; i < branchElectron->GetEntries(); i++) {
                    Electron *electron = (Electron*) branchElectron->At(i);
            
                    if(electron->P4().P() < 10. || electron->IsolationVar > 0.4) continue;
                    if(electron->P4().DeltaR(chargedHadron->P4()) < 0.5){closeToLep = true; break;}
                }

                for (Int_t i = 0; i <branchMuon->GetEntries(); i++) {
                    Muon *muon = (Muon*) branchMuon->At(i);
            
                    if(muon->P4().P() < 10. || muon->IsolationVar > 0.4) continue;
                    if(muon->P4().DeltaR(chargedHadron->P4()) < 0.5){closeToLep = true; break;}
                }

                if(closeToLep) continue;
                
                if(!track) track = chargedHadron;
                else if(chargedHadron->P4().P() > track->P4().P()) track = chargedHadron;

            }

            if(track) orderedChargedHadrons.push_back(track);
            else break;

        }

        // Look for tracks that are likely from a tau decay of interest (no other tracks and exactly two photons within cone of 0.4)
        vector<Track*> tauCandidates;
        for(UInt_t i = 0; i < orderedChargedHadrons.size(); i++){

            Track *track = orderedChargedHadrons.at(i);

            Int_t isolated = 0;
            for(Int_t n = 0; n < branchChargedHadron->GetEntries(); n++){
                Track *testTrack = (Track*) branchChargedHadron->At(n);
                if(track->P4().DeltaR(testTrack->P4()) < 0.4) isolated++;
            }

            Int_t pi0 = 0;
            for(Int_t n = 0; n < branchPhoton->GetEntries(); n++){
                Photon *photon = (Photon*) branchPhoton->At(n);
                if(photon->P4().P() < 0.5) continue;
                if(track->P4().DeltaR(photon->P4()) < 0.4) pi0++;
            }

            if(isolated == 1 && pi0 == 2) tauCandidates.push_back(track);

        }

        // Require at least two tau candidates
        if(tauCandidates.size() < 2) continue;

        // Select most likely taus from tau candidates
        // It was found that the number of tau candidates is generally low, so the two with highest pt are selected
        // Additionally opposite charge and a minimum separation is required, an appropriate value was determined using gen-level info
        for(UInt_t i = 0; i < tauCandidates.size(); i++){
            
            for(UInt_t j = i+1; j < tauCandidates.size(); j++) {
                
                Track *tempTrack1 = tauCandidates.at(i), *tempTrack2 = tauCandidates.at(j);
                
                if(tempTrack1->Charge == tempTrack2->Charge) continue;
                if(tempTrack1->P4().DeltaR(tempTrack2->P4()) < 2.25) continue;

                if(tempTrack1->Charge > 0){
                    CPion1 = tempTrack1;
                    CPion2 = tempTrack2;
                }
                else{
                    CPion1 = tempTrack2;
                    CPion2 = tempTrack1;
                }
                break;
            }
        }

        // Check whether there are were two selected tau candidates
        if(!CPion1 || !CPion2) continue;

        // Store the corresponding pair of photons from each tau
        for (Int_t i = 0; i < branchPhoton->GetEntries(); i++) {
            Photon *photon = (Photon*) branchPhoton->At(i);
            
            if(photon->P4().P() < 0.5 || Abs(photon->Eta) > 2.5) continue;

            if(photon->P4().DeltaR(CPion1->P4()) < 0.4){
                if(!NPion11){
                    NPion11 = photon;
                }
                else if(!NPion12){
                    NPion12 = photon;
                }
            }

            else if(photon->P4().DeltaR(CPion2->P4()) < 0.4){
                if(!NPion21){
                    NPion21 = photon;
                }
                else if(!NPion22){
                    NPion22 = photon;
                }
            }

        }

        // Check whether there are two photons from each tau
        // This should already be the case, but is is checked to avoid crashes
        if(!NPion11 || !NPion12 || !NPion21 || !NPion22){cout << "Error finding photons, check the code" << endl; return;}

        // Check whether the event contains a Higgs
        // This is done because some background events might contain a Higgs, so they are actually signal events
        for(Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) {
            GenParticle *particle = (GenParticle*) branchParticle->At(iParticle);
            if(abs(particle->PID)==25){hasH = 1; break;}
        }

        // Find the two jets that reconstruct a mass closest to the Z mass
        // It was found that a better selection is achieved by using a mass of 85 GeV
        Double_t minDeltaM = 999;
        for (Int_t i = 0; i < branchJet->GetEntries(); i++) {
            Jet *jet1 = (Jet*) branchJet->At(i);
            
            if(jet1->P4().P() < 20. || Abs(jet1->Eta) > 2.5) continue;
            
            NJets++;
            if(jet1->TauTag) NTauJets++;

            for(Int_t j = i+1; j < branchJet->GetEntries(); j++){
                Jet *jet2 = (Jet*) branchJet->At(j);

                if(jet2->P4().P() < 20. || Abs(jet2->Eta) > 2.5) continue;
                if(jet1->P4().DeltaR(jet2->P4()) < 0.4) continue;

                TLorentzVector vTempZ = jet1->P4() + jet2->P4();

                Double_t tempDeltaM = fabs(vTempZ.M() - 85. /*or 91.2*/);

                if(tempDeltaM < minDeltaM){
                    minDeltaM = tempDeltaM;
                    ZJet1 = jet1;
                    ZJet2 = jet2;
                }
            }
        }

        // As an alternative method of reconstructing the Z, all the objects that were not taken as coming
        // from the Higgs are added up

        for(Int_t i = 0; i < branchEFlowTrack->GetEntries(); i++) {
            Track *track = (Track*) branchEFlowTrack->At(i);

            if(track->P4().DeltaR(CPion1->P4()) < 0.4) continue;
            if(track->P4().DeltaR(CPion2->P4()) < 0.4) continue;

            ZReco += track->P4();
        }

        for(Int_t i = 0; i < branchEFlowNeutralHadron->GetEntries(); i++) {
            Tower *neutralHadron = (Tower*) branchEFlowNeutralHadron->At(i);

            if(neutralHadron->P4().DeltaR(CPion1->P4()) < 0.4) continue;
            if(neutralHadron->P4().DeltaR(CPion2->P4()) < 0.4) continue;

            ZReco += neutralHadron->P4();

        }

        for(Int_t i = 0; i < branchPhoton->GetEntries(); i++) {
            Photon *photon = (Photon*) branchPhoton->At(i);

            if(photon->P4().DeltaR(CPion1->P4()) < 0.4) continue;
            if(photon->P4().DeltaR(CPion2->P4()) < 0.4) continue;

            ZReco += photon->P4();
        }

        // Selection ends here

        // Set event weight (in this case all events are weighted equally)
        eventWeight = 1000000.*XSec; // Cross section in pb

        // Set all tree variables
        CPion1_Pt   = CPion1->PT;
        CPion1_Eta  = CPion1->Eta;
        CPion1_Phi  = CPion1->Phi;
        CPion1_Mass = 0.139570;
        CPion2_Pt   = CPion2->PT;
        CPion2_Eta  = CPion2->Eta;
        CPion2_Phi  = CPion2->Phi;
        CPion2_Mass = 0.139570;

        NPion11_Pt   = NPion11->PT;
        NPion11_Eta  = NPion11->Eta;
        NPion11_Phi  = NPion11->Phi;
        NPion11_Mass = NPion11->P4().M();
        NPion12_Pt   = NPion12->PT;
        NPion12_Eta  = NPion12->Eta;
        NPion12_Phi  = NPion12->Phi;
        NPion12_Mass = NPion12->P4().M();
        NPion21_Pt   = NPion21->PT;
        NPion21_Eta  = NPion21->Eta;
        NPion21_Phi  = NPion21->Phi;
        NPion21_Mass = NPion21->P4().M();
        NPion22_Pt   = NPion22->PT;
        NPion22_Eta  = NPion22->Eta;
        NPion22_Phi  = NPion22->Phi;
        NPion22_Mass = NPion22->P4().M();

        if(ZJet1 && ZJet2){
            ZJet1_Pt   = ZJet1->PT;
            ZJet1_Eta  = ZJet1->Eta;
            ZJet1_Phi  = ZJet1->Phi;
            ZJet1_Mass = ZJet1->Mass;
            ZJet2_Pt   = ZJet2->PT;
            ZJet2_Eta  = ZJet2->Eta;
            ZJet2_Phi  = ZJet2->Phi;
            ZJet2_Mass = ZJet2->Mass;
        }
        else{
            ZJet1_Pt   = 0.;
            ZJet1_Eta  = 0.;
            ZJet1_Phi  = 0.;
            ZJet1_Mass = 0.;
            ZJet2_Pt   = 0.;
            ZJet2_Eta  = 0.;
            ZJet2_Phi  = 0.;
            ZJet2_Mass = 0.;
        }

        if(ZElectron1 && ZElectron2){
            ZLepton1_Pt   = ZElectron1->PT;
            ZLepton1_Eta  = ZElectron1->Eta;
            ZLepton1_Phi  = ZElectron1->Phi;
            ZLepton1_Mass = 0.000510999;
            ZLepton2_Pt   = ZElectron2->PT;
            ZLepton2_Eta  = ZElectron2->Eta;
            ZLepton2_Phi  = ZElectron2->Phi;
            ZLepton2_Mass = 0.000510999;
            sameCharge = (ZElectron1->Charge == ZElectron2->Charge ? 1 : 0);
        }
        else if(ZMuon1 && ZMuon2){
            if(!ZElectron1 || !ZElectron2){
                ZLepton1_Pt   = ZMuon1->PT;
                ZLepton1_Eta  = ZMuon1->Eta;
                ZLepton1_Phi  = ZMuon1->Phi;
                ZLepton1_Mass = 0.105658;
                ZLepton2_Pt   = ZMuon2->PT;
                ZLepton2_Eta  = ZMuon2->Eta;
                ZLepton2_Phi  = ZMuon2->Phi;
                ZLepton2_Mass = 0.105658;
                sameCharge = (ZMuon1->Charge == ZMuon2->Charge ? 1 : 0);
            }
            else if(fabs((ZMuon1->P4() + ZMuon2->P4()).M() - 91.2) < fabs((ZElectron1->P4() + ZElectron2->P4()).M() - 91.2)){
                ZLepton1_Pt   = ZMuon1->PT;
                ZLepton1_Eta  = ZMuon1->Eta;
                ZLepton1_Phi  = ZMuon1->Phi;
                ZLepton1_Mass = 0.105658;
                ZLepton2_Pt   = ZMuon2->PT;
                ZLepton2_Eta  = ZMuon2->Eta;
                ZLepton2_Phi  = ZMuon2->Phi;
                ZLepton2_Mass = 0.105658;
                sameCharge = (ZMuon1->Charge == ZMuon2->Charge ? 1 : 0);
            }
        }
        else{
            ZLepton1_Pt   = 0.;
            ZLepton1_Eta  = 0.;
            ZLepton1_Phi  = 0.;
            ZLepton1_Mass = 0.;
            ZLepton2_Pt   = 0.;
            ZLepton2_Eta  = 0.;
            ZLepton2_Phi  = 0.;
            ZLepton2_Mass = 0.;
            sameCharge = -1;
        }

        ZReco_Pt   = ZReco.Pt();
        ZReco_Eta  = ZReco.Eta();
        ZReco_Phi  = ZReco.Phi();
        ZReco_Mass = ZReco.M();     

        // Fill output tree
        outtree->Fill();

    }

    // Write tree to file and save file
    outfile->Write();
    outfile->Save();

}
