#if !defined(__CINT__) || defined(__MAKECINT__)
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
#include <TString.h>

#include <TLorentzVector.h>
#include <TVector3.h>

#include <TRandom3.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include "modules/Delphes.h"                   // delphes
#include "ExRootAnalysis/ExRootTreeReader.h"   // delphes
#include "classes/DelphesClasses.h"            // delphes

#endif

using namespace TMath;
using namespace std;

void selectAlt(const TString tempinput="zh_delphes_0pi12_1.root", const Double_t xsec = 1,
    const Int_t eosflag = 1){

    // read input input file
    TChain chain("Delphes");
    if(eosflag==0) chain.Add("/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples/" + tempinput);
    else if(eosflag==1) chain.Add("/afs/cern.ch/work/a/ariostas/private/CP-Higgs_Samples_temp/" + tempinput);
    else if(eosflag==2) chain.Add("root://eoscms//store/user/arapyan/mc/Delphes/" + tempinput);
    else chain.Add(tempinput);
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    // set up branches to read in from file
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchChargedHadron = treeReader->UseBranch("Track");
    TClonesArray *branchNeutralHadron = treeReader->UseBranch("Tower");
    
    if (!(branchJet)) {
        cout << "  file broken" << endl;
        return;
    }

    // set up loop variables
    Jet *jet=0;
    Electron *electron=0;
    Photon *photon=0;
    Muon *muon=0;
    GenParticle *particle=0;
    Track *chargedHadron=0;
    Tower *neutralHadron=0;
    
    // set up storage variables
    Jet *ZJet1=0, *ZJet2=0;
    Jet *JetTau1=0, *JetTau2=0;
    Electron *ZElectron1=0, *ZElectron2=0;
    Muon *ZMuon1=0, *ZMuon2=0;
    Track *CPion1=0, *CPion2=0;

    // Output file and trees
    TString output = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_smallAlt/out_" + tempinput;
    TFile *outfile = new TFile(output, "RECREATE");

    TTree *infoTree = new TTree("Count", "Count");
    infoTree->Branch("numberOfEntries", &numberOfEntries, "numberOfEntries/i");
    infoTree->Fill();

    TTree *outtree = new TTree("Events", "Events");

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

    for (Int_t iEntry=0; iEntry<chain.GetEntries(); iEntry++){
        
        treeReader->ReadEntry(iEntry);

        // Reset variables
        NLeptons=NJets=NTauJets=hasH=ZFromLep=oppositeTrackCharge=0;
        ZJet1=ZJet2=JetTau1=JetTau2=0;
        ZElectron1=ZElectron2=0;
        ZMuon1=ZMuon2=0;
        CPion1=CPion2=0;


        for(Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) {
            particle = (GenParticle*) branchParticle->At(iParticle);

            if(abs(particle->PID)==25) hasH = 1;
        }

        vector<Jet*> TauCandidates;
        vector<Track*> CPionCandidates;
        vector<TLorentzVector> NPionCandidates, TrackSums, PhotonSums;

        for(Int_t x=0; x<branchJet->GetEntries(); x++){
            jet = (Jet*) branchJet->At(x);
            
            if(jet->PT < 20) continue;

            NJets++;

            if(jet->TauTag) NTauJets++;

            Track *pi=0;

            for(Int_t iChargedHadron=0; iChargedHadron<branchChargedHadron->GetEntries(); iChargedHadron++) {
                chargedHadron = (Track*) branchChargedHadron->At(iChargedHadron);
                
                if(chargedHadron->PT < 0.5) continue;
                if(chargedHadron->P4().DeltaR(jet->P4()) > 0.2) continue;
                
                if(!pi) pi = chargedHadron;
                else if(chargedHadron->P4().P() > pi->P4().P()) pi = chargedHadron;
            }

            if(!pi) continue;

            TLorentzVector TracksP4;
            for(Int_t iChargedHadron=0; iChargedHadron<branchChargedHadron->GetEntries(); iChargedHadron++){
                chargedHadron = (Track*) branchChargedHadron->At(iChargedHadron);
                
                if(chargedHadron->P4().DeltaR(jet->P4()) > 0.4) continue;
                if(chargedHadron == pi) continue;

                TracksP4 += chargedHadron->P4();
            }

            if(TracksP4.P() > 5.) continue;

            Photon *phot1=0, *phot2=0;
            for(Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) {
                photon = (Photon*) branchPhoton->At(iPhoton);
                if(photon->P4().DeltaR(jet->P4()) > 0.4 || photon->PT < 0.5 || Abs(photon->Eta) > 2.5) continue;
                
                if(!phot1){
                    phot1 = photon;
                }
                else if(photon->P4().P() > phot1->P4().P()){
                    phot2 = phot1;
                    phot1 = photon;
                }
                else if(!phot2){
                    phot2 = photon;
                }
                else if(photon->P4().P() > phot2->P4().P()){
                    phot2 = photon;
                }
            }

            if(!phot1 || !phot2) continue;

            TLorentzVector pi0 = phot1->P4() + phot2->P4();

            TLorentzVector PhotonsP4;
            for(Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) {
                photon = (Photon*) branchPhoton->At(iPhoton);
                // cout << photon->P4().P() << endl;
                if(photon->P4().DeltaR(jet->P4()) > 0.4) continue;
                if(photon == phot1 || photon == phot2) continue;

                PhotonsP4 += photon->P4();
            }
            if(PhotonsP4.P() > 5.) continue;

            TauCandidates.push_back(jet);
            CPionCandidates.push_back(pi);
            NPionCandidates.push_back(pi0);
            TrackSums.push_back(TracksP4);
            PhotonSums.push_back(PhotonsP4);

        }

        if(TauCandidates.size() < 2) continue;

        TLorentzVector NPion1, NPion2, TracksTau1, TracksTau2, PhotonsTau1, PhotonsTau2;

        for(Int_t n = 0; n < TauCandidates.size(); n++){

            if(!CPion1){
                CPion1 = CPionCandidates.at(n);
                NPion1 = NPionCandidates.at(n);
                JetTau1 = TauCandidates.at(n);
                TracksTau1 = TrackSums.at(n);
                PhotonsTau1 = PhotonSums.at(n);
            }
            else if(TauCandidates.at(n)->P4().P() > CPion1->P4().P()){
                CPion2 = CPion1;
                NPion2 = NPion1;
                JetTau2 = JetTau1;
                TracksTau2 = TracksTau1;
                PhotonsTau2 = PhotonsTau1;
                CPion1 = CPionCandidates.at(n);
                NPion1 = NPionCandidates.at(n);
                JetTau1 = TauCandidates.at(n);
                TracksTau1 = TrackSums.at(n);
                PhotonsTau1 = PhotonSums.at(n);
            }
            else if(!CPion2){
                CPion2 = CPionCandidates.at(n);
                NPion2 = NPionCandidates.at(n);
                JetTau2 = TauCandidates.at(n);
                TracksTau2 = TrackSums.at(n);
                PhotonsTau2 = PhotonSums.at(n);
            }
            else if(TauCandidates.at(n)->P4().P() > CPion2->P4().P()){
                CPion2 = CPionCandidates.at(n);
                NPion2 = NPionCandidates.at(n);
                JetTau2 = TauCandidates.at(n);
                TracksTau2 = TrackSums.at(n);
                PhotonsTau2 = PhotonSums.at(n);
            }
        }

        TauCandidates.clear();
        CPionCandidates.clear();
        NPionCandidates.clear();
        TrackSums.clear();
        PhotonSums.clear();

        Track *tempT1 = CPion1, *tempT2 = CPion2;
        TLorentzVector tempV1 = NPion1, tempV2 = NPion2;
        Jet *tempJ1 = JetTau1, *tempJ2 = JetTau2;
        TLorentzVector tempTs1 = TracksTau1, tempTs2 = TracksTau2;
        TLorentzVector tempP1 = PhotonsTau1, tempP2 = PhotonsTau2;
        if(tempT1->Charge > 0){
            CPion1 = tempT1;
            NPion1 = tempV1;
            JetTau1 = tempJ1;
            TracksTau1 = tempTs1;
            PhotonsTau1 = tempP1;
            CPion2 = tempT2;
            NPion2 = tempV2;
            JetTau2 = tempJ2;
            TracksTau2 = tempTs2;
            PhotonsTau2 = tempP2;
        }
        else{
            CPion1 = tempT2;
            NPion1 = tempV2;
            JetTau1 = tempJ2;
            TracksTau1 = tempTs2;
            PhotonsTau1 = tempP2;
            CPion2 = tempT1;
            NPion2 = tempV1;
            JetTau2 = tempJ1;
            TracksTau2 = tempTs1;
            PhotonsTau2 = tempP1;
        }

        if(CPion1->Charge == CPion2->Charge) oppositeTrackCharge = 0;
        else oppositeTrackCharge = 1;

        for(Int_t iElectron=0; iElectron<branchElectron->GetEntries(); iElectron++){
            electron = (Electron*) branchElectron->At(iElectron);
            
            if(electron->PT < 10. || electron->IsolationVar > 0.4) continue;
            
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

        if(!ZElectron1 || !ZElectron2) ZElectron1 = ZElectron2 = 0;
        else if(ZElectron1->Charge == ZElectron2->Charge) continue;

        for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) {
            muon = (Muon*) branchMuon->At(iMuon);
            
            if(muon->PT < 10. || muon->IsolationVar > 0.4) continue;
            
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

        if(!ZMuon1 || !ZMuon2) ZMuon1 = ZMuon2 = 0;
        else if(ZMuon1->Charge == ZMuon2->Charge) continue;

        for(Int_t x=0; x<branchJet->GetEntries(); x++){
            jet = (Jet*) branchJet->At(x);
            
            if(jet->PT < 20) continue;

            if(jet->P4().DeltaR(JetTau1->P4()) < 0.8 || jet->P4().DeltaR(JetTau1->P4()) < 0.8) continue;

            if(!ZJet1){
                ZJet1 = jet;
            }
            else if(jet->P4().P() > ZJet1->P4().P()){
                ZJet2 = ZJet1;
                ZJet1 = jet;
            }
            else if(ZJet2){
                ZJet2 = jet;
            }
            else if(jet->P4().P() > ZJet1->P4().P()){
                ZJet2 = jet;
            }

        }

        Double_t minDR = 0.2;
        TLorentzVector recoZ; recoZ.SetPxPyPzE(0,0,0,0);

        for(Int_t iChargedHadron=0; iChargedHadron<branchChargedHadron->GetEntries(); iChargedHadron++) {
            chargedHadron = (Track*) branchChargedHadron->At(iChargedHadron);

            if(chargedHadron->PT < .5) continue;
            if(Abs(chargedHadron->Eta) > 2.5) continue;
            if(chargedHadron->P4().DeltaR(JetTau1->P4()) < minDR) continue;
            if(chargedHadron->P4().DeltaR(JetTau2->P4()) < minDR) continue;

            recoZ += chargedHadron->P4();

        }

        for(Int_t iNeutralHadron=0; iNeutralHadron<branchNeutralHadron->GetEntries(); iNeutralHadron++) {
            neutralHadron = (Tower*) branchNeutralHadron->At(iNeutralHadron);

            if(neutralHadron->P4().Pt() < .5) continue;
            if(Abs(neutralHadron->Eta) > 2.5) continue;
            if(neutralHadron->P4().DeltaR(JetTau1->P4()) < minDR) continue;
            if(neutralHadron->P4().DeltaR(JetTau2->P4()) < minDR) continue;

            recoZ += neutralHadron->P4();

        }

        for(Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) {
            photon = (Photon*) branchPhoton->At(iPhoton);

            if(photon->PT < .5) continue;
            if(Abs(photon->Eta) > 2.5) continue;
            if(photon->P4().DeltaR(JetTau1->P4()) < minDR+0.2) continue;
            if(photon->P4().DeltaR(JetTau2->P4()) < minDR+0.2) continue;

            recoZ += photon->P4();

        }

        for(Int_t iElectron=0; iElectron<branchElectron->GetEntries(); iElectron++) {
            electron = (Electron*) branchElectron->At(iElectron);

            if(electron->PT < .5) continue;
            if(chargedHadron->P4().DeltaR(JetTau1->P4()) < minDR) continue;
            if(chargedHadron->P4().DeltaR(JetTau2->P4()) < minDR) continue;
            TLorentzVector tempEletron = electron->P4();
            tempEletron.SetPtEtaPhiM(tempEletron.Pt(), tempEletron.Eta(), tempEletron.Phi(), 0.000510999);

            recoZ += tempEletron;

        }

        for(Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) {
            muon = (Muon*) branchMuon->At(iMuon);

            if(muon->PT < .5) continue;
            if(chargedHadron->P4().DeltaR(JetTau1->P4()) < minDR) continue;
            if(chargedHadron->P4().DeltaR(JetTau2->P4()) < minDR) continue;
            TLorentzVector tempMuon = muon->P4();
            tempMuon.SetPtEtaPhiM(tempMuon.Pt(), tempMuon.Eta(), tempMuon.Phi(), 0.105658);

            recoZ += tempMuon;

        }

        eventWeight = 1000000.*xsec;

        CPion1_Pt   = CPion1->PT;
        CPion1_Eta  = CPion1->Eta;
        CPion1_Phi  = CPion1->Phi;
        CPion1_Mass = 0.139570;
        CPion2_Pt   = CPion2->PT;
        CPion2_Eta  = CPion2->Eta;
        CPion2_Phi  = CPion2->Phi;
        CPion2_Mass = 0.139570;

        NPion1_Pt   = NPion1.Pt();
        NPion1_Eta  = NPion1.Eta();
        NPion1_Phi  = NPion1.Phi();
        NPion1_Mass = NPion1.M();
        NPion2_Pt   = NPion2.Pt();
        NPion2_Eta  = NPion2.Eta();
        NPion2_Phi  = NPion2.Phi();
        NPion2_Mass = NPion2.M();

        JetTau1_Pt   = JetTau1->PT;
        JetTau1_Eta  = JetTau1->Eta;
        JetTau1_Phi  = JetTau1->Phi;
        JetTau1_Mass = JetTau1->Mass;
        JetTau2_Pt   = JetTau2->PT;
        JetTau2_Eta  = JetTau2->Eta;
        JetTau2_Phi  = JetTau2->Phi;
        JetTau2_Mass = JetTau2->Mass;

        TracksTau1_Pt   = TracksTau1.Pt();
        TracksTau1_Eta  = TracksTau1.Eta();
        TracksTau1_Phi  = TracksTau1.Phi();
        TracksTau1_Mass = TracksTau1.M();
        TracksTau2_Pt   = TracksTau2.Pt();
        TracksTau2_Eta  = TracksTau2.Eta();
        TracksTau2_Phi  = TracksTau2.Phi();
        TracksTau2_Mass = TracksTau2.M();

        PhotonsTau1_Pt   = PhotonsTau1.Pt();
        PhotonsTau1_Eta  = PhotonsTau1.Eta();
        PhotonsTau1_Phi  = PhotonsTau1.Phi();
        PhotonsTau1_Mass = PhotonsTau1.M();
        PhotonsTau2_Pt   = PhotonsTau2.Pt();
        PhotonsTau2_Eta  = PhotonsTau2.Eta();
        PhotonsTau2_Phi  = PhotonsTau2.Phi();
        PhotonsTau2_Mass = PhotonsTau2.M();
        
        if(ZElectron1 && ZElectron2){
            ZFromLep = 1;
            ZParticle1_Pt   = ZElectron1->PT;
            ZParticle1_Eta  = ZElectron1->Eta;
            ZParticle1_Phi  = ZElectron1->Phi;
            ZParticle1_Mass = 0.000510999;
            ZParticle2_Pt   = ZElectron2->PT;
            ZParticle2_Eta  = ZElectron2->Eta;
            ZParticle2_Phi  = ZElectron2->Phi;
            ZParticle2_Mass = 0.000510999;
        }
        else if(ZMuon1 && ZMuon2){
            ZFromLep = 2;
            ZParticle1_Pt   = ZMuon1->PT;
            ZParticle1_Eta  = ZMuon1->Eta;
            ZParticle1_Phi  = ZMuon1->Phi;
            ZParticle1_Mass = 0.105658;
            ZParticle2_Pt   = ZMuon2->PT;
            ZParticle2_Eta  = ZMuon2->Eta;
            ZParticle2_Phi  = ZMuon2->Phi;
            ZParticle2_Mass = 0.105658;
        }
        else if(ZJet1 && ZJet2){
            ZFromLep = 0;
            ZParticle1_Pt   = ZJet1->PT;
            ZParticle1_Eta  = ZJet1->Eta;
            ZParticle1_Phi  = ZJet1->Phi;
            ZParticle1_Mass = ZJet1->Mass;
            ZParticle2_Pt   = ZJet2->PT;
            ZParticle2_Eta  = ZJet2->Eta;
            ZParticle2_Phi  = ZJet2->Phi;
            ZParticle2_Mass = ZJet2->Mass;
        }
        else{
            ZFromLep = -1;
            ZParticle1_Pt   = 0;
            ZParticle1_Eta  = 0;
            ZParticle1_Phi  = 0;
            ZParticle1_Mass = 0;
            ZParticle2_Pt   = 0;
            ZParticle2_Eta  = 0;
            ZParticle2_Phi  = 0;
            ZParticle2_Mass = 0;
        }

        if(recoZ.Pt() > 0.){
            ZReco_Pt   = recoZ.Pt();
            ZReco_Eta  = recoZ.Eta();
            ZReco_Phi  = recoZ.Phi();
            ZReco_Mass = recoZ.M();
        }
        else{
            ZReco_Pt   = 0;
            ZReco_Eta  = 0;
            ZReco_Phi  = 0;
            ZReco_Mass = 0;
        }

        outtree->Fill();

    }

    outfile->Write();
    outfile->Save();

}
