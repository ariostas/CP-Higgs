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

// Declare functions
Double_t deltaR(TLorentzVector, TLorentzVector);

void select(const TString sample="", const TString tempinput="bacon_sample.root", const Double_t xsec = 1,
    const Int_t eosflag = 2){

    // read input input file
    TChain chain("Delphes");
    if(eosflag==0) chain.Add("/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples/" + tempinput);
    else chain.Add(tempinput);
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    // set up branches to read in from file
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchChargedHadron = treeReader->UseBranch("ChargedHadron");
    TClonesArray *branchNeutralHadron = treeReader->UseBranch("NeutralHadron");
    
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
    Jet *genJet=0;
    Track *chargedHadron=0;
    Tower *neutralHadron=0;
    
    // set up storage variables
    Jet *ZJet1=0, *ZJet2=0;
    Electron *ZElectron1=0, *ZElectron2=0;
    Muon *ZMuon1=0, *ZMuon2=0;
    Track *CPion1=0, *CPion2=0;
    //Tower *NPion1=0, *NPion2=0;
    Photon *NPion11=0, *NPion12=0, *NPion21=0, *NPion22=0;

    // Output file and trees
    TString output = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/out_" + tempinput;
    if(eosflag != 0) output = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small_autogen/out_" + sample + ".root";
    TFile *outfile = new TFile(output, "RECREATE");

    TTree *infoTree = new TTree("Count", "Count");
    infoTree->Branch("numberOfEntries", &numberOfEntries, "numberOfEntries/i");
    infoTree->Fill();

    TTree *outtree = new TTree("Events", "Events");

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

    for (Int_t iEntry=0; iEntry<chain.GetEntries(); iEntry++){
        
        treeReader->ReadEntry(iEntry);

        // Reset variables
        NLeptons=NJets=NTauJets=hasH=sameCharge=NCPions1=NCPions2=NNPions1=NNPions2=0;
        ZJet1=ZJet2=0;
        ZElectron1=ZElectron2=0;
        ZMuon1=ZMuon2=0;
        CPion1=CPion2=0;
        NPion11=NPion12=NPion21=NPion22=0;

        for (Int_t iElectron=0; iElectron<branchElectron->GetEntries(); iElectron++) {
            electron = (Electron*) branchElectron->At(iElectron);
            
            if(electron->P4().P() < 10. || electron->IsolationVar > 0.4) continue;
            
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
        //else if(ZElectron1->Charge == ZElectron2->Charge) ZElectron1 = ZElectron2 = 0;
        //else if(fabs((ZElectron1->P4()+ZElectron2->P4()).M() - 91.2) > 5.) ZElectron1 = ZElectron2 = 0;

        for (Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) {
            muon = (Muon*) branchMuon->At(iMuon);
            
            if(muon->P4().P() < 10. || muon->IsolationVar > 0.4) continue;
            
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
        //else if(ZMuon1->Charge == ZMuon2->Charge) ZMuon1 = ZMuon2 = 0;
        //else if(fabs((ZMuon1->P4()+ZMuon2->P4()).M() - 91.2) > 5.) ZMuon1 = ZMuon2 = 0;

        vector<Track*> orderedChargedHadrons;
        for(Int_t iChargedHadron=0; iChargedHadron<branchChargedHadron->GetEntries(); iChargedHadron++) {
            
            Track *tempTrack=0;
            for(Int_t iChargedHadron2=0; iChargedHadron2<branchChargedHadron->GetEntries(); iChargedHadron2++) {
                chargedHadron = (Track*) branchChargedHadron->At(iChargedHadron2);

                if(chargedHadron->P4().P() < 0.5) continue;

                bool isDifferent = true;
                for(Int_t x = 0; x < orderedChargedHadrons.size(); x++){
                    if(chargedHadron == orderedChargedHadrons.at(x)) isDifferent = false;
                }
                if(!isDifferent) continue;

                bool closeToLep = false;
                if(ZElectron1 && ZElectron2){
                    if(deltaR(ZElectron1->P4(), chargedHadron->P4()) < 0.4) closeToLep = true;
                    if(deltaR(ZElectron2->P4(), chargedHadron->P4()) < 0.4) closeToLep = true;
                }
                if(ZMuon1 && ZMuon2){
                    if(deltaR(ZMuon1->P4(), chargedHadron->P4()) < 0.4) closeToLep = true;
                    if(deltaR(ZMuon2->P4(), chargedHadron->P4()) < 0.4) closeToLep = true;
                }

                if(closeToLep) continue;
                
                if(!tempTrack) tempTrack = chargedHadron;
                else if(chargedHadron->P4().P() > tempTrack->P4().P()) tempTrack = chargedHadron;

            }

            if(tempTrack) orderedChargedHadrons.push_back(tempTrack);
            else break;

        }

        Double_t maxMass=0;
        for(Int_t iChargedHadron=0; iChargedHadron<orderedChargedHadrons.size(); iChargedHadron++){
            
            for(Int_t iChargedHadron2=iChargedHadron+1; iChargedHadron2<orderedChargedHadrons.size(); iChargedHadron2++) {
                
                Track *tempTrack1 = orderedChargedHadrons.at(iChargedHadron), *tempTrack2 = orderedChargedHadrons.at(iChargedHadron2);
                
                if(tempTrack1->Charge == tempTrack2->Charge) continue;
                if(deltaR(tempTrack1->P4(), tempTrack2->P4()) < 2.25) continue;
                if((tempTrack1->P4() + tempTrack2->P4()).M() < maxMass) continue;
                maxMass = (tempTrack1->P4() + tempTrack2->P4()).M();

                if(tempTrack1->Charge == 1){
                    CPion1 = tempTrack1;
                    CPion2 = tempTrack2;
                }
                else{
                    CPion1 = tempTrack2;
                    CPion2 = tempTrack1;
                }

            }

        }

        if(!CPion1 || !CPion2) continue;

        TLorentzVector V4CPions1, V4CPions2;
        for(Int_t iChargedHadron=0; iChargedHadron<branchChargedHadron->GetEntries(); iChargedHadron++) {
            Track *tempTrack = (Track*) branchChargedHadron->At(iChargedHadron);

            if(tempTrack->P4().P() < .5) continue;

            if(deltaR(tempTrack->P4(), CPion1->P4()) < 0.3){
                NCPions1++;
                V4CPions1+=tempTrack->P4();
            }
            if(deltaR(tempTrack->P4(), CPion2->P4()) < 0.3){
                NCPions2++;
                V4CPions2+=tempTrack->P4();
            }
        }

        if(V4CPions1.Pt() == 0.) cout << NCPions1 << " " << CPion1->P4().P() << endl;

        Int_t nPhotons1=0, nPhotons2=0;
        TLorentzVector V4NPions1, V4NPions2;
        for(Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) {
            photon = (Photon*) branchPhoton->At(iPhoton);

            if(photon->P4().P() < .5) continue;

            if(deltaR(photon->P4(), CPion1->P4()) < 0.4){
                nPhotons1++;
                V4NPions1+=photon->P4();
            }
            if(deltaR(photon->P4(), CPion2->P4()) < 0.4){
                nPhotons2++;
                V4NPions2+=photon->P4();
            }
        }

        if(nPhotons1 == 2 || nPhotons1 == 4 || nPhotons1 == 6) NNPions1 = nPhotons1/2;
        else NNPions1 = -nPhotons1;
        if(nPhotons2 == 2 || nPhotons2 == 4 || nPhotons2 == 6) NNPions2 = nPhotons2/2;
        else NNPions2 = -nPhotons2;

        for (Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) {
            photon = (Photon*) branchPhoton->At(iPhoton);
            
            if(photon->P4().P() < .5) continue;

            if(deltaR(photon->P4(), CPion1->P4()) < 0.4){
                if(!NPion11){
                    NPion11 = photon;
                }
                else if(photon->P4().P() > NPion11->P4().P()){
                    NPion12 = NPion11;
                    NPion11 = photon;
                }
                else if(!NPion12){
                    NPion12 = photon;
                }
                else if(photon->P4().P() > NPion12->P4().P()){
                    NPion12 = photon;
                }
            }

            else if(deltaR(photon->P4(), CPion2->P4()) < 0.4){
                if(!NPion21){
                    NPion21 = photon;
                }
                else if(photon->P4().P() > NPion21->P4().P()){
                    NPion22 = NPion21;
                    NPion21 = photon;
                }
                else if(!NPion22){
                    NPion22 = photon;
                }
                else if(photon->P4().P() > NPion22->P4().P()){
                    NPion22 = photon;
                }
            }

        }

        if(!NPion11 || !NPion12 || !NPion21 || !NPion22) continue;

        TLorentzVector vnpion1, vnpion2;
        vnpion1 = NPion11->P4() + NPion12->P4(); vnpion2 = NPion21->P4() + NPion22->P4();

        if(deltaR(vnpion1, CPion1->P4()) > 0.4 || deltaR(vnpion2, CPion2->P4()) > 0.4) continue;

        for(Int_t iParticle=0; iParticle<branchParticle->GetEntries(); iParticle++) {
            particle = (GenParticle*) branchParticle->At(iParticle);

            if(abs(particle->PID)==25) hasH = 1;
        }

        TLorentzVector recoZ; recoZ.SetPxPyPzE(0,0,0,0);

        for(Int_t iChargedHadron=0; iChargedHadron<branchChargedHadron->GetEntries(); iChargedHadron++) {
            chargedHadron = (Track*) branchChargedHadron->At(iChargedHadron);

            if(chargedHadron->P4().P() < .1) continue;
            if(deltaR(chargedHadron->P4(), CPion1->P4()) < 0.2) continue;
            if(deltaR(chargedHadron->P4(), CPion2->P4()) < 0.2) continue;
            if(deltaR(chargedHadron->P4(), vnpion1) < 0.2) continue;
            if(deltaR(chargedHadron->P4(), vnpion2) < 0.2) continue;

            recoZ += chargedHadron->P4();

        }

        for(Int_t iNeutralHadron=0; iNeutralHadron<branchNeutralHadron->GetEntries(); iNeutralHadron++) {
            neutralHadron = (Tower*) branchNeutralHadron->At(iNeutralHadron);

            if(neutralHadron->P4().P() < .1) continue;
            if(deltaR(neutralHadron->P4(), CPion1->P4()) < 0.2) continue;
            if(deltaR(neutralHadron->P4(), CPion2->P4()) < 0.2) continue;
            if(deltaR(neutralHadron->P4(), vnpion1) < 0.2) continue;
            if(deltaR(neutralHadron->P4(), vnpion2) < 0.2) continue;

            recoZ += neutralHadron->P4();

        }

        for(Int_t iPhoton=0; iPhoton<branchPhoton->GetEntries(); iPhoton++) {
            photon = (Photon*) branchPhoton->At(iPhoton);

            if(photon->P4().P() < .1) continue;
            if(deltaR(photon->P4(), CPion1->P4()) < 0.4) continue;
            if(deltaR(photon->P4(), CPion2->P4()) < 0.4) continue;
            //if(deltaR(photon->P4(), vnpion1) < 0.5) continue;
            //if(deltaR(photon->P4(), vnpion2) < 0.5) continue;

            recoZ += photon->P4();

        }

        for(Int_t iElectron=0; iElectron<branchElectron->GetEntries(); iElectron++) {
            electron = (Electron*) branchElectron->At(iElectron);

            if(electron->P4().P() < .1) continue;
            if(deltaR(electron->P4(), CPion1->P4()) < 0.2) continue;
            if(deltaR(electron->P4(), CPion2->P4()) < 0.2) continue;
            if(deltaR(electron->P4(), vnpion1) < 0.2) continue;
            if(deltaR(electron->P4(), vnpion2) < 0.2) continue;
            TLorentzVector tempEletron = electron->P4();
            tempEletron.SetPtEtaPhiM(tempEletron.Pt(), tempEletron.Eta(), tempEletron.Phi(), 0.000510999);

            recoZ += tempEletron;

        }

        for(Int_t iMuon=0; iMuon<branchMuon->GetEntries(); iMuon++) {
            muon = (Muon*) branchMuon->At(iMuon);

            if(muon->P4().P() < .1) continue;
            if(deltaR(muon->P4(), CPion1->P4()) < 0.2) continue;
            if(deltaR(muon->P4(), CPion2->P4()) < 0.2) continue;
            if(deltaR(muon->P4(), vnpion1) < 0.2) continue;
            if(deltaR(muon->P4(), vnpion2) < 0.2) continue;
            TLorentzVector tempMuon = muon->P4();
            tempMuon.SetPtEtaPhiM(tempMuon.Pt(), tempMuon.Eta(), tempMuon.Phi(), 0.105658);

            recoZ += tempMuon;

        }

        Double_t minDeltaM = 999;
        for (Int_t iJet=0; iJet<branchJet->GetEntries(); iJet++) {
            jet = (Jet*) branchJet->At(iJet);
            
            if(jet->P4().P() < 30.) continue;
            
            NJets++;
            if(jet->TauTag == 1) NTauJets++;

            for(Int_t iJet2=iJet+1; iJet2<branchJet->GetEntries(); iJet2++){
                Jet *jet2 = (Jet*) branchJet->At(iJet2);

                if(jet2->P4().P() < 30.) continue;

                TLorentzVector vTempZ = jet->P4() + jet2->P4();

                Double_t tempDeltaM = fabs(vTempZ.M() - 91.2);

                if(tempDeltaM < minDeltaM){
                    minDeltaM = tempDeltaM;
                    ZJet1 = jet;
                    ZJet2 = jet2;
                }
            }
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
        CPions1_Pt   = V4CPions1.Pt();
        CPions1_Eta  = V4CPions1.Eta();
        CPions1_Phi  = V4CPions1.Phi();
        CPions1_Mass = V4CPions1.M();
        CPions2_Pt   = V4CPions2.Pt();
        CPions2_Eta  = V4CPions2.Eta();
        CPions2_Phi  = V4CPions2.Phi();
        CPions2_Mass = V4CPions2.M();

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
        NPions1_Pt   = V4NPions1.Pt();
        NPions1_Eta  = V4NPions1.Eta();
        NPions1_Phi  = V4NPions1.Phi();
        NPions1_Mass = V4NPions1.M();
        NPions2_Pt   = V4NPions2.Pt();
        NPions2_Eta  = V4NPions2.Eta();
        NPions2_Phi  = V4NPions2.Phi();
        NPions2_Mass = V4NPions2.M();

        if(ZJet1 && ZJet2){
            ZJet1_Pt   = ZJet1->PT;
            ZJet1_Eta  = ZJet1->Eta;
            ZJet1_Phi  = ZJet1->Phi;
            ZJet1_Mass = ZJet1->P4().M();
            ZJet2_Pt   = ZJet2->PT;
            ZJet2_Eta  = ZJet2->Eta;
            ZJet2_Phi  = ZJet2->Phi;
            ZJet2_Mass = ZJet2->P4().M();
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

Double_t deltaR(TLorentzVector vector1, TLorentzVector vector2){

  const Double_t pi = 3.14159265358979;
  const Double_t eta1 = vector1.Eta(), eta2 = vector2.Eta();
  const Double_t phi1 = vector1.Phi(), phi2 = vector2.Phi();

  Double_t etaDiff = (eta1-eta2);
  Double_t phiDiff = fabs(phi1-phi2);
  while (phiDiff>pi) phiDiff = fabs(phiDiff-2.0*pi);

  Double_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

  return TMath::Sqrt(deltaRSquared);

}