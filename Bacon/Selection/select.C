#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // Access to gROOT, entry point to ROOT system
#include <TSystem.h>                // Interface to OS
#include <TFile.h>                  // File handle class
#include <TTree.h>                  // Class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // Class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // Standard I/O
#include <sstream>                  // Standard I/O
#include <iomanip>                  // Functions to format standard I/O
#include <fstream>                  // Functions for file I/O
#include <TChain.h>
#include <TH1.h>
#include "TLorentzVector.h"         // 4-vector class

#include "TGenEventInfo.hh"         //
#include "TGenJet.hh"               // Bacon classes
#include "TGenParticle.hh"          //

#endif

Double_t deltaR( const Double_t eta1, const Double_t eta2, const Double_t phi1, const Double_t phi2 );

using namespace baconhep;

void select(const TString sample="", const TString tempinput="bacon_sample.root", const Double_t xsec = 1,
    const Int_t eosflag = 2)
{
    TString input;
    if(eosflag==2) input = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs/" + tempinput;
    else if(eosflag==1) input = "root://eoscms.cern.ch//store/user/arapyan/higgs_cp_samples/" + sample + "/Bacon/" + tempinput;
    else input = "/afs/cern.ch/work/a/arapyan/public/forMarkus/higgs_cp_samples/" + sample + "/" + tempinput;

    TString output = "/afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/" + sample + "/" + tempinput;

    TChain chain("Events");
    chain.Add(input);

    const Double_t DR=0.4;

    // Data structures to store info from TTrees
    TGenEventInfo *info = new TGenEventInfo();
    TClonesArray *jet = new TClonesArray("baconhep::TGenJet");
    TClonesArray *part = new TClonesArray("baconhep::TGenParticle");

    chain.SetBranchAddress("GenEvtInfo",  &info);
    TBranch *infoBr = chain.GetBranch("GenEvtInfo");
    chain.SetBranchAddress("GenJet" ,     &jet );
    TBranch *jetBr = chain.GetBranch("GenJet");
    chain.SetBranchAddress("GenParticle", &part);
    TBranch *partBr = chain.GetBranch("GenParticle");

    // Output file and trees
    TFile *outfile = new TFile(output, "RECREATE");

    TTree *infoTree = new TTree("Count", "Count");
    Long64_t n = chain.GetEntries();
    infoTree->Branch("n", &n, "n/i");
    infoTree->Fill();

    TTree *outtree = new TTree("Events", "Events");

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

    // Reco Z
    Double_t recoz_pt, recoz_eta, recoz_phi, recoz_mass;

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

    outtree->Branch("recoz_pt",         &recoz_pt,          "recoz_pt/D");
    outtree->Branch("recoz_eta",        &recoz_eta,         "recoz_eta/D");
    outtree->Branch("recoz_phi",        &recoz_phi,         "recoz_phi/D");
    outtree->Branch("recoz_mass",       &recoz_mass,        "recoz_mass/D");

    for (Int_t i=0; i<chain.GetEntries(); i++)
    {
        if(i==0 || (i+1)%1000==0) cout << "Processing event " << i+1 << endl;

        infoBr->GetEntry(i);

        part->Clear();
        partBr->GetEntry(i);
        jet->Clear();
        //jetBr->GetEntry(i);

        if (part->GetEntries()==0) continue;

        Int_t iTau1=-1, iTau2=-1;
        Int_t iT1=-1, iT2=-1;

        Int_t iZ1=-1, iZ2=-1;
        Int_t ijetZ1=-1, ijetZ2=-1;

        Int_t iZ=-1;

        // Reset variables
        nLeptons=ncpions1=ncpions2=nnpions1=nnpions2=zToLep=0;

        jetTau1_pt=jetTau1_eta=jetTau1_phi=jetTau1_mass=0;
        jetTau2_pt=jetTau2_eta=jetTau2_phi=jetTau2_mass=0;

        genTau1_pt=genTau1_eta=genTau1_phi=genTau1_mass=0;
        genTau2_pt=genTau2_eta=genTau2_phi=genTau2_mass=0;

        visTau1_pt=visTau1_eta=visTau1_phi=visTau1_mass=0;
        visTau2_pt=visTau2_eta=visTau2_phi=visTau2_mass=0;

        cpions1_pt=cpions1_eta=cpions1_phi=cpions1_mass=0;
        cpions2_pt=cpions2_eta=cpions2_phi=cpions2_mass=0;

        npions1_pt=npions1_eta=npions1_phi=npions1_mass=0;
        npions2_pt=npions2_eta=npions2_phi=npions2_mass=0;

        z1_pt=z1_eta=z1_phi=z1_mass=0;
        z2_pt=z2_eta=z2_phi=z2_mass=0;

        z_pt=z_eta=z_phi=z_mass=0;
        recoz_pt=recoz_eta=recoz_phi=recoz_mass=0;

        eventWeight=xsec*1000.0;

        bool hasH=false, hasZ=false;

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            if(abs(genloop->pdgId)==25) hasH = true;
            if(abs(genloop->pdgId)==23) hasZ = true;
        }

        Int_t tausParent=-1;
        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            Int_t parentPdg=abs(dynamic_cast<TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId);

            if(hasH){
                if (abs(genloop->pdgId)==15)
                {
                    if (iTau1==-1 && parentPdg==25)
                    {
                        iTau1=j;
                        tausParent=genloop->parent;
                    }
                    else if (genloop->parent==iTau1)
                    {
                        iTau1=j;
                    }
                    else if (iTau2==-1 && parentPdg==25 && tausParent==genloop->parent)
                    {
                        iTau2=j;
                    }
                    else if (genloop->parent==iTau2)
                    {
                        iTau2=j;
                    }
                }
            }

            else if(hasZ){
                if (abs(genloop->pdgId)==15)
                {
                    if (iTau1==-1 && parentPdg==23)
                    {
                        iTau1=j;
                        tausParent=genloop->parent;
                    }
                    else if (genloop->parent==iTau1)
                    {
                        iTau1=j;
                    }
                    else if (iTau2==-1 && parentPdg==23 && tausParent==genloop->parent)
                    {
                        iTau2=j;
                    }
                    else if (genloop->parent==iTau2)
                    {
                        iTau2=j;
                    }
                }
            }

            else{
                if (abs(genloop->pdgId)==15)
                {
                    if (iTau1==-1)
                    {
                        iTau1=j;
                        tausParent=genloop->parent;
                    }
                    else if (genloop->parent==iTau1)
                    {
                        iTau1=j;
                    }
                    else if (iTau2==-1 && tausParent==genloop->parent)
                    {
                        iTau2=j;
                    }
                    else if (genloop->parent==iTau2)
                    {
                        iTau2=j;
                    }
                }
            }
        }

        if (iTau1==-1 || iTau2==-1) continue;

        Int_t isLep=0;
        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            Int_t parentPdg=abs(dynamic_cast<TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId);

            if (genloop->parent != iTau1 && genloop->parent != iTau2) continue;
            if ((abs(genloop->pdgId)==13||abs(genloop->pdgId)==11)){
                isLep=1;
            }
        }

        if (isLep!=0) continue;

        if(dynamic_cast<TGenParticle *>(part->At(iTau1))->pdgId==15 &&  dynamic_cast<TGenParticle *>(part->At(iTau2))->pdgId==-15)
        {
            Int_t temp=iTau1;
            iTau1=iTau2;
            iTau2=temp;
        }
        else if(dynamic_cast<TGenParticle *>(part->At(iTau1))->pdgId==-15 &&  dynamic_cast<TGenParticle *>(part->At(iTau2))->pdgId==15)
        {

        }
        else
        {
            cout << "Ignoring event." << endl;
            continue;
        }

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            TGenParticle* genloop = (TGenParticle*) ((*part)[j]);

            Int_t pdg=abs(genloop->pdgId);
            Int_t parent=genloop->parent;

            if(pdg==23){
                if(iZ==-1){
                    iZ=j;
                }
                else if(parent==iZ){
                    iZ=j;
                }
            }
        }

        TLorentzVector cpions1, cpions2, npions1, npions2;
        cpions1.SetPxPyPzE(0,0,0,0); cpions2.SetPxPyPzE(0,0,0,0);
        npions1.SetPxPyPzE(0,0,0,0); npions2.SetPxPyPzE(0,0,0,0);

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            TGenParticle* genloop = (TGenParticle*) ((*part)[j]);

            Int_t pdg=abs(genloop->pdgId);
            Int_t parent=genloop->parent;

            if(pdg!=211 && pdg!=111) continue;

            Int_t parentTau=0;

            if(parent==iTau1) parentTau=1;
            else if(parent==iTau2) parentTau=2;

            if(parentTau==0) continue;

            TLorentzVector vTemp;
            vTemp.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);

            if(parentTau==1)
            {
                if(pdg==211){
                    cpions1+=vTemp;
                    ncpions1++;
                }
                else{
                    npions1+=vTemp;
                    nnpions1++;
                }
            }
            else
            {
                if(pdg==211){
                    cpions2+=vTemp;
                    ncpions2++;
                }
                else{
                    npions2+=vTemp;
                    nnpions2++;
                }
            }

        }


        TLorentzVector vRecoZ;
        vRecoZ.SetPxPyPzE(0,0,0,0);
        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            if (fabs(genloop->pdgId)==12 || fabs(genloop->pdgId)==14 || fabs(genloop->pdgId)==14) continue;
            if(fabs(genloop->eta) > 5. || genloop->status!=1) continue;
            if(deltaR(genloop->eta, cpions1.Eta(), genloop->phi, cpions1.Phi()) < DR) continue;
            if(deltaR(genloop->eta, cpions2.Eta(), genloop->phi, cpions2.Phi()) < DR) continue;
            if(deltaR(genloop->eta, npions1.Eta(), genloop->phi, npions1.Phi()) < DR) continue;
            if(deltaR(genloop->eta, npions2.Eta(), genloop->phi, npions2.Phi()) < DR) continue;
            TLorentzVector vTemp;
            vTemp.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);
            vRecoZ+=vTemp;
        }

        const TGenParticle* genTau1 = (TGenParticle*) ((*part)[iTau1]);
        TLorentzVector vGenTau1;
        vGenTau1.SetPtEtaPhiM(genTau1->pt, genTau1->eta, genTau1->phi, genTau1->mass);

        const TGenParticle* genTau2 = (TGenParticle*) ((*part)[iTau2]);
        TLorentzVector vGenTau2;
        vGenTau2.SetPtEtaPhiM(genTau2->pt, genTau2->eta, genTau2->phi, genTau2->mass);

        TLorentzVector vVisTau1, vVisTau2;

        if (deltaR(genTau1->eta, genTau2->eta, genTau1->phi, genTau2->phi)<DR) continue;

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            if (fabs(genloop->pdgId)!=16) continue;
            TLorentzVector vTemp;
            vTemp.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, 0);
            if (genloop->parent==iTau1)
            {
                vVisTau1=vGenTau1-vTemp;
            }
            if (genloop->parent==iTau2)
            {
                vVisTau2=vGenTau2-vTemp;
            }
        }


        for (Int_t j=0; j<jet->GetEntries(); j++)
        {
            const TGenJet* loop = (TGenJet*) ((*jet)[j]);

            if(loop->pt < 10) continue;

            if (iT1==-1)
            {
                if(deltaR(loop->eta,vVisTau1.Eta(),loop->phi,vVisTau1.Phi())<DR){
                    iT1=j;
                    continue;
                }
            }
            if(iT1!=-1 && ((TGenJet*) ((*jet)[iT1]))->pt < loop->pt){
                if(deltaR(loop->eta,vVisTau1.Eta(),loop->phi,vVisTau1.Phi())<DR){
                    iT1=j;
                    continue;
                }
            }
            if (iT2==-1)
            {
                if(deltaR(loop->eta,vVisTau2.Eta(),loop->phi,vVisTau2.Phi())<DR){
                    iT2=j;
                    continue;
                }
            }
            if(iT2!=-1 && ((TGenJet*) ((*jet)[iT2]))->pt < loop->pt){
                if(deltaR(loop->eta,vVisTau2.Eta(),loop->phi,vVisTau2.Phi())<DR){
                  iT2=j;
                  continue;
                }
            }

        }

        for (Int_t j=0; j<jet->GetEntries(); j++)
        {
            const TGenJet* loop = (TGenJet*) ((*jet)[j]);

            if(loop->pt < 20) continue;

            if(deltaR(vVisTau1.Eta(), loop->eta, vVisTau1.Phi(), loop->phi)<DR+0.5) continue;
            if(deltaR(vVisTau2.Eta(), loop->eta, vVisTau2.Phi(), loop->phi)<DR+0.5) continue;
            
            if (ijetZ1==-1)
            {
                ijetZ1=j;
            }
            else if(((TGenJet*) ((*jet)[ijetZ1]))->pt < loop->pt){
                ijetZ2=ijetZ1;
                ijetZ1=j;
            }
            else if (ijetZ2==-1)
            {
                ijetZ2=j;
            }
            else if(((TGenJet*) ((*jet)[ijetZ2]))->pt < loop->pt){
                ijetZ2=j;
            }
        }

        if(ijetZ1==-1 || ijetZ2==-1) zToLep=1;

        if(zToLep==1){
            for (Int_t j=0; j<part->GetEntries(); j++)
            {
                TGenParticle *genloop = (TGenParticle*) ((*part)[j]);

                Int_t pdg=abs(genloop->pdgId);
                Int_t parent=genloop->parent;

                if(pdg!=11 && pdg!= 13) continue;
                if(genloop->pt < 20) continue;

                if(deltaR(vVisTau1.Eta(), genloop->eta, vVisTau1.Phi(), genloop->phi)<DR+0.5) continue;
                if(deltaR(vVisTau2.Eta(), genloop->eta, vVisTau2.Phi(), genloop->phi)<DR+0.5) continue;

                if (iZ1==-1)
                {
                    iZ1=j;
                }
                else if(((TGenParticle*) ((*part)[iZ1]))->pt < genloop->pt){
                    iZ1=j;
                }
            }

            if(iZ1!=-1){
                for (Int_t j=0; j<part->GetEntries(); j++)
                {
                    TGenParticle* genloop = (TGenParticle*) ((*part)[j]);

                    Int_t pdg=abs(genloop->pdgId);
                    Int_t parent=genloop->parent;

                    if(pdg!=11 && pdg!= 13) continue;
                    if(genloop->pt < 30) continue;
                    if(iZ1==j) continue;
                    if(pdg!=abs(((TGenParticle*) ((*part)[iZ1]))->pdgId)) continue;
                    if(deltaR(vVisTau1.Eta(), genloop->eta, vVisTau1.Phi(), genloop->phi)<DR+0.5) continue;
                    if(deltaR(vVisTau2.Eta(), genloop->eta, vVisTau2.Phi(), genloop->phi)<DR+0.5) continue;

                    if (iZ2==-1)
                    {
                        iZ2=j;
                    }
                    else if(((TGenParticle*) ((*part)[iZ2]))->pt < genloop->pt){
                        iZ2=j;
                    }
                }
            }
        }
        
        if(zToLep==1 && (iZ1==-1 || iZ2==-1)) zToLep=-1;

        if(iZ==-1 && (iZ1==-1 || iZ2==-1) && (ijetZ1==-1 || ijetZ2==-1)) continue;

        if(zToLep==0){
            const TGenJet *jetZ1, *jetZ2;

            jetZ1 = (TGenJet*) ((*jet)[ijetZ1]);
            jetZ2 = (TGenJet*) ((*jet)[ijetZ2]);

            z1_pt=jetZ1->pt;
            z1_eta=jetZ1->eta;
            z1_phi=jetZ1->phi;
            z1_mass=jetZ1->mass;

            z2_pt=jetZ2->pt;
            z2_eta=jetZ2->eta;
            z2_phi=jetZ2->phi;
            z2_mass=jetZ2->mass;
        }

        if(zToLep==1){
            const TGenParticle *z1, *z2;

            z1 = (TGenParticle*) ((*part)[iZ1]);
            z2 = (TGenParticle*) ((*part)[iZ2]);

            z1_pt=z1->pt;
            z1_eta=z1->eta;
            z1_phi=z1->phi;
            z1_mass=z1->mass;

            z2_pt=z2->pt;
            z2_eta=z2->eta;
            z2_phi=z2->phi;
            z2_mass=z2->mass;
        }

        if(iZ!=-1){
            const TGenParticle *parZ = (TGenParticle*) ((*part)[iZ]);

            z_pt = parZ->pt;
            z_eta = parZ->eta;
            z_phi = parZ->phi;
            z_mass = parZ->mass;
        }

        recoz_pt = vRecoZ.Pt();
        recoz_eta = vRecoZ.Eta();
        recoz_phi = vRecoZ.Phi();
        recoz_mass = vRecoZ.M();

        if(iT1!=-1 && iT2!=-1){
            const TGenJet *taujet1, *taujet2;

            taujet1 = (TGenJet*) ((*jet)[iT1]);
            taujet2 = (TGenJet*) ((*jet)[iT2]);

            jetTau1_pt=taujet1->pt;
            jetTau1_eta=taujet1->eta;
            jetTau1_phi=taujet1->phi;
            jetTau1_mass=taujet1->mass;

            jetTau2_pt=taujet2->pt;
            jetTau2_eta=taujet2->eta;
            jetTau2_phi=taujet2->phi;
            jetTau2_mass=taujet2->mass;
        } 

        genTau1_pt=genTau1->pt;
        genTau1_eta=genTau1->eta;
        genTau1_phi=genTau1->phi;
        genTau1_mass=genTau1->mass;

        genTau2_pt=genTau2->pt;
        genTau2_eta=genTau2->eta;
        genTau2_phi=genTau2->phi;
        genTau2_mass=genTau2->mass;

        visTau1_pt=vVisTau1.Pt();
        visTau1_eta=vVisTau1.Eta();
        visTau1_phi=vVisTau1.Phi();
        visTau1_mass=vVisTau1.M();

        visTau2_pt=vVisTau2.Pt();
        visTau2_eta=vVisTau2.Eta();
        visTau2_phi=vVisTau2.Phi();
        visTau2_mass=vVisTau2.M();

        cpions1_pt=cpions1.Pt();
        cpions1_eta=cpions1.Eta();
        cpions1_phi=cpions1.Phi();
        cpions1_mass=cpions1.M();

        cpions2_pt=cpions2.Pt();
        cpions2_eta=cpions2.Eta();
        cpions2_phi=cpions2.Phi();
        cpions2_mass=cpions2.M();

        npions1_pt=npions1.Pt();
        npions1_eta=npions1.Eta();
        npions1_phi=npions1.Phi();
        npions1_mass=npions1.M();

        npions2_pt=npions2.Pt();
        npions2_eta=npions2.Eta();
        npions2_phi=npions2.Phi();
        npions2_mass=npions2.M();

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            if ( genloop->pt<20) continue;
            if ( abs(genloop->pdgId)!=13 && abs(genloop->pdgId)!=11 ) continue;
            nLeptons++;
        }

        outtree->Fill();

    }

    outfile->Write();
    outfile->Save();

}

Double_t deltaR( const Double_t eta1, const Double_t eta2, const Double_t phi1, const Double_t phi2 )
{

    const Double_t pi = 3.14159265358979;

    Double_t etaDiff = (eta1-eta2);
    Double_t phiDiff = fabs(phi1-phi2);
    while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

    Double_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

    return TMath::Sqrt(deltaRSquared);

}