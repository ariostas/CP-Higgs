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

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 );

using namespace baconhep;

void select(const TString sample="ee", const TString tempinput="/afs/cern.ch/work/a/arapyan/public/forMarkus/higgs_cp_samples/ZH_tautau_rhorho_CP0_mad.root", const Float_t xsec = 1,
    const Int_t eosflag = 1)
{
    TString input;
    if(eosflag==1) input = "root://eoscms.cern.ch//store/user/arapyan/higgs_cp_samples/" + sample + "/Bacon/" + tempinput;
    else input = "/afs/cern.ch/work/a/arapyan/public/forMarkus/" + sample + "/" + tempinput;

    TString output = "/afs/cern.ch/work/a/ariostas/private/CP-Higgs_temp/" + sample + "/" + tempinput;

    TChain chain("Events");
    chain.Add(input);

    const Float_t DR=0.4;

    // Data structures to store info from TTrees
    TGenEventInfo *info = new TGenEventInfo();
    TClonesArray *jet             = new TClonesArray("baconhep::TGenJet");
    TClonesArray *part            = new TClonesArray("baconhep::TGenParticle");

    chain.SetBranchAddress("GenEvtInfo",  &info);
    TBranch *infoBr     = chain.GetBranch("GenEvtInfo");
    chain.SetBranchAddress("GenJet" ,     &jet );
    TBranch *jetBr      = chain.GetBranch("GenJet");
    chain.SetBranchAddress("GenParticle", &part);
    TBranch *partBr     = chain.GetBranch("GenParticle");

    // Output file and trees
    TFile *outfile = new TFile(output, "RECREATE");

    TTree *infoTree = new TTree("Count", "Count");
    Long64_t n = chain.GetEntries();
    infoTree->Branch("n", &n, "n/i");
    infoTree->Fill();

    TTree *outtree = new TTree("Events", "Events");

    Float_t eventWeight;

    UInt_t nProngTau1=0, nProngTau2=0, nLeptons=0;

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

    outtree->Branch("eventWeight",      &eventWeight,       "eventWeight/f");   // event weight from cross-section and Event->Weight

    outtree->Branch("nLeptons",         &nLeptons,          "nLeptons/i");      // number of leptons
    outtree->Branch("nProngTau1",       &nProngTau1,        "nProngTau1/i");
    outtree->Branch("nProngTau2",       &nProngTau2,        "nProngTau2/i");

    outtree->Branch("jetTau1_pt",       &jetTau1_pt,        "jetTau1_pt/f");    // pt(Tau1)
    outtree->Branch("jetTau1_eta",      &jetTau1_eta,       "jetTau1_eta/f");   // eta(Tau1)
    outtree->Branch("jetTau1_phi",      &jetTau1_phi,       "jetTau1_phi/f");   // phi(Tau1)
    outtree->Branch("jetTau1_mass",     &jetTau1_mass,      "jetTau1_mass/f");  // m(Tau1)

    outtree->Branch("jetTau2_pt",       &jetTau2_pt,        "jetTau2_pt/f");    // pt(Tau2)
    outtree->Branch("jetTau2_eta",      &jetTau2_eta,       "jetTau2_eta/f");   // eta(Tau2)
    outtree->Branch("jetTau2_phi",      &jetTau2_phi,       "jetTau2_phi/f");   // phi(Tau2)
    outtree->Branch("jetTau2_mass",     &jetTau2_mass,      "jetTau2_mass/f");  // m(Tau2)

    outtree->Branch("genTau1_pt",       &genTau1_pt,        "genTau1_pt/f");    // pt(Tau1)
    outtree->Branch("genTau1_eta",      &genTau1_eta,       "genTau1_eta/f");   // eta(Tau1)
    outtree->Branch("genTau1_phi",      &genTau1_phi,       "genTau1_phi/f");   // phi(Tau1)
    outtree->Branch("genTau1_mass",     &genTau1_mass,      "genTau1_mass/f");  // m(Tau1)

    outtree->Branch("genTau2_pt",       &genTau2_pt,        "genTau2_pt/f");    // pt(Tau2)
    outtree->Branch("genTau2_eta",      &genTau2_eta,       "genTau2_eta/f");   // eta(Tau2)
    outtree->Branch("genTau2_phi",      &genTau2_phi,       "genTau2_phi/f");   // phi(Tau2)
    outtree->Branch("genTau2_mass",     &genTau2_mass,      "genTau2_mass/f");  // m(Tau2)

    outtree->Branch("visTau1_pt",       &visTau1_pt,        "visTau1_pt/f");    // pt(Tau1)
    outtree->Branch("visTau1_eta",      &visTau1_eta,       "visTau1_eta/f");   // eta(Tau1)
    outtree->Branch("visTau1_phi",      &visTau1_phi,       "visTau1_phi/f");   // phi(Tau1)
    outtree->Branch("visTau1_mass",     &visTau1_mass,      "visTau1_mass/f");  // m(Tau1)

    outtree->Branch("visTau2_pt",       &visTau2_pt,        "visTau2_pt/f");    // pt(Tau2)
    outtree->Branch("visTau2_eta",      &visTau2_eta,       "visTau2_eta/f");   // eta(Tau2)
    outtree->Branch("visTau2_phi",      &visTau2_phi,       "visTau2_phi/f");   // phi(Tau2)
    outtree->Branch("visTau2_mass",     &visTau2_mass,      "visTau2_mass/f");  // m(Tau2)

    outtree->Branch("pions1_pt",        &pions1_pt,          "pions1_pt/f");
    outtree->Branch("pions1_eta",       &pions1_eta,         "pions1_eta/f");
    outtree->Branch("pions1_phi",       &pions1_phi,         "pions1_phi/f");
    outtree->Branch("pions1_mass",      &pions1_mass,        "pions1_mass/f");

    outtree->Branch("pions2_pt",        &pions2_pt,          "pions2_pt/f");
    outtree->Branch("pions2_eta",       &pions2_eta,         "pions2_eta/f");
    outtree->Branch("pions2_phi",       &pions2_phi,         "pions2_phi/f");
    outtree->Branch("pions2_mass",      &pions2_mass,        "pions2_mass/f");

    outtree->Branch("neutpions1_pt",    &neutpions1_pt,      "neutpions1_pt/f");
    outtree->Branch("neutpions1_eta",   &neutpions1_eta,     "neutpions1_eta/f");
    outtree->Branch("neutpions1_phi",   &neutpions1_phi,     "neutpions1_phi/f");
    outtree->Branch("neutpions1_mass",  &neutpions1_mass,    "neutpions1_mass/f");

    outtree->Branch("neutpions2_pt",    &neutpions2_pt,      "neutpions2_pt/f");
    outtree->Branch("neutpions2_eta",   &neutpions2_eta,     "neutpions2_eta/f");
    outtree->Branch("neutpions2_phi",   &neutpions2_phi,     "neutpions2_phi/f");
    outtree->Branch("neutpions2_mass",  &neutpions2_mass,    "neutpions2_mass/f");

    for (Int_t i=0; i<chain.GetEntries(); i++)
    {
        if(i%1000==0) cout << "Processing event " << i << endl;

        infoBr->GetEntry(i);

        part->Clear();
        partBr->GetEntry(i);
        jet->Clear();
        jetBr->GetEntry(i);

        if (part->GetEntries()==0) continue;

        Int_t iTau1=-1, iTau2=-1;
        Int_t iT1=-1, iT2=-1;

        // Reset variables
        nLeptons=nProngTau1=nProngTau2=0;

        jetTau1_pt=jetTau1_eta=jetTau1_phi=jetTau1_mass=0;
        jetTau2_pt=jetTau2_eta=jetTau2_phi=jetTau2_mass=0;

        genTau1_pt=genTau1_eta=genTau1_phi=genTau1_mass=0;
        genTau2_pt=genTau2_eta=genTau2_phi=genTau2_mass=0;

        visTau1_pt=visTau1_eta=visTau1_phi=visTau1_mass=0;
        visTau2_pt=visTau2_eta=visTau2_phi=visTau2_mass=0;

        pions1_pt=pions1_eta=pions1_phi=pions1_mass=0;
        pions2_pt=pions2_eta=pions2_phi=pions2_mass=0;

        neutpions1_pt=neutpions1_eta=neutpions1_phi=neutpions1_mass=0;
        neutpions2_pt=neutpions2_eta=neutpions2_phi=neutpions2_mass=0;

        eventWeight=xsec*3000000;

        Int_t isLep=0;
        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            Int_t parentPdg=dynamic_cast<TGenParticle *>(part->At(genloop->parent>-1 ? genloop->parent : 0))->pdgId;

            if (!(fabs(genloop->pdgId)==15)&&!(fabs(parentPdg)==15)) continue;
            if ( (fabs(genloop->pdgId)==13||fabs(genloop->pdgId)==11||fabs(genloop->pdgId)==24) && (fabs(parentPdg)==15) )
            {
                isLep=1;
                continue;
            }

            if (fabs(genloop->pdgId)==15)
            {
                if (iTau1==-1)
                {
                    iTau1=j;
                }
                else if (genloop->parent==iTau1)
                {
                    iTau1=j;
                }
                else if (iTau2==-1)
                {
                    iTau2=j;
                }
                else if (genloop->parent==iTau2)
                {
                    iTau2=j;
                }
            }
        }

        if (iTau2==-1) continue;
        if (isLep!=0) continue;

        if(dynamic_cast<TGenParticle *>(part->At(iTau1))->pdgId==15 &&  dynamic_cast<TGenParticle *>(part->At(iTau2))->pdgId==-15)
        {
            Int_t temp=iTau1;
            iTau1=iTau2;
            iTau2=temp;
        }
        else if(dynamic_cast<TGenParticle *>(part->At(iTau1))->pdgId==-15 &&  dynamic_cast<TGenParticle *>(part->At(iTau2))->pdgId==15)
        {
            Int_t temp=iTau1;
        }
        else
        {
            cout << "Ignoring event." << endl;
            continue;
        }

        vector<Int_t> tau1par, tau2par;
        tau1par.clear(); tau2par.clear();

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            TGenParticle* genloop = (TGenParticle*) ((*part)[j]);

            Int_t pdg=fabs(genloop->pdgId);
            Int_t parent=genloop->parent;

            if((pdg==213 || pdg==20213) && parent==iTau1)
            {
                tau1par.push_back(j);
            }
            else if((pdg==213 || pdg==20213) && parent==iTau2)
            {
                tau2par.push_back(j);
            }

        }

        TLorentzVector pions1, pions2, neutpions1, neutpions2;
        pions1.SetPtEtaPhiE(0,0,0,0); pions2.SetPtEtaPhiE(0,0,0,0);
        neutpions1.SetPtEtaPhiE(0,0,0,0); neutpions2.SetPtEtaPhiE(0,0,0,0);

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            TGenParticle* genloop = (TGenParticle*) ((*part)[j]);

            Int_t pdg=abs(genloop->pdgId);
            Int_t parent=genloop->parent;

            if(pdg!=211 && pdg!=111) continue;

            Int_t parentTau=0;

            if(parent==iTau1) parentTau=1;
            else if(parent==iTau2) parentTau=2;
            else
            {
                for(UInt_t x=0; x<tau1par.size(); x++)
                {
                    if(parent==tau1par.at(x)) parentTau=1;
                }
                for(UInt_t x=0; x<tau2par.size(); x++)
                {
                    if(parent==tau2par.at(x)) parentTau=2;
                }
            }

            if(parentTau==0) continue;

            TLorentzVector vTemp;
            vTemp.SetPtEtaPhiM(genloop->pt, genloop->eta, genloop->phi, genloop->mass);

            if(parentTau==1)
            {
                if(pdg==211){
                    pions1+=vTemp;
                    nProngTau1++;
                }
                else{
                    neutpions1+=vTemp;
                }
            }
            else
            {
                if(pdg==211){
                    pions2+=vTemp;
                    nProngTau2++;
                }
                else{
                    neutpions2+=vTemp;
                }
            }

        }

        const TGenParticle* genTau1 = (TGenParticle*) ((*part)[iTau1]);
        TLorentzVector vGenTau1;
        vGenTau1.SetPtEtaPhiM(genTau1->pt,genTau1->eta,genTau1->phi,genTau1->mass);

        const TGenParticle* genTau2 = (TGenParticle*) ((*part)[iTau2]);
        TLorentzVector vGenTau2;
        vGenTau2.SetPtEtaPhiM(genTau2->pt,genTau2->eta,genTau2->phi,genTau2->mass);

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

        /*for (Int_t j=0; j<jet->GetEntries(); j++)
        {
            const TGenJet* loop = (TGenJet*) ((*jet)[j]);

            if (iT1==-1)
            {
                if (deltaR(loop->eta,vGenTau1.Eta(),loop->phi,vGenTau1.Phi())<DR) iT1=j;
            }
            if (iT2==-1)
            {
                if (deltaR(loop->eta,vGenTau2.Eta(),loop->phi,vGenTau2.Phi())<DR) iT2=j;
            }
        }

        if (iT1==-1|| iT2==-1) continue;

        const TGenJet *taujet1 = (TGenJet*) ((*jet)[iT1]);
        const TGenJet *taujet2 = (TGenJet*) ((*jet)[iT2]);

        jetTau1_pt=taujet1->pt;
        jetTau1_eta=taujet1->eta;
        jetTau1_phi=taujet1->phi;
        jetTau1_mass=taujet1->mass;

        jetTau2_pt=taujet2->pt;
        jetTau2_eta=taujet2->eta;
        jetTau2_phi=taujet2->phi;
        jetTau2_mass=taujet2->mass;*/

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

        pions1_pt=pions1.Pt();
        pions1_eta=pions1.Eta();
        pions1_phi=pions1.Phi();
        pions1_mass=pions1.M();

        pions2_pt=pions2.Pt();
        pions2_eta=pions2.Eta();
        pions2_phi=pions2.Phi();
        pions2_mass=pions2.M();

        neutpions1_pt=neutpions1.Pt();
        neutpions1_eta=neutpions1.Eta();
        neutpions1_phi=neutpions1.Phi();
        neutpions1_mass=neutpions1.M();

        neutpions2_pt=neutpions2.Pt();
        neutpions2_eta=neutpions2.Eta();
        neutpions2_phi=neutpions2.Phi();
        neutpions2_mass=neutpions2.M();

        for (Int_t j=0; j<part->GetEntries(); j++)
        {
            const TGenParticle* genloop = (TGenParticle*) ((*part)[j]);
            if ( genloop->pt<20) continue;
            if ( fabs(genloop->pdgId)!=13&&fabs(genloop->pdgId)!=11 ) continue;
            nLeptons++;
        }

        outtree->Fill();

    }

    outfile->Write();
    outfile->Save();

}

Float_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 )
{

    const Float_t pi = 3.14159265358979;

    Float_t etaDiff = (eta1-eta2);
    Float_t phiDiff = fabs(phi1-phi2);
    while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

    Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

    return TMath::Sqrt(deltaRSquared);

}