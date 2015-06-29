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
#include <TString.h>

#include <TLorentzVector.h>
#include <TVector3.h>

#include <TRandom3.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include "Dataset.h"
#include "CalcT.h"

#include "neutrinos.h"

using namespace std;
using namespace TMath;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(vector<TH1D*>, vector<TString>, TCanvas*, const TString, const TString, const TString);
void histogramS(vector<TH1D*>, vector<TString>, TCanvas*, const TString, const TString, const TString);
void saveResults();
void analyze(TString, Double_t, Int_t);
Double_t deltaR(const Float_t, const Float_t, const Float_t, const Float_t);
Double_t getTheta(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
Double_t getPhi(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
Double_t getDelta(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
TLorentzVector getNeut1(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, Int_t);
TLorentzVector getNeut2(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);

// Initialize histograms
TH1D *Theta_obs0 = new TH1D("Theta_obs0", "Theta_obs0", 20, -3.1416, 3.1416);
TH1D *Theta_obspi4 = new TH1D("Theta_obspi4", "Theta_obspi4", 20, -3.1416, 3.1416);
TH1D *Theta_obspi2 = new TH1D("Theta_obspi2", "Theta_obspi2", 20, -3.1416, 3.1416);
TH1D *Theta_obs3pi4 = new TH1D("Theta_obs3pi4", "Theta_obs3pi4", 20, -3.1416, 3.1416);
TH1D *Theta_sig0 = new TH1D("Theta_sig0", "Theta_sig0", 20, -3.1416, 3.1416);
TH1D *Theta_sigpi2 = new TH1D("Theta_sigpi2", "Theta_sigpi2", 20, -3.1416, 3.1416);
TH1D *Theta_dy = new TH1D("Theta_dy", "Theta_dy", 20, -3.1416, 3.1416);
TH1D *Theta_WW = new TH1D("Theta_WW", "Theta_WW", 20, -3.1416, 3.1416);
TH1D *Theta_ZZ = new TH1D("Theta_ZZ", "Theta_ZZ", 20, -3.1416, 3.1416);
TH1D *Theta_ZZee = new TH1D("Theta_ZZee", "Theta_ZZee", 20, -3.1416, 3.1416);

TH1D *hThetaS0 = new TH1D("hThetaS1", "hThetaS1", 20, -3.1416, 3.1416);
TH1D *hThetaSpi4 = new TH1D("hThetaS2", "hThetaS2", 20, -3.1416, 3.1416);
TH1D *hThetaSpi2 = new TH1D("hThetaS3", "hThetaS3", 20, -3.1416, 3.1416);
TH1D *hThetaS3pi4 = new TH1D("hThetaS4", "hThetaS4", 20, -3.1416, 3.1416);
TH1D *hThetaB = new TH1D("hThetaB", "hThetaB", 20, -3.1416, 3.1416);

TH1D *hPhiS0 = new TH1D("hPhiS0", "hPhiS0", 20, -3.1416, 3.1416);
TH1D *hPhiSpi4 = new TH1D("hPhiSpi4", "hPhiSpi4", 20, -3.1416, 3.1416);
TH1D *hPhiSpi2 = new TH1D("hPhiSpi2", "hPhiSpi2", 20, -3.1416, 3.1416);
TH1D *hPhiS3pi4 = new TH1D("hPhiS3pi4", "hPhiS3pi4", 20, -3.1416, 3.1416);
TH1D *hPhiB = new TH1D("hPhiB", "hPhiB", 20, -3.1416, 3.1416);

TH1D *genhistoS1 = new TH1D("histoS1", "histoS1", 50, 0, 150);
TH1D *genhistoS2 = new TH1D("histoS2", "histoS2", 50, 0, 150);
TH1D *genhistoS3 = new TH1D("histoS3", "histoS3", 50, 0, 150);
TH1D *genhistoS4 = new TH1D("histoS4", "histoS4", 50, 0, 150);
TH1D *genhistoB = new TH1D("histoB", "histoB", 50, 0, 150);

TH1D *genhisto1 = new TH1D("genhisto1", "genhisto1", 50, 50, 110);
TH1D *genhisto2 = new TH1D("genhisto2", "genhisto2", 50, 50, 110);

// Initialize data sets
vector<vector<Dataset> > datasets;
const Int_t nDatasets = 10;

// Initialize storage variables
vector<Double_t> total, selection, kinematicCuts, massCuts;
vector<Double_t> totalError, selectionError, kinematicCutsError, massCutsError;
vector<Int_t> signalFlags;
vector<TString> sampleNames;
vector<TH1D*> hTheta,hThetaTemp, genhistos, hPhi, genhistos2;
vector<TString> histogramNames, histogramNames2;

/*
 * MAIN FUNCTION
 */

 void cphiggs(TString sample = "all", TString inputFile = "xsec.txt"){
    
    cout << "\n\nStarting process...\n\n";

    hThetaTemp.push_back(hThetaS0);
    hThetaTemp.push_back(hThetaSpi4);
    hThetaTemp.push_back(hThetaSpi2);
    hThetaTemp.push_back(hThetaS3pi4);
    hThetaTemp.push_back(Theta_dy);
    hThetaTemp.push_back(Theta_WW);
    hThetaTemp.push_back(Theta_ZZ);
    hThetaTemp.push_back(Theta_ZZee);

    hTheta.push_back(hThetaB);
    hTheta.push_back(hThetaS0);
    hTheta.push_back(hThetaSpi4);
    hTheta.push_back(hThetaSpi2);
    hTheta.push_back(hThetaS3pi4);

    hPhi.push_back(hPhiB);
    hPhi.push_back(hPhiS0);
    hPhi.push_back(hPhiSpi4);
    hPhi.push_back(hPhiSpi2);
    hPhi.push_back(hPhiS3pi4);

    genhistos.push_back(genhistoB);
    genhistos.push_back(genhistoS1);
    genhistos.push_back(genhistoS2);
    genhistos.push_back(genhistoS3);
    genhistos.push_back(genhistoS4);

    genhistos2.push_back(genhisto1);
    genhistos2.push_back(genhisto2);

    histogramNames2.push_back("Real");
    histogramNames2.push_back("From jets/leptons");

    histogramNames.push_back("Backgrounds");
    histogramNames.push_back("Signal (#Delta=0)");
    histogramNames.push_back("Signal (#Delta=#pi/4)");
    histogramNames.push_back("Signal (#Delta=#pi/2)");
    histogramNames.push_back("Signal (#Delta=3#pi/4)");

    vector<Dataset> tempVectorDataset;
    for(Int_t i = 0; i < nDatasets; i++){
        Dataset tempDataset(1,10000);
        tempVectorDataset.push_back(tempDataset);
    }

    for(Int_t i = 0; i < 4; i++){
        datasets.push_back(tempVectorDataset);
    }
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}
    
    Int_t nSamples = 0;
    TString sampleName, sampleBin, crossSection, sampleNumber;
    
    while(ifs >> sampleName >> sampleBin >> crossSection >> sampleNumber){
        
        if(sampleName == "#") continue;
        
        if(sample != sampleName && sample != "all") continue;
        
        Int_t matchedSample = sampleNames.size();
        for(UInt_t x=0;x<sampleNames.size();x++){
            
            if(sampleName == sampleNames.at(x)) matchedSample = x;
            
        }
        
        nSamples = matchedSample;
        
        if(nSamples == sampleNames.size()){
            
            sampleNames.push_back(sampleName);
            signalFlags.push_back(atof(string(sampleNumber).c_str()));
            total.push_back(0);         totalError.push_back(0);
            selection.push_back(0);     selectionError.push_back(0);
            kinematicCuts.push_back(0); kinematicCutsError.push_back(0);
            massCuts.push_back(0);      massCutsError.push_back(0);
            
        }
        
        analyze(sampleBin, atof(string(crossSection).c_str()), nSamples);
        
    }
    
    // Save results
    saveResults();
    
}

void analyze(TString inputfile, Double_t crossSection, Int_t samp)
{
    const TString inputFileTemp = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs/" + inputfile + ".root";

    inputfile = "Reading " + inputfile + " events... ";

    inputfile.Resize(60);

    cout << inputfile << endl;

    vector<vector<Double_t> > finalEvents;

    TRandom3 *random = new TRandom3();
    random->SetSeed(0);

    // Set up storage variables
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

    TLorentzVector vTau1Sol1, vTau1Sol2, vTau2Sol1, vTau2Sol2, vcpion1, vcpion2, vnpion1, vnpion2, vHiggs, vrho1, vrho2;

    Double_t thetaSol1, thetaSol2;
    Double_t phiSol1, phiSol2;

    TFile* infile = new TFile(inputFileTemp);
    assert(infile);
    TTree* intree = (TTree*) infile->Get("Events");
    assert(intree);

    intree->SetBranchAddress("nLeptons",        &nLeptons);
    intree->SetBranchAddress("eventWeight",     &eventWeight);
    intree->SetBranchAddress("ncpions1",        &ncpions1);
    intree->SetBranchAddress("ncpions2",        &ncpions2);
    intree->SetBranchAddress("nnpions1",        &nnpions1);
    intree->SetBranchAddress("nnpions2",        &nnpions2);
    intree->SetBranchAddress("zToLep",          &zToLep);
    intree->SetBranchAddress("genTau1_pt",      &genTau1_pt);
    intree->SetBranchAddress("genTau1_eta",     &genTau1_eta);
    intree->SetBranchAddress("genTau1_phi",     &genTau1_phi);
    intree->SetBranchAddress("genTau1_mass",    &genTau1_mass);
    intree->SetBranchAddress("genTau2_pt",      &genTau2_pt);
    intree->SetBranchAddress("genTau2_eta",     &genTau2_eta);
    intree->SetBranchAddress("genTau2_phi",     &genTau2_phi);
    intree->SetBranchAddress("genTau2_mass",    &genTau2_mass);
    intree->SetBranchAddress("cpions1_pt",      &cpions1_pt);
    intree->SetBranchAddress("cpions1_eta",     &cpions1_eta);
    intree->SetBranchAddress("cpions1_phi",     &cpions1_phi);
    intree->SetBranchAddress("cpions1_mass",    &cpions1_mass);
    intree->SetBranchAddress("cpions2_pt",      &cpions2_pt);
    intree->SetBranchAddress("cpions2_eta",     &cpions2_eta);
    intree->SetBranchAddress("cpions2_phi",     &cpions2_phi);
    intree->SetBranchAddress("cpions2_mass",    &cpions2_mass);
    intree->SetBranchAddress("npions1_pt",      &npions1_pt);
    intree->SetBranchAddress("npions1_eta",     &npions1_eta);
    intree->SetBranchAddress("npions1_phi",     &npions1_phi);
    intree->SetBranchAddress("npions1_mass",    &npions1_mass);
    intree->SetBranchAddress("npions2_pt",      &npions2_pt);
    intree->SetBranchAddress("npions2_eta",     &npions2_eta);
    intree->SetBranchAddress("npions2_phi",     &npions2_phi);
    intree->SetBranchAddress("npions2_mass",    &npions2_mass);
    intree->SetBranchAddress("z1_pt",           &z1_pt);
    intree->SetBranchAddress("z1_eta",          &z1_eta);
    intree->SetBranchAddress("z1_phi",          &z1_phi);
    intree->SetBranchAddress("z1_mass",         &z1_mass);
    intree->SetBranchAddress("z2_pt",           &z2_pt);
    intree->SetBranchAddress("z2_eta",          &z2_eta);
    intree->SetBranchAddress("z2_phi",          &z2_phi);
    intree->SetBranchAddress("z2_mass",         &z2_mass);
    intree->SetBranchAddress("z_pt",            &z_pt);
    intree->SetBranchAddress("z_eta",           &z_eta);
    intree->SetBranchAddress("z_phi",           &z_phi);
    intree->SetBranchAddress("z_mass",          &z_mass);

    Double_t tempSelection=0, tempSelectionError=0;

    for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++)   // Event loop
    {
        intree->GetEntry(iEntry);

        if(cpions1_pt==0 || cpions2_pt==0 || npions1_pt==0 || npions2_pt==0 || ncpions1!=1 || ncpions2!=1 || nnpions1!=1 || nnpions2!=1) continue;

        // Charged pions 4-vectors
        vcpion1.SetPtEtaPhiM(cpions1_pt, cpions1_eta, cpions1_phi, cpions1_mass);
        vcpion2.SetPtEtaPhiM(cpions2_pt, cpions2_eta, cpions2_phi, cpions2_mass);

        // Neutral pions 4-vectors
        vnpion1.SetPtEtaPhiM(npions1_pt, npions1_eta, npions1_phi, npions1_mass);
        vnpion2.SetPtEtaPhiM(npions2_pt, npions2_eta, npions2_phi, npions2_mass);

        // Rho 4-vectors
        vrho1=vcpion1+vnpion1;
        vrho2=vcpion2+vnpion2;

        TLorentzVector v1, v2, vZ;

        v1.SetPtEtaPhiM(z1_pt, z1_eta, z1_phi, z1_mass);
        v2.SetPtEtaPhiM(z2_pt, z2_eta, z2_phi, z2_mass);

        TLorentzVector vInit;
        vInit.SetPxPyPzE(0,0,0,240);
        vZ = v1 + v2;

        // If the event has a generated Z use its real 4-vector
        if(z_pt>0.01) vZ.SetPtEtaPhiM(z_pt, z_eta, z_phi, z_mass);

        vHiggs = vInit-vZ;

        // Scale down the backgrounds to simulate a cut on the mass of the Higgs.
        if(signalFlags.at(samp)>0) eventWeight/=2.0;
        eventWeight*=0.66*0.66;

        //if(vHiggs.M() < 120) continue;
        //if(vHiggs.Pt() > 55) continue;

        //if(deltaR(vrho1.Eta(), vrho2.Eta(), vrho1.Phi(), vrho2.Phi()) < 2.25) continue;
        //if(deltaR(v1.Eta(), v2.Eta(), v1.Phi(), v2.Phi()) < 2.0) continue;

        TLorentzVector vNeutrino1Sol1 = getNeut1(vZ, vInit, vrho1, vrho2, 1);
        TLorentzVector vNeutrino2Sol1 = getNeut2(vrho1, vrho2, vNeutrino1Sol1, vZ, vInit);

        TLorentzVector vNeutrino1Sol2 = getNeut1(vZ, vInit, vrho1, vrho2, 2);
        TLorentzVector vNeutrino2Sol2 = getNeut2(vrho1, vrho2, vNeutrino1Sol2, vZ, vInit);

        //if(vNeutrino1Sol1.M() < -20 || vNeutrino2Sol1.M() < -20 || vNeutrino1Sol2.M() < -20 || vNeutrino2Sol2.M() < -20) continue;

        vTau1Sol1 = vrho1 + vNeutrino1Sol1;
        vTau2Sol1 = vrho2 + vNeutrino2Sol1;

        vTau1Sol2 = vrho1 + vNeutrino1Sol2;
        vTau2Sol2 = vrho2 + vNeutrino2Sol2;

        //if(vTau1Sol1.M() < -10 || vTau2Sol1.M() < -10 || vTau1Sol2.M() < -10 || vTau2Sol2.M() < -10) continue;
        if(vTau1Sol1.Pt() < 20 || vTau2Sol1.Pt() < 20 || vTau1Sol2.Pt() < 20 || vTau2Sol2.Pt() < 20) continue;
        if(vTau1Sol1.Pt() > 90 || vTau2Sol1.Pt() > 90 || vTau1Sol2.Pt() > 90 || vTau2Sol2.Pt() > 90) continue;

        // Compute theta variable
        thetaSol1 = getTheta(vcpion1, vnpion1, vcpion2, vnpion2, vTau1Sol1, vTau2Sol1);
        thetaSol2 = getTheta(vcpion1, vnpion1, vcpion2, vnpion2, vTau1Sol2, vTau2Sol2);

        // Compute phi variable
        phiSol1 = getPhi(vcpion1, vnpion1, vcpion2, vnpion2, vTau1Sol1, vTau2Sol1);
        phiSol2 = getPhi(vcpion1, vnpion1, vcpion2, vnpion2, vTau1Sol2, vTau2Sol2);

        // Check whether theta variable is well defined
        if(thetaSol1 != thetaSol1 || thetaSol2 != thetaSol2) continue;

        // Fill histogram for theta variable
        hThetaTemp.at(signalFlags.at(samp)+3)->Fill(thetaSol1, eventWeight/2.0);
        hThetaTemp.at(signalFlags.at(samp)+3)->Fill(thetaSol2, eventWeight/2.0);

        // Check whether phi variable is well defined
        //if(phiSol1 != phiSol1 || phiSol2 != phiSol2) continue;

        // Fill histogram for phi variable
        //hPhi.at(signalFlags.at(samp))->Fill(phiSol1);
        //hPhi.at(signalFlags.at(samp))->Fill(phiSol2);

        //genhistos.at(signalFlags.at(samp))->Fill(vNeutrino1Sol1.Pt());

        vector<double> vars1, vars2;
        vars1.push_back(thetaSol1);
        //vars1.push_back(eventWeight/10.0);
        vars2.push_back(thetaSol2);
        //vars2.push_back(eventWeight/10.0);

        finalEvents.push_back(vars1);
        finalEvents.push_back(vars2);

        tempSelection+=eventWeight;
        tempSelectionError++;

    } // end event loop

    for(Int_t i = 0; i < nDatasets; i++){
        for(Int_t x = 0; x < Nint(tempSelection)/5; x++){
            Int_t randEntry = Nint(random->Uniform(finalEvents.size()/2-1));

            if(signalFlags.at(samp)>0){
                datasets.at(0).at(i).add(finalEvents.at(randEntry));
                datasets.at(0).at(i).add(finalEvents.at(randEntry+1));
                datasets.at(1).at(i).add(finalEvents.at(randEntry));
                datasets.at(1).at(i).add(finalEvents.at(randEntry+1));
                datasets.at(2).at(i).add(finalEvents.at(randEntry));
                datasets.at(2).at(i).add(finalEvents.at(randEntry+1));
                datasets.at(3).at(i).add(finalEvents.at(randEntry));
                datasets.at(3).at(i).add(finalEvents.at(randEntry+1));
            }
            else{
                datasets.at(signalFlags.at(samp)+3).at(i).add(finalEvents.at(randEntry));
                datasets.at(signalFlags.at(samp)+3).at(i).add(finalEvents.at(randEntry+1));
            }
        }
    } 

    TString out = "";
    out += tempSelection;
    out.Resize(9);

    selection.at(samp) += tempSelection;
    if(tempSelectionError > 0) selectionError.at(samp) += sqrtf(tempSelectionError)*tempSelection/tempSelectionError;

    cout << "\e[A";
    cout << inputfile << (signalFlags.at(samp)<=0 ? "\033[1;32m" : (out == "0       " ?  "\033[1;34m": "\033[1;31m")) << out << "\033[0m" << " events passed all cuts" << endl;

    infile->Close();
}

/*
 * FUNCTION FOR PRINTING AND SAVING THE RESULTS
 */

void saveResults()
{
    cout << endl;

    hThetaB->Add(Theta_dy);
    hThetaB->Add(Theta_WW);
    hThetaB->Add(Theta_ZZ);
    hThetaB->Add(Theta_ZZee);
    Theta_sig0->Add(hThetaS0);
    Theta_sigpi2->Add(hThetaSpi2);

    Theta_obs0->Add(Theta_dy);
    Theta_obs0->Add(Theta_WW);
    Theta_obs0->Add(Theta_ZZ);
    Theta_obs0->Add(Theta_ZZee);
    Theta_obs0->Add(hThetaS0);

    Theta_obspi4->Add(Theta_dy);
    Theta_obspi4->Add(Theta_WW);
    Theta_obspi4->Add(Theta_ZZ);
    Theta_obspi4->Add(Theta_ZZee);
    Theta_obspi4->Add(hThetaSpi4);

    Theta_obspi2->Add(Theta_dy);
    Theta_obspi2->Add(Theta_WW);
    Theta_obspi2->Add(Theta_ZZ);
    Theta_obspi2->Add(Theta_ZZee);
    Theta_obspi2->Add(hThetaSpi2);

    Theta_obs3pi4->Add(Theta_dy);
    Theta_obs3pi4->Add(Theta_WW);
    Theta_obs3pi4->Add(Theta_ZZ);
    Theta_obs3pi4->Add(Theta_ZZee);
    Theta_obs3pi4->Add(hThetaS3pi4);

    TH1D *hp1 = new TH1D("hp1", "hp1", 100, 0, 1);
    TH1D *hp2 = new TH1D("hp2", "hp2", 100, 0, 1);
    TH1D *hp3 = new TH1D("hp3", "hp3", 100, 0, 1);
    TH1D *hp4 = new TH1D("hp4", "hp4", 100, 0, 1);

    vector<TH1D*> hp;
    hp.push_back(hp1); hp.push_back(hp2); hp.push_back(hp3); hp.push_back(hp4);
    vector<TString> hpNames;
    hpNames.push_back("Delta=0 and Delta=0"); hpNames.push_back("Delta=0 and Delta=pi/4"); hpNames.push_back("Delta=0 and Delta=pi/2"); hpNames.push_back("Delta=0 and Delta=3pi/4");

    srand(22);
    for(Int_t i=0; i<4; i++){
        Double_t averageT=0, averageP=0, averagePError=0;

        for(Int_t x =0; x < nDatasets; x++){

            Double_t tval = calcT(datasets.at(0).at(x),datasets.at(i).at(x));

            Double_t ptval;
            Int_t pval=0;
            Int_t nperm=100;
            for(Int_t p=0; p<nperm; p++){
                ptval = permCalcT(datasets.at(0).at(x),datasets.at(i).at(x));
                //hp.at(i)->Fill(ptval);
                if(ptval > tval) pval++;
                if(p!=0) cout << "\e[A";
                cout << "[" << string((p+1)/(nperm/50),'-') << string(50-(p+1)/(nperm/50),' ') << "]  " << 100*(p+1)/nperm << "\% completed.   Permutation " << p+1 << " of " << nperm;
                cout << "   Dataset " << x+1 << " of " << nDatasets << endl;
            }
            Double_t p = pval/Double_t(nperm), pError = Sqrt(p*(1-p)/Double_t(nperm));

            hp.at(i)->Fill(p);

            averageT+=tval;
            averageP+=p;
            averagePError+=pError;
        }
        averageT/=Double_t(nDatasets);
        averageP/=Double_t(nDatasets);
        averagePError/=Double_t(nDatasets);

        cout << "Average T of " << hpNames.at(i) << ": " << averageT << endl;
        cout << "Average p of " << hpNames.at(i) << ": " << averageP  << " +- " << averagePError << endl << endl;

    }

    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);

    gStyle->SetOptStat(kFALSE);

    histogram(hTheta, histogramNames, c1, "#Theta variable", "Fraction", "Theta");
    histogram(hPhi, histogramNames, c1, "#varphi^{*} variable", "Fraction", "Phi");
    histogram(genhistos, histogramNames, c1, "H mass", "Fraction", "histo");
    histogram(genhistos2, histogramNames2, c1, "Z mass", "Fraction", "h3");

    for(Int_t i=0; i<4; i++){
        TString filename="p_"; filename+=i; filename+=".jpg";
        histogram(hp.at(i), hpNames.at(i), c1, "p distribution", "Fraction", filename);
    }




    TFile f("Datacards/Histograms.root","new");
    Theta_obs0->Write();
    Theta_obspi4->Write();
    Theta_obspi2->Write();
    Theta_obs3pi4->Write();
    Theta_sig0->Write();
    Theta_sigpi2->Write();
    Theta_dy->Write();
    Theta_WW->Write();
    Theta_ZZ->Write();
    Theta_ZZee->Write();
    f.Close();

    cout << "\n\n\nProcess finished\nPrinting results...\n\n" << "\033[1;34mResults\033[0m\n\n";
    
    for(UInt_t x = 0; x<sampleNames.size();x++){
        
        cout << (signalFlags.at(x)<=0 ? "\033[1;32m" : "\033[1;31m") << sampleNames.at(x) << "\033[0m\n";
        cout << "Events after selection: " << selection.at(x) << " +- " << selectionError.at(x) << endl << endl;
        
    }

    cout << "Done\nExiting...\n\n\n";

}

/*
 * FUNCTION FOR SAVING ONE HISTOGRAM
 */

void histogram(TH1D *histo, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){
    Double_t norm=1;
    norm/=histo->Integral();
    histo->Scale(norm);
    histo->SetLineWidth(3);
    histo->Draw();
    // add axis labels
    histo->GetXaxis()->SetTitle(xTitle);
    histo->GetYaxis()->SetTitle(yTitle);
    histo->SetTitle(histName); // title on top

    can->SaveAs(name);
}

/*
 * FUNCTION FOR SAVING MULTIPLE HISTOGRAMS
 */

void histogram(vector<TH1D*> histos, vector<TString> histNames, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){   
    
    if(histos.size()!=histNames.size()){ cout << "Number of histograms and names don't match." << endl; return;}

    Double_t max=0, min=1;

    for(Int_t i=0; i<histos.size(); i++){
        histos.at(i)->Scale(Double_t(1)/histos.at(i)->Integral());
        if(histos.at(i)->GetMaximum()>max) max=histos.at(i)->GetMaximum();
        if(histos.at(i)->GetMinimum()<min) min=histos.at(i)->GetMinimum();
    }

    max*=1.1;
    min*=0.9;

    vector<Int_t> colors;
    colors.push_back(kRed); colors.push_back(kBlue); colors.push_back(kGreen); colors.push_back(kGreen+3); colors.push_back(kMagenta+2); colors.push_back(kBlack); 

    for(Int_t i=0; i<histos.size(); i++){
        histos.at(i)->SetMaximum(max);
        histos.at(i)->SetMinimum(min);
        histos.at(i)->SetLineWidth(3);
        histos.at(i)->SetLineColor(colors.at(i%6));
        if(i==0){
            histos.at(i)->GetXaxis()->SetTitle(xTitle);
            histos.at(i)->GetYaxis()->SetTitle(yTitle);
            histos.at(i)->SetTitle("");
            histos.at(i)->Draw();
        }
        else histos.at(i)->Draw("same");
        
    }

    TLegend *leg = new TLegend(0.605,0.675,0.885,0.875);
    leg->SetTextFont(72);
    leg->SetTextSize(0.04);
    for(Int_t i=0; i<histos.size(); i++){
        leg->AddEntry(histos.at(i),histNames.at(i),"l");
    }
    leg->Draw();

    TString filename=name; filename+=".jpg";

    can->SaveAs(filename);
}

/*
 * FUNCTION FOR SAVING MULTIPLE HISTOGRAMS SEPARATELY
 */

void histogramS(vector<TH1D*> histos, vector<TString> histNames, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){   
    
    if(histos.size()!=histNames.size()){ cout << "Number of histograms and names don't match." << endl; return;}

    Double_t max0=0, min0=1;

    histos.at(0)->Scale(Double_t(1)/histos.at(0)->Integral());
    max0=histos.at(0)->GetMaximum();
    min0=histos.at(0)->GetMinimum();

    for(Int_t i=1; i<histos.size(); i++){
        Double_t max=max0, min=min0;
        //histos.at(i)->Scale(Double_t(1)/histos.at(i)->Integral());
        if(histos.at(i)->GetMaximum()>max0) max=histos.at(i)->GetMaximum();
        if(histos.at(i)->GetMinimum()<min0) min=histos.at(i)->GetMinimum();

        max*=1.1;
        min*=0.9; 

        histos.at(0)->SetMaximum(max);
        histos.at(0)->SetMinimum(min);
        histos.at(0)->SetLineWidth(3);
        histos.at(0)->SetLineColor(kRed);

        histos.at(i)->SetMaximum(max);
        histos.at(i)->SetMinimum(min);
        histos.at(i)->SetLineWidth(3);
        histos.at(i)->SetLineColor(kBlue);

        histos.at(0)->GetXaxis()->SetTitle(xTitle);
        histos.at(0)->GetYaxis()->SetTitle(yTitle);
        histos.at(0)->SetTitle("");
        histos.at(0)->Draw();

        histos.at(i)->Draw("same");

        TLegend *leg = new TLegend(0.605,0.675,0.885,0.875);
        leg->SetTextFont(72);
        leg->SetTextSize(0.04);
        leg->AddEntry(histos.at(0),histNames.at(0),"l");
        leg->AddEntry(histos.at(i),histNames.at(i),"l");
        leg->Draw();

        TString filename=name; filename+="_"; filename+=i; filename+=".jpg";

        can->SaveAs(filename);

    }
}

/*
 * FUNCTION FOR dR calculation
 */
Double_t deltaR( const Float_t eta1, const Float_t eta2, const Float_t phi1, const Float_t phi2 )
{

    const Float_t pi = 3.14159265358979;

    Float_t etaDiff = (eta1-eta2);
    Float_t phiDiff = fabs(phi1-phi2);
    while ( phiDiff>pi ) phiDiff = fabs(phiDiff-2.0*pi);

    Float_t deltaRSquared = etaDiff*etaDiff + phiDiff*phiDiff;

    return TMath::Sqrt(deltaRSquared);

}

/*
 * FUNCTION FOR COMPUTING THE THETA VARIABLE
 */

Double_t getTheta(TLorentzVector vcpion1, TLorentzVector vnpion1, TLorentzVector vcpion2, TLorentzVector vnpion2, TLorentzVector vTau1, TLorentzVector vTau2){

        Double_t r=0.14;

        TLorentzVector vHiggs = vTau1 + vTau2;

        TLorentzVector vRho1 = vcpion1 + vnpion1;
        TLorentzVector vRho2 = vcpion2 + vnpion2;

        TVector3 v3Higgs;
        v3Higgs.SetXYZ(vHiggs.Px()/vHiggs.E(), vHiggs.Py()/vHiggs.E(), vHiggs.Pz()/vHiggs.E());

        TLorentzVector vq1=vcpion1-vnpion1;
        TLorentzVector vq2=vcpion2-vnpion2;

        Double_t y1=(vq1.Dot(vTau1))/(vRho1.Dot(vTau1));
        Double_t y2=(vq2.Dot(vTau2))/(vRho2.Dot(vTau2));

        // Boost vectors to Higgs frame
        vcpion1.Boost(-v3Higgs);
        vcpion2.Boost(-v3Higgs);
        vnpion1.Boost(-v3Higgs);
        vnpion2.Boost(-v3Higgs);
        vTau1.Boost(-v3Higgs);
        vTau2.Boost(-v3Higgs);

        // Set up 3-vectors in Higgs frame
        TVector3 v3cpion1, v3cpion2, v3npion1, v3npion2;

        TVector3 v3tau1, v3tau2;
        v3tau1.SetXYZ(vTau1.Px()/vTau1.E(), vTau1.Py()/vTau1.E(), vTau1.Pz()/vTau1.E());
        v3tau2.SetXYZ(vTau2.Px()/vTau2.E(), vTau2.Py()/vTau2.E(), vTau2.Pz()/vTau2.E());

        v3cpion1.SetXYZ(vcpion1.Px()/vcpion1.E(), vcpion1.Py()/vcpion1.E(), vcpion1.Pz()/vcpion1.E());
        v3cpion2.SetXYZ(vcpion2.Px()/vcpion2.E(), vcpion2.Py()/vcpion2.E(), vcpion2.Pz()/vcpion2.E());
        v3npion1.SetXYZ(vnpion1.Px()/vnpion1.E(), vnpion1.Py()/vnpion1.E(), vnpion1.Pz()/vnpion1.E());
        v3npion2.SetXYZ(vnpion2.Px()/vnpion2.E(), vnpion2.Py()/vnpion2.E(), vnpion2.Pz()/vnpion2.E());

        TVector3 tempE1, tempE2;

        tempE1=(y1-r)*v3cpion1-(y1+r)*v3npion1;
        tempE2=(y2-r)*v3cpion2-(y2+r)*v3npion2;

        TVector3 vE1=tempE1-(tempE1.Dot(v3tau1))*v3tau1;
        TVector3 vE2=tempE2-(tempE2.Dot(v3tau2))*v3tau2;

        // Compute theta variable
        return TMath::Sign(Double_t(1),v3tau1.Dot(vE2.Cross(vE1)))*TMath::ACos(vE1.Dot(vE2)/(vE1.Mag()*vE2.Mag()));

 }

 /*
 * FUNCTION FOR COMPUTING THE PHI VARIABLE
 */

Double_t getPhi(TLorentzVector vcpion1, TLorentzVector vnpion1, TLorentzVector vcpion2, TLorentzVector vnpion2, TLorentzVector vTau1, TLorentzVector vTau2){

        TVector3 v3Tau1, v3Tau2;
        v3Tau1.SetXYZ(vTau1.Px()/vTau1.E(), vTau1.Py()/vTau1.E(), vTau1.Pz()/vTau1.E());
        v3Tau2.SetXYZ(vTau2.Px()/vTau2.E(), vTau2.Py()/vTau2.E(), vTau2.Pz()/vTau2.E());

        // Boost vectors to respective tau frames
        vcpion1.Boost(-v3Tau1);
        vcpion2.Boost(-v3Tau2);
        vnpion1.Boost(-v3Tau1);
        vnpion2.Boost(-v3Tau2);

        Double_t y1 = (vcpion1.E()-vnpion1.E())/(vcpion1.E()+vnpion1.E());
        Double_t y2 = (vcpion2.E()-vnpion2.E())/(vcpion2.E()+vnpion2.E());

        bool isPositive;
        isPositive = (y1*y2 > 0 ? true : false);

        // Boost vectors back to lab frame
        vcpion1.Boost(v3Tau1);
        vcpion2.Boost(v3Tau2);
        vnpion1.Boost(v3Tau1);
        vnpion2.Boost(v3Tau2);

        TLorentzVector vHiggs = vTau1 + vTau2;  

        TVector3 v3Higgs;
        v3Higgs.SetXYZ(vHiggs.Px()/vHiggs.E(), vHiggs.Py()/vHiggs.E(), vHiggs.Pz()/vHiggs.E());      

        // Boost vectors to Higgs frame
        vcpion1.Boost(-v3Higgs);
        vcpion2.Boost(-v3Higgs);
        vnpion1.Boost(-v3Higgs);
        vnpion2.Boost(-v3Higgs);
        vTau1.Boost(-v3Higgs);
        vTau2.Boost(-v3Higgs);

        // Set up 3-vectors in Higgs frame
        TVector3 v3cpion1, v3cpion2, v3npion1, v3npion2;

        v3Tau1.SetXYZ(vTau1.Px(), vTau1.Py(), vTau1.Pz());
        v3Tau2.SetXYZ(vTau2.Px(), vTau2.Py(), vTau2.Pz());

        v3cpion1.SetXYZ(vcpion1.Px(), vcpion1.Py(), vcpion1.Pz());
        v3cpion2.SetXYZ(vcpion2.Px(), vcpion2.Py(), vcpion2.Pz());
        v3npion1.SetXYZ(vnpion1.Px(), vnpion1.Py(), vnpion1.Pz());
        v3npion2.SetXYZ(vnpion2.Px(), vnpion2.Py(), vnpion2.Pz());

        TVector3 n1, n2;
        n1 = v3cpion1.Cross(v3Tau1)*(1.0/v3cpion1.Cross(v3Tau1).Mag());
        n2 = v3cpion2.Cross(v3Tau2)*(1.0/v3cpion2.Cross(v3Tau2).Mag());

        Double_t phi = ACos(n1.Dot(n2));

        // Compute phi variable
        return (isPositive ? phi : phi-3.14);

 }

 /*
 * FUNCTION FOR COMPUTING THE DELTA VARIABLE
 */

Double_t getDelta(TLorentzVector vcpion1, TLorentzVector vcpion2, TLorentzVector vTau1, TLorentzVector vTau2){

        TLorentzVector vHiggs = vTau1 + vTau2;

        TVector3 v3Higgs;
        v3Higgs.SetXYZ(vHiggs.Px()/vHiggs.E(), vHiggs.Py()/vHiggs.E(), vHiggs.Pz()/vHiggs.E());

        // Boost vectors to Higgs frame
        vcpion1.Boost(-v3Higgs);
        vcpion2.Boost(-v3Higgs);

        // Set up 3-vectors in Higgs frame
        TVector3 v3cpion1, v3cpion2;

        v3cpion1.SetXYZ(vcpion1.Px()/vcpion1.E(), vcpion1.Py()/vcpion1.E(), vcpion1.Pz()/vcpion1.E());
        v3cpion2.SetXYZ(vcpion2.Px()/vcpion2.E(), vcpion2.Py()/vcpion2.E(), vcpion2.Pz()/vcpion2.E());

        v3cpion1 *= 1.0/v3cpion1.Mag();
        v3cpion2 *= 1.0/v3cpion2.Mag();

        // Compute theta variable
        return ACos(v3cpion1.Dot(v3cpion2));

 }