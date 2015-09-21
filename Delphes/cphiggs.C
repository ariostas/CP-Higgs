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
#include <TH2.h>
#include <TF1.h>
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

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(TH2D*, const TString, TCanvas*, const TString);
void histogram(vector<TH2D*>, const TString, TCanvas*, const TString);
void histogram(vector<TH1D*>, vector<TString>, TCanvas*, const TString, const TString, const TString);
void histogramS(vector<TH1D*>, vector<TString>, TCanvas*, const TString, const TString, const TString);
void histogram(TH1D*, TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void saveResults(TString, TString, TString);
void analyze(TString, Double_t, Int_t);
Double_t deltaR(TLorentzVector, TLorentzVector);
Double_t getTheta(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
TLorentzVector getNeut1(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, Int_t);
TLorentzVector getNeut2(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
Double_t getChiSquared(TH1D*, TH1D*);
Double_t getLikelihood(TH1D*, TH1D*);

// Initialize histograms
TH1D *Theta_obs0 = new TH1D("Theta_obs0", "Theta_obs0", 20, -3.1416, 3.1416);
TH1D *Theta_obspi4 = new TH1D("Theta_obspi4", "Theta_obspi4", 20, -3.1416, 3.1416);
TH1D *Theta_obspi2 = new TH1D("Theta_obspi2", "Theta_obspi2", 20, -3.1416, 3.1416);
TH1D *Theta_obs3pi4 = new TH1D("Theta_obs3pi4", "Theta_obs3pi4", 20, -3.1416, 3.1416);
TH1D *Theta_sig0 = new TH1D("Theta_sig0", "Theta_sig0", 20, -3.1416, 3.1416);
TH1D *Theta_sigpi4 = new TH1D("Theta_sigpi4", "Theta_sigpi4", 20, -3.1416, 3.1416);
TH1D *Theta_sigpi2 = new TH1D("Theta_sigpi2", "Theta_sigpi2", 20, -3.1416, 3.1416);
TH1D *Theta_sig3pi4 = new TH1D("Theta_sig3pi4", "Theta_sig3pi4", 20, -3.1416, 3.1416);
TH1D *Theta_jj = new TH1D("Theta_jj", "Theta_jj", 20, -3.1416, 3.1416);
TH1D *Theta_dy = new TH1D("Theta_dy", "Theta_dy", 20, -3.1416, 3.1416);
TH1D *Theta_ZZ = new TH1D("Theta_ZZ", "Theta_ZZ", 20, -3.1416, 3.1416);
TH1D *Theta_Zll = new TH1D("Theta_Zll", "Theta_Zll", 20, -3.1416, 3.1416);

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

TH1D *genhistoS1 = new TH1D("histoS1", "histoS1", 100, -.5, 9.5);
TH1D *genhistoS2 = new TH1D("histoS2", "histoS2", 100, -.5, 9.5);
TH1D *genhistoS3 = new TH1D("histoS3", "histoS3", 100, -.5, 9.5);
TH1D *genhistoS4 = new TH1D("histoS4", "histoS4", 100, -.5, 9.5);
TH1D *genhistoB = new TH1D("histoB", "histoB", 100, -.5, 9.5);

TH1D *genhisto1 = new TH1D("genhisto1", "genhisto1", 50, -1, 150);
TH1D *genhisto2 = new TH1D("genhisto2", "genhisto2", 50, -1, 150);

TH2D *h2Theta0 = new TH2D("h2Theta0", "h2Theta0", 20, -3.1416, 3.1416, 20, -3.1416, 3.1416);
TH2D *h2Thetapi4 = new TH2D("h2Thetapi4", "h2Thetapi4", 20, -3.1416, 3.1416, 20, -3.1416, 3.1416);
TH2D *h2Thetapi2 = new TH2D("h2Thetapi2", "h2Thetapi2", 20, -3.1416, 3.1416, 20, -3.1416, 3.1416);
TH2D *h2Theta3pi4 = new TH2D("h2Theta3pi4", "h2Theta3pi4", 20, -3.1416, 3.1416, 20, -3.1416, 3.1416);

// Initialize data sets
vector<vector<Dataset> > datasets;
const Int_t nDatasets = 1000;

// Initialize storage variables
vector<Double_t> total, selection, kinematicCuts, massCuts;
vector<Double_t> totalError, selectionError, kinematicCutsError, massCutsError;
vector<Int_t> signalFlags;
vector<TString> sampleNames;
vector<TH1D*> hTheta, hThetaTemp, genhistos, hPhi, genhistos2, hThetaObs;
vector<TString> histogramNames, histogramNames2;
vector<TH2D*> h2Theta;

/*
 * MAIN FUNCTION
 */

void cphiggs(TString calcP = "false", TString calcChi = "true", TString calcL = "true", TString sample = "all", TString inputFile = "xsec.txt"){
    
    cout << "\n\nStarting process...\n\n";

    hThetaObs.push_back(Theta_obs0);
    hThetaObs.push_back(Theta_obspi4);
    hThetaObs.push_back(Theta_obspi2);
    hThetaObs.push_back(Theta_obs3pi4);

    h2Theta.push_back(h2Theta0);
    h2Theta.push_back(h2Thetapi4);
    h2Theta.push_back(h2Thetapi2);
    h2Theta.push_back(h2Theta3pi4);

    hThetaTemp.push_back(hThetaS0);
    hThetaTemp.push_back(hThetaSpi4);
    hThetaTemp.push_back(hThetaSpi2);
    hThetaTemp.push_back(hThetaS3pi4);
    hThetaTemp.push_back(Theta_jj);
    hThetaTemp.push_back(Theta_dy);
    hThetaTemp.push_back(Theta_ZZ);
    hThetaTemp.push_back(Theta_Zll);

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
        Dataset tempDataset(2,10000);
        tempVectorDataset.push_back(tempDataset);
    }

    for(Int_t i = 0; i < 4; i++){
        datasets.push_back(tempVectorDataset);
    }
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}
    
    Int_t nSamples = 0;
    TString sampleName, sampleBin, crossSection, sampleNumber;
    
    while(ifs >> sampleName >> sampleBin >> crossSection >> sampleNumber){
        
        if(sampleName.Contains("#")) continue;
        
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
    saveResults(calcP,calcChi, calcL);
    
}

void analyze(TString inputfile, Double_t xsec, Int_t samp){
    
    TString inputFile = inputfile;

    vector<vector<Double_t> > finalEvents;

    TRandom3 *random = new TRandom3();
    random->SetSeed(0);

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

    const TString inputFileTemp = "/afs/cern.ch/work/a/ariostas/public/CP-Higgs_Samples_small/" + inputfile + ".root";
    inputfile = "Reading " + inputfile + " events... ";
    inputfile.Resize(60);
    cout << inputfile << endl;

    TFile* infile = new TFile(inputFileTemp);
    assert(infile);
    TTree* intree = (TTree*) infile->Get("Events");
    assert(intree);
    Long64_t numberOfEntries = intree->GetEntries();

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

    Double_t tempSelection=0, tempSelectionError=0, tempMassCuts=0, tempMassCutsError=0;

    for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) { // entry loop
        intree->GetEntry(iEntry);

        if(hasH == 1 && signalFlags.at(samp) > 0) continue;

        if(NCPions1 != 1 || NCPions2 != 1) continue;
        if(NNPions1 != 1 || NNPions2 != 1) continue;
        //if(NNPions1 < 1 || NNPions2 < 1) continue;

        TLorentzVector V4Z, V4H, V4CPion1, V4CPion2, V4NPion11, V4NPion12, V4NPion21, V4NPion22, V4NPion1, V4NPion2, V4Init;

        V4Init.SetPxPyPzE(0, 0, 0, 240);

        V4CPion1.SetPtEtaPhiM(CPion1_Pt, CPion1_Eta, CPion1_Phi, CPion1_Mass);
        V4CPion2.SetPtEtaPhiM(CPion2_Pt, CPion2_Eta, CPion2_Phi, CPion2_Mass);

        V4NPion11.SetPtEtaPhiM(NPion11_Pt, NPion11_Eta, NPion11_Phi, NPion11_Mass);
        V4NPion12.SetPtEtaPhiM(NPion12_Pt, NPion12_Eta, NPion12_Phi, NPion12_Mass);
        V4NPion21.SetPtEtaPhiM(NPion21_Pt, NPion21_Eta, NPion21_Phi, NPion21_Mass);
        V4NPion22.SetPtEtaPhiM(NPion22_Pt, NPion22_Eta, NPion22_Phi, NPion22_Mass);

        V4NPion1 = V4NPion11 + V4NPion12;
        V4NPion2 = V4NPion21 + V4NPion22;

        //if(NNPions1 != 1) V4NPion1.SetPtEtaPhiM(NPions1_Pt, NPions1_Eta, NPions1_Phi, NPions1_Mass);
        //if(NNPions2 != 1) V4NPion2.SetPtEtaPhiM(NPions2_Pt, NPions2_Eta, NPions2_Phi, NPions2_Mass);

        Int_t ZFromLep;
        if(ZLepton1_Pt != 0. && ZLepton2_Pt != 0.){
            TLorentzVector V4Temp1, V4Temp2;
            V4Temp1.SetPtEtaPhiM(ZLepton1_Pt, ZLepton1_Eta, ZLepton1_Phi, ZLepton1_Mass);
            V4Temp2.SetPtEtaPhiM(ZLepton2_Pt, ZLepton2_Eta, ZLepton2_Phi, ZLepton2_Mass);
            V4Z = V4Temp1 + V4Temp2;
            ZFromLep = 1;
        }
        else if(ZReco_Pt > .1 && ZReco_Mass > 0.1){
            V4Z.SetPtEtaPhiM(ZReco_Pt, ZReco_Eta, ZReco_Phi, ZReco_Mass);
            ZFromLep = -1;
        }
        // else if(ZJet1_Pt != 0. && ZJet2_Pt != 0.){
        //    TLorentzVector V4Temp1, V4Temp2;
        //     V4Temp1.SetPtEtaPhiM(ZJet1_Pt, ZJet1_Eta, ZJet1_Phi, ZJet1_Mass);
        //     V4Temp2.SetPtEtaPhiM(ZJet2_Pt, ZJet2_Eta, ZJet2_Phi, ZJet2_Mass);
        //     V4Z = V4Temp1 + V4Temp2;
        //     ZFromLep = 0;
        // }
        else{
            continue;
        }

        //if(ZFromLep != -1) continue;

        //tempSelection += eventWeight;
        tempSelectionError++;

        V4H = V4Init - V4Z;

        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(V4NPion1.M(), eventWeight);
        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(V4NPion2.M(), eventWeight);
        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(V4H.M(), eventWeight);

        if(fabs(V4Z.M() - 91.2) > 5.) continue;

        if(fabs(V4H.M() - 125) > 5.) continue;

        // Rho 4-vectors
        TLorentzVector V4Rho1, V4Rho2, V4Neutrino1Sol1, V4Neutrino1Sol2, V4Neutrino2Sol1, V4Neutrino2Sol2, V4Tau1Sol1, V4Tau1Sol2, V4Tau2Sol1, V4Tau2Sol2;

        V4Rho1=V4CPion1+V4NPion1;
        V4Rho2=V4CPion2+V4NPion2;
        V4Neutrino1Sol1 = getNeut1(V4Z, V4Init, V4Rho1, V4Rho2, 1);
        V4Neutrino2Sol1 = getNeut2(V4Rho1, V4Rho2, V4Neutrino1Sol1, V4Z, V4Init);
        V4Neutrino1Sol2 = getNeut1(V4Z, V4Init, V4Rho1, V4Rho2, 2);
        V4Neutrino2Sol2 = getNeut2(V4Rho1, V4Rho2, V4Neutrino1Sol2, V4Z, V4Init);
        V4Tau1Sol1 = V4Rho1 + V4Neutrino1Sol1;
        V4Tau1Sol2 = V4Rho1 + V4Neutrino1Sol2;
        V4Tau2Sol1 = V4Rho2 + V4Neutrino2Sol1;
        V4Tau2Sol2 = V4Rho2 + V4Neutrino2Sol2;

        if(V4Neutrino1Sol1.M() < -5 || V4Neutrino2Sol1.M() < -5 || V4Neutrino1Sol2.M() < -5 || V4Neutrino2Sol2.M() < -5) continue;

        tempSelection += eventWeight;

        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill((V4CPion1 + V4NPion1).M(), eventWeight);
        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill((V4CPion2 + V4NPion2).M(), eventWeight);

        // genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion1, V4NPion11), eventWeight);
        // genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion1, V4NPion12), eventWeight);
        // genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion2, V4NPion21), eventWeight);
        // genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion2, V4NPion22), eventWeight);

        TLorentzVector V4Temp1, V4Temp2;
        V4Temp1.SetPtEtaPhiM(ZLepton1_Pt, ZLepton1_Eta, ZLepton1_Phi, ZLepton1_Mass);
        V4Temp2.SetPtEtaPhiM(ZLepton2_Pt, ZLepton2_Eta, ZLepton2_Phi, ZLepton2_Mass);

        if(ZFromLep == 1){
            if(sameCharge != 0 || NLeptons != 2) continue;
            if(deltaR(V4CPion1, V4Temp1) < 0.4 || deltaR(V4CPion1, V4Temp2) < 0.4 || deltaR(V4CPion2, V4Temp1) < 0.4 || deltaR(V4CPion2, V4Temp2) < 0.4) continue;
            if(fabs(V4Z.P()-51.6) > 2.) continue;
        }
        else if(ZFromLep == -1){
            //if(fabs(V4Z.P()-51.6) > 3.) continue;
            if(NLeptons != 0) continue;
            //if(NJets < 2) continue;
            if(fabs(V4Z.P()-51.6) > 3. || V4Z.Pt() < 10) continue;
            if(V4Tau1Sol1.M() < 0 || V4Tau2Sol1.M() < 0 || V4Tau1Sol2.M() < 0 || V4Tau2Sol2.M() < 0) continue;
        }
        
        genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(sameCharge, eventWeight);

        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion1, V4NPion11), eventWeight);
        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion1, V4NPion12), eventWeight);
        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion2, V4NPion21), eventWeight);
        //genhistos.at(signalFlags.at(samp)>0?0:signalFlags.at(samp)+4)->Fill(deltaR(V4CPion2, V4NPion22), eventWeight);

        //if(V4Tau1Sol1.Pt() < 20 || V4Tau2Sol1.Pt() < 20 || V4Tau1Sol2.Pt() < 20 || V4Tau2Sol2.Pt() < 20) continue;
        //if(V4Tau1Sol1.Pt() > 90 || V4Tau2Sol1.Pt() > 90 || V4Tau1Sol2.Pt() > 90 || V4Tau2Sol2.Pt() > 90) continue;

        // Compute theta variable
        Double_t thetaSol1, thetaSol2;
        thetaSol1 = getTheta(V4CPion1, V4NPion1, V4CPion2, V4NPion2, V4Tau1Sol1, V4Tau2Sol1);
        thetaSol2 = getTheta(V4CPion1, V4NPion1, V4CPion2, V4NPion2, V4Tau1Sol2, V4Tau2Sol2);

        hThetaTemp.at(signalFlags.at(samp)+3)->Fill(thetaSol1, eventWeight/2.);
        hThetaTemp.at(signalFlags.at(samp)+3)->Fill(thetaSol2, eventWeight/2.);

        vector<double> vars1, vars2;
        vars1.push_back(thetaSol1);
        vars1.push_back(thetaSol2);
        //vars2.push_back(thetaSol2);

        finalEvents.push_back(vars1);
        //finalEvents.push_back(vars2);

        tempMassCuts += eventWeight;
        tempMassCutsError++;

        if(signalFlags.at(samp)<=0) h2Theta.at(signalFlags.at(samp)+3)->Fill(thetaSol1, thetaSol2, 5);

    } // end event loop

    for(Int_t i = 0; i < nDatasets; i++){
        for(Int_t x = 0; x < Nint(tempMassCuts); x++){
            Int_t randEntry = Nint(random->Uniform(finalEvents.size()-1));

            if(signalFlags.at(samp)>0){
                datasets.at(0).at(i).add(finalEvents.at(randEntry));
                //datasets.at(0).at(i).add(finalEvents.at(randEntry+1));
                datasets.at(1).at(i).add(finalEvents.at(randEntry));
                //datasets.at(1).at(i).add(finalEvents.at(randEntry+1));
                datasets.at(2).at(i).add(finalEvents.at(randEntry));
                //datasets.at(2).at(i).add(finalEvents.at(randEntry+1));
                datasets.at(3).at(i).add(finalEvents.at(randEntry));
                //datasets.at(3).at(i).add(finalEvents.at(randEntry+1));
            }
            else{
                datasets.at(signalFlags.at(samp)+3).at(i).add(finalEvents.at(randEntry));
                //datasets.at(signalFlags.at(samp)+3).at(i).add(finalEvents.at(randEntry+1));
            }
        }
    }

    selection.at(samp) += tempSelection;
    if(tempSelectionError > 0) selectionError.at(samp) += sqrtf(tempSelectionError)*tempSelection/tempSelectionError;
    massCuts.at(samp) += tempMassCuts;
    if(tempMassCutsError > 0) massCutsError.at(samp) += sqrtf(tempMassCutsError)*tempMassCuts/tempMassCutsError;

    TString out = TString::Format("%4.2f", tempMassCuts);
    out.Resize(9);

    cout << "\e[A";
    cout << inputfile << (signalFlags.at(samp)<=0 ? "\033[1;32m" : (out == "0       " ?  "\033[1;34m": "\033[1;31m")) << out << "\033[0m" << " events passed" << endl;

}

void saveResults(TString calcP, TString calcChi, TString calcL)
{
    cout << endl << endl;

    hThetaB->Add(Theta_jj);
    hThetaB->Add(Theta_dy);
    hThetaB->Add(Theta_ZZ);
    hThetaB->Add(Theta_Zll);
    Theta_sig0->Add(hThetaS0);
    Theta_sigpi4->Add(hThetaSpi4);
    Theta_sigpi2->Add(hThetaSpi2);
    Theta_sig3pi4->Add(hThetaS3pi4);

    Theta_obs0->Add(Theta_jj);
    Theta_obs0->Add(Theta_dy);
    Theta_obs0->Add(Theta_ZZ);
    Theta_obs0->Add(Theta_Zll);
    Theta_obs0->Add(hThetaS0);

    Theta_obspi4->Add(Theta_jj);
    Theta_obspi4->Add(Theta_dy);
    Theta_obspi4->Add(Theta_ZZ);
    Theta_obspi4->Add(Theta_Zll);
    Theta_obspi4->Add(hThetaSpi4);

    Theta_obspi2->Add(Theta_jj);
    Theta_obspi2->Add(Theta_dy);
    Theta_obspi2->Add(Theta_ZZ);
    Theta_obspi2->Add(Theta_Zll);
    Theta_obspi2->Add(hThetaSpi2);

    Theta_obs3pi4->Add(Theta_jj);
    Theta_obs3pi4->Add(Theta_dy);
    Theta_obs3pi4->Add(Theta_ZZ);
    Theta_obs3pi4->Add(Theta_Zll);
    Theta_obs3pi4->Add(hThetaS3pi4);

    TH1D *pS0andS0 = new TH1D("pS0andS0", "pS0andS0", 10, -0.005, 1.005);
    TH1D *pS0andSpi4 = new TH1D("pS0andSpi4", "pS0andSpi4", 10, -0.005, 1.005);
    TH1D *pS0andSpi2 = new TH1D("pS0andSpi2", "pS0andSpi2", 10, -0.005, 1.005);
    TH1D *pS0andS3pi4 = new TH1D("pS0andS3pi4", "pS0andS3pi4", 10, -0.005, 1.005);
    TH1D *chiS0andS0 = new TH1D("chiS0andS0", "chiS0andS0", 10, -0.005, 1.005);
    TH1D *chiS0andSpi4 = new TH1D("chiS0andSpi4", "chiS0andSpi4", 10, -0.005, 1.005);
    TH1D *chiS0andSpi2 = new TH1D("chiS0andSpi2", "chiS0andSpi2", 10, -0.005, 1.005);
    TH1D *chiS0andS3pi4 = new TH1D("chiS0andS3pi4", "chiS0andS3pi4", 10, -0.005, 1.005);
    TH1D *lS0andS0 = new TH1D("lS0andS0", "lS0andS0", 100, -0.0005, 0.005);
    TH1D *lS0andSpi4 = new TH1D("lS0andSpi4", "lS0andSpi4", 100, -0.0005, 0.005);
    TH1D *lS0andSpi2 = new TH1D("lS0andSpi2", "lS0andSpi2", 100, -0.0005, 0.005);
    TH1D *lS0andS3pi4 = new TH1D("lS0andS3pi4", "lS0andS3pi4", 100, -0.0005, 0.005);

    vector<TH1D*> hp, hl, hchi;
    hp.push_back(pS0andS0); hp.push_back(pS0andSpi4); hp.push_back(pS0andSpi2); hp.push_back(pS0andS3pi4);
    hchi.push_back(chiS0andS0); hchi.push_back(chiS0andSpi4); hchi.push_back(chiS0andSpi2); hchi.push_back(chiS0andS3pi4);
    hl.push_back(lS0andS0); hl.push_back(lS0andSpi4); hl.push_back(lS0andSpi2); hl.push_back(lS0andS3pi4);
    vector<TString> hpNames;
    hpNames.push_back("Delta=0 and Delta=0"); hpNames.push_back("Delta=0 and Delta=pi/4"); hpNames.push_back("Delta=0 and Delta=pi/2"); hpNames.push_back("Delta=0 and Delta=3pi/4");

    if(calcP == "true"){

        cout << "\033[1;34mMike's method\033[0m\n\n";

        TRandom3 *randGen = new TRandom3();
        randGen->SetSeed(0);
        for(Int_t i=0; i<4; i++){
            Double_t averageT=0, averageP=0, averagePError=0;

            for(Int_t x =0; x < nDatasets; x++){

                if(x!=0) cout << "\e[A" << "\e[A";
                cout << "[" << string((x+1)/(nDatasets/50),'#') << string(50-(x+1)/(nDatasets/50),' ') << "]  " << 100*(x+1)/nDatasets;
                cout << "\% completed.   Dataset " << x+1 << " of " << nDatasets << "  " << endl;

                Double_t tval;
                if(i == 0) tval = calcT(datasets.at(0).at(x),datasets.at(i).at( (x+1==nDatasets?0:x+1) ));
                else tval = calcT(datasets.at(0).at(x),datasets.at(i).at(x));

                Double_t ptval;
                Int_t pval=0;
                Int_t nperm=200;
                for(Int_t p=0; p<nperm; p++){
                    if(i == 0) ptval = permCalcT(datasets.at(0).at(x),datasets.at(i).at( (x+1==nDatasets?0:x+1) ), randGen);
                    else ptval = permCalcT(datasets.at(0).at(x),datasets.at(i).at(x), randGen);
                    //hp.at(i)->Fill(ptval);
                    if(ptval > tval) pval++;
                    if(p!=0) cout << "\e[A";
                    cout << "[" << string((p+1)/(nperm/50),'#') << string(50-(p+1)/(nperm/50),' ') << "]  " << 100*(p+1)/nperm << "\% completed.   Permutation " << p+1 << " of " << nperm;
                    cout << "  " << endl;
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
            cout << "Average p of " << hpNames.at(i) << ": " << averageP  << " +- " << averagePError << endl << endl << endl;

        }

    }

    TH1D *tempS01 = new TH1D("tempS01", "tempS01", 20, -3.1416, 3.1416);
    TH1D *tempS02 = new TH1D("tempS02", "tempS02", 20, -3.1416, 3.1416);
    TH1D *tempSpi4 = new TH1D("tempSpi4", "tempSpi4", 20, -3.1416, 3.1416);
    TH1D *tempSpi2 = new TH1D("tempSpi2", "tempSpi2", 20, -3.1416, 3.1416);
    TH1D *tempS3pi4 = new TH1D("tempS3pi4", "tempS3pi4", 20, -3.1416, 3.1416);

    if(calcChi == "true"){

        cout << "\033[1;34mChi squared test\033[0m\n\n";   

        for(Int_t x = 0; x < nDatasets; x++){        

            for(Int_t i = 0; i < datasets.at(0).at(x).size(); i++){
                tempS01->Fill(datasets.at(0).at(x).get(i,0), 0.5);
                tempS01->Fill(datasets.at(0).at(x).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(0).at( (x+1==nDatasets?0:x+1) ).size(); i++){
                tempS02->Fill(datasets.at(0).at( (x+1==nDatasets?0:x+1) ).get(i,0), 0.5);
                tempS02->Fill(datasets.at(0).at( (x+1==nDatasets?0:x+1) ).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(1).at(x).size(); i++){
                tempSpi4->Fill(datasets.at(1).at(x).get(i,0), 0.5);
                tempSpi4->Fill(datasets.at(1).at(x).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(2).at(x).size(); i++){
                tempSpi2->Fill(datasets.at(2).at(x).get(i,0), 0.5);
                tempSpi2->Fill(datasets.at(2).at(x).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(3).at(x).size(); i++){
                tempS3pi4->Fill(datasets.at(3).at(x).get(i,0), 0.5);
                tempS3pi4->Fill(datasets.at(3).at(x).get(i,1), 0.5);
            }

            Double_t chiS01andS02=0, chiS01andSpi4=0, chiS01andSpi2=0, chiS01andS3pi4=0;

            chiS01andS02 = getChiSquared(tempS01, tempS02);
            chiS01andSpi4 = getChiSquared(tempS01, tempSpi4);
            chiS01andSpi2 = getChiSquared(tempS01, tempSpi2);
            chiS01andS3pi4 = getChiSquared(tempS01, tempS3pi4);
            
            Double_t pvalS01andS02=0, pvalS01andSpi4=0, pvalS01andSpi2=0, pvalS01andS3pi4=0;

            pvalS01andS02 = TMath::Prob(chiS01andS02, 19);
            pvalS01andSpi4 = TMath::Prob(chiS01andSpi4, 19);
            pvalS01andSpi2 = TMath::Prob(chiS01andSpi2, 19);
            pvalS01andS3pi4 = TMath::Prob(chiS01andS3pi4, 19);

            chiS0andS0->Fill(pvalS01andS02);
            chiS0andSpi4->Fill(pvalS01andSpi4);
            chiS0andSpi2->Fill(pvalS01andSpi2);
            chiS0andS3pi4->Fill(pvalS01andS3pi4);

            tempS01->Reset("M");
            tempS02->Reset("M");
            tempSpi4->Reset("M");
            tempSpi2->Reset("M");
            tempS3pi4->Reset("M");

        }

    }

    if(calcL == "true"){

        cout << "\033[1;34mLikelihood\033[0m\n\n";

        for(Int_t x = 0; x < nDatasets; x++){        

            for(Int_t i = 0; i < datasets.at(0).at(x).size(); i++){
                tempS01->Fill(datasets.at(0).at(x).get(i,0), 0.5);
                tempS01->Fill(datasets.at(0).at(x).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(0).at( (x+1==nDatasets?0:x+1) ).size(); i++){
                tempS02->Fill(datasets.at(0).at( (x+1==nDatasets?0:x+1) ).get(i,0), 0.5);
                tempS02->Fill(datasets.at(0).at( (x+1==nDatasets?0:x+1) ).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(1).at(x).size(); i++){
                tempSpi4->Fill(datasets.at(1).at(x).get(i,0), 0.5);
                tempSpi4->Fill(datasets.at(1).at(x).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(2).at(x).size(); i++){
                tempSpi2->Fill(datasets.at(2).at(x).get(i,0), 0.5);
                tempSpi2->Fill(datasets.at(2).at(x).get(i,1), 0.5);
            }
            for(Int_t i = 0; i < datasets.at(3).at(x).size(); i++){
                tempS3pi4->Fill(datasets.at(3).at(x).get(i,0), 0.5);
                tempS3pi4->Fill(datasets.at(3).at(x).get(i,1), 0.5);
            }

            Double_t likeS01andS02=0, likeS01andSpi4=0, likeS01andSpi2=0, likeS01andS3pi4=0;

            likeS01andS02 = getLikelihood(tempS01, tempS02);
            likeS01andSpi4 = getLikelihood(tempS01, tempSpi4);
            likeS01andSpi2 = getLikelihood(tempS01, tempSpi2);
            likeS01andS3pi4 = getLikelihood(tempS01, tempS3pi4);

            lS0andS0->Fill(likeS01andS02, 2);
            lS0andSpi4->Fill(likeS01andSpi4, 2);
            lS0andSpi2->Fill(likeS01andSpi2, 2);
            lS0andS3pi4->Fill(likeS01andS3pi4, 2);

            //cout << likeS01andS02 << " " << likeS01andSpi4 << " " << likeS01andSpi2 << " " << likeS01andS3pi4 << endl;

            tempS01->Reset("M");
            tempS02->Reset("M");
            tempSpi4->Reset("M");
            tempSpi2->Reset("M");
            tempS3pi4->Reset("M");

        }

        for(Int_t i=0; i<4; i++){
            Double_t Likelihood=1.;

            for(Int_t x = 1; x <= 20; x++){

                Likelihood*=Poisson(Nint(hThetaObs.at(0)->GetBinContent(x)), Nint(hThetaObs.at(i)->GetBinContent(x)));
                Likelihood/=Poisson(Nint(hThetaObs.at(0)->GetBinContent(x)), Nint(hThetaObs.at(0)->GetBinContent(x)));
            }

            //cout << "Likelihood of " << hpNames.at(i) << " is " << Likelihood << endl;
        }

    }

    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);

    gStyle->SetOptStat(kFALSE);

    histogram(hTheta, histogramNames, c1, "#Theta variable", "Normalized yield", "Theta");

    TF1 *myfit = new TF1("myfit","[0] + [1]*cos(x + [2])", 0, 2.);

    myfit->SetParName(0,"C");
    myfit->SetParName(1,"A");
    myfit->SetParName(2,"W");
    myfit->SetParameter(0, 0.5);
    myfit->SetParameter(1, 0.1);
    myfit->SetParameter(2, 0);

    hTheta.at(1)->Fit("myfit");
    hTheta.at(2)->Fit("myfit");
    hTheta.at(3)->Fit("myfit");
    hTheta.at(4)->Fit("myfit");

    //histogram(hPhi, histogramNames, c1, "#varphi^{*} variable", "Normalized yield", "Phi");
    histogram(genhistos, histogramNames, c1, "H mass", "Normalized yield", "Zmass");
    //histogram(genhistos2, histogramNames2, c1, "Z mass", "Fraction", "Zmass");
    histogram(hTheta, histogramNames, c1, "#Theta variable", "Normalized yield", "ThetaFit");
    histogram(h2Theta, "Theta corr", c1, "ThetaCorr");

    if(calcP == "true" && calcChi == "true"){
        for(Int_t i=0; i<4; i++){
            TString filename="p_"; filename+=i;
            histogram(hp.at(i), hchi.at(i), hpNames.at(i), c1, "p values", "Toys", filename);
        }
    }
    else if(calcP == "true"){
        for(Int_t i=0; i<4; i++){
            TString filename="p_"; filename+=i;
            histogram(hp.at(i), hpNames.at(i), c1, "p distribution (Mike)", "Fraction", filename);
        }
    }
    else if(calcChi == "true"){
        for(Int_t i=0; i<4; i++){
            TString filename="p_"; filename+=i;
            histogram(hl.at(i), hpNames.at(i), c1, "p distribution (Chi squared)", "Fraction", filename);
        }
    }

    if(calcL == "true"){
        for(Int_t i=0; i<4; i++){
            TString filename="l_"; filename+=i;
            histogram(hl.at(i), hpNames.at(i), c1, "Likelihood values", "Toys", filename);
        }
    }

    TFile f("Histograms.root","RECREATE");
    Theta_obs0->Write();
    Theta_obspi4->Write();
    Theta_obspi2->Write();
    Theta_obs3pi4->Write();
    Theta_sig0->Write();
    Theta_sigpi4->Write();
    Theta_sigpi2->Write();
    Theta_sig3pi4->Write();
    Theta_jj->Write();
    Theta_dy->Write();
    Theta_ZZ->Write();
    Theta_Zll->Write();
    pS0andS0->Write();
    pS0andSpi4->Write();
    pS0andSpi2->Write();
    pS0andS3pi4->Write();
    f.Close();

    cout << "\n\n\nProcess finished\nPrinting results...\n\n" << "\033[1;34mResults\033[0m\n\n";
    
    for(UInt_t x = 0; x<sampleNames.size();x++){
        
        cout << (signalFlags.at(x)<=0 ? "\033[1;32m" : "\033[1;31m") << sampleNames.at(x) << "\033[0m\n";
        cout << "Events after selection: " << selection.at(x) << " +- " << selectionError.at(x) << endl << endl;
        cout << "Events after mass cuts: " << massCuts.at(x) << " +- " << massCutsError.at(x) << endl << endl;
        
    }

    cout << "Done\nExiting...\n\n\n";

}

Double_t getChiSquared(TH1D *histo1, TH1D *histo2){

    Double_t chiSquared=0;

    for(Int_t i = 1; i <= 20; i++){
        chiSquared += Power(histo1->Integral()*histo2->GetBinContent(i) - histo2->Integral()*histo1->GetBinContent(i) , 2)/(histo1->Integral()*histo2->Integral()*(histo1->GetBinContent(i) + histo2->GetBinContent(i)));
    }

    return chiSquared;
}

Double_t getLikelihood(TH1D *histo1, TH1D *histo2){

    Double_t Likelihood=1.;

    for(Int_t x = 1; x <= 20; x++){

        Likelihood*=Poisson(Nint(histo1->GetBinContent(x)), Nint(histo2->GetBinContent(x)));
        Likelihood/=Poisson(Nint(histo1->GetBinContent(x)), Nint(histo1->GetBinContent(x)));
    }

    return Likelihood;
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

    can->SaveAs(name + ".jpg");
}

/*
 * FUNCTION FOR SAVING TWO (P-value) HISTOGRAMs
 */

void histogram(TH1D *histo1, TH1D *histo2, const TString histName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){
    
    Double_t max = histo1->GetMaximum();
    if(histo2->GetMaximum() > max) max = histo2->GetMaximum();

    histo1->SetMaximum(max*1.2);
    histo2->SetMaximum(max*1.2);
    histo1->SetMinimum(0);
    histo2->SetMinimum(0);

    histo1->SetLineWidth(3);
    histo2->SetLineWidth(3);
    histo2->SetLineColor(kRed);
    histo1->Draw();
    histo2->Draw("same");
    // add axis labels
    histo1->GetXaxis()->SetTitle(xTitle);
    histo1->GetYaxis()->SetTitle(yTitle);
    histo1->SetTitle(histName); // title on top

    TLegend *leg = new TLegend(0.605,0.675,0.885,0.875);
    leg->SetTextFont(72);
    leg->SetTextSize(0.04);
    leg->AddEntry(histo1, "Mike's method","l");
    leg->AddEntry(histo2, "Chi-squared test","l");
    leg->Draw();

    can->SaveAs(name + ".jpg");
}

/*
 * FUNCTION FOR SAVING ONE 2D HISTOGRAM
 */

void histogram(TH2D *histo, const TString histName, TCanvas *can, const TString name){
 
    histo->SetLineWidth(3);
    histo->Draw();
    // add axis labels
    histo->SetTitle(histName); // title on top

    can->SaveAs(name + ".jpg");
}

/*
 * FUNCTION FOR SAVING MULTIPLE 2D HISTOGRAM
 */

void histogram(vector<TH2D*> histos, const TString histName, TCanvas *can, const TString name){
 
    can->Clear();
    can->Divide(2,2);
    for(Int_t x=0; x<=3; x++){
        can->cd(x+1);
        histos.at(x)->Draw();
        //histos.at(x)->SetTitle(TString::Format(histName + " %i", x));
        can->cd();
    }

    can->SaveAs(name + ".jpg");
}

/*
 * FUNCTION FOR SAVING MULTIPLE HISTOGRAMS
 */

void histogram(vector<TH1D*> histos, vector<TString> histNames, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){   
    
    if(histos.size()!=histNames.size()){ cout << "Number of histograms and names don't match." << endl; return;}

    Double_t max=0, min=1;

    for(Int_t i=0; i<histos.size(); i++){
        Double_t integral = histos.at(i)->Integral();
        if(integral != 0.){
            histos.at(i)->Scale(Double_t(1)/integral);
            if(histos.at(i)->GetMaximum()>max) max=histos.at(i)->GetMaximum();
            if(histos.at(i)->GetMinimum()<min) min=histos.at(i)->GetMinimum();
        }
    }

    max*=1.1;
    min*=0.9;
    //if(name == "Theta"){max=.1; min=0;}

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
