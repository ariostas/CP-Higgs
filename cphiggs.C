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
TH1D *hThetaS0 = new TH1D("hThetaS1", "hThetaS1", 10, -3.1416, 3.1416);
TH1D *hThetaSpi4 = new TH1D("hThetaS2", "hThetaS2", 10, -3.1416, 3.1416);
TH1D *hThetaSpi2 = new TH1D("hThetaS3", "hThetaS3", 10, -3.1416, 3.1416);
TH1D *hThetaS3pi4 = new TH1D("hThetaS4", "hThetaS4", 10, -3.1416, 3.1416);
TH1D *hThetaB = new TH1D("hThetaB", "hThetaB", 10, -3.1416, 3.1416);

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
Dataset data0(2,10000), datapi4(2,10000), datapi2(2,10000), data3pi4(2,10000);
vector<Dataset> datasets;

// Initialyze storage variables
vector<Double_t> total, selection, kinematicCuts, massCuts;
vector<Double_t> totalError, selectionError, kinematicCutsError, massCutsError;
vector<Int_t> signalFlags;
vector<TString> sampleNames;
vector<TH1D*> hTheta, genhistos, hPhi, genhistos2;
vector<TString> histogramNames, histogramNames2;

/*
 * MAIN FUNCTION
 */

 void cphiggs(TString sample = "all", TString inputFile = "xsec.txt"){
    
    cout << "\n\nStarting process...\n\n";

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

    datasets.push_back(data0);
    datasets.push_back(datapi4);
    datasets.push_back(datapi2);
    datasets.push_back(data3pi4);
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}
    
    Int_t nSamples = 0;
    TString sampleName, sampleBin, crossSection, signalFlag;
    
    while(ifs >> sampleName >> sampleBin >> crossSection >> signalFlag){
        
        if(sampleName == "#") continue;
        
        if(sample != sampleName && sample != "all") continue;
        
        Int_t matchedSample = sampleNames.size();
        for(UInt_t x=0;x<sampleNames.size();x++){
            
            if(sampleName == sampleNames.at(x)) matchedSample = x;
            
        }
        
        nSamples = matchedSample;
        
        if(nSamples == sampleNames.size()){
            
            sampleNames.push_back(sampleName);
            signalFlags.push_back(atof(string(signalFlag).c_str()));
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

    cout << inputfile;

    vector<vector<Double_t> > finalEvents;

    // Set up storage variables
    Float_t genTau1_pt, genTau2_pt;
    Float_t genTau1_eta, genTau2_eta;
    Float_t genTau1_phi, genTau2_phi;
    Float_t genTau1_mass, genTau2_mass;

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

    // Particles from Z
    Float_t z1_pt, z2_pt;
    Float_t z1_eta, z2_eta;
    Float_t z1_phi, z2_phi;
    Float_t z1_mass, z2_mass;

    // Jets from Z
    Float_t jetz1_pt, jetz2_pt;
    Float_t jetz1_eta, jetz2_eta;
    Float_t jetz1_phi, jetz2_phi;
    Float_t jetz1_mass, jetz2_mass;

    // Z
    Float_t z_pt, z_eta, z_phi, z_mass;

    Float_t eventWeight;

    TLorentzVector vTau1Sol1, vTau1Sol2, vTau2Sol1, vTau2Sol2, vcpion1, vcpion2, vnpion1, vnpion2, vHiggs, vrho1, vrho2;
    UInt_t nProngTau1=0, nProngTau2=0, nLeptons=0, zToLep=0;

    Double_t thetaSol1=0, thetaSol2=0;

    TFile* infile = new TFile(inputFileTemp);
    assert(infile);
    TTree* intree = (TTree*) infile->Get("Events");
    assert(intree);

    intree->SetBranchAddress("nLeptons",		&nLeptons);
    intree->SetBranchAddress("eventWeight",		&eventWeight);
    intree->SetBranchAddress("nProngTau1",		&nProngTau1);
    intree->SetBranchAddress("nProngTau2",		&nProngTau2);
    intree->SetBranchAddress("zToLep",          &zToLep);
    intree->SetBranchAddress("genTau1_pt",		&genTau1_pt);
    intree->SetBranchAddress("genTau1_eta",		&genTau1_eta);
    intree->SetBranchAddress("genTau1_phi",		&genTau1_phi);
    intree->SetBranchAddress("genTau1_mass",	&genTau1_mass);
    intree->SetBranchAddress("genTau2_pt",		&genTau2_pt);
    intree->SetBranchAddress("genTau2_eta",		&genTau2_eta);
    intree->SetBranchAddress("genTau2_phi",		&genTau2_phi);
    intree->SetBranchAddress("genTau2_mass",	&genTau2_mass);
    intree->SetBranchAddress("pions1_pt",		&pions1_pt);
    intree->SetBranchAddress("pions1_eta",		&pions1_eta);
    intree->SetBranchAddress("pions1_phi",		&pions1_phi);
    intree->SetBranchAddress("pions1_mass",		&pions1_mass);
    intree->SetBranchAddress("pions2_pt",		&pions2_pt);
    intree->SetBranchAddress("pions2_eta",		&pions2_eta);
    intree->SetBranchAddress("pions2_phi",		&pions2_phi);
    intree->SetBranchAddress("pions2_mass",		&pions2_mass);
    intree->SetBranchAddress("neutpions1_pt",	&neutpions1_pt);
    intree->SetBranchAddress("neutpions1_eta",	&neutpions1_eta);
    intree->SetBranchAddress("neutpions1_phi",	&neutpions1_phi);
    intree->SetBranchAddress("neutpions1_mass",	&neutpions1_mass);
    intree->SetBranchAddress("neutpions2_pt",	&neutpions2_pt);
    intree->SetBranchAddress("neutpions2_eta",	&neutpions2_eta);
    intree->SetBranchAddress("neutpions2_phi",	&neutpions2_phi);
    intree->SetBranchAddress("neutpions2_mass",	&neutpions2_mass);
    intree->SetBranchAddress("z1_pt",           &z1_pt);
    intree->SetBranchAddress("z1_eta",          &z1_eta);
    intree->SetBranchAddress("z1_phi",          &z1_phi);
    intree->SetBranchAddress("z1_mass",         &z1_mass);
    intree->SetBranchAddress("z2_pt",           &z2_pt);
    intree->SetBranchAddress("z2_eta",          &z2_eta);
    intree->SetBranchAddress("z2_phi",          &z2_phi);
    intree->SetBranchAddress("z2_mass",         &z2_mass);
    intree->SetBranchAddress("jetz1_pt",        &jetz1_pt);
    intree->SetBranchAddress("jetz1_eta",       &jetz1_eta);
    intree->SetBranchAddress("jetz1_phi",       &jetz1_phi);
    intree->SetBranchAddress("jetz1_mass",      &jetz1_mass);
    intree->SetBranchAddress("jetz2_pt",        &jetz2_pt);
    intree->SetBranchAddress("jetz2_eta",       &jetz2_eta);
    intree->SetBranchAddress("jetz2_phi",       &jetz2_phi);
    intree->SetBranchAddress("jetz2_mass",      &jetz2_mass);
    intree->SetBranchAddress("z_pt",            &z_pt);
    intree->SetBranchAddress("z_eta",           &z_eta);
    intree->SetBranchAddress("z_phi",           &z_phi);
    intree->SetBranchAddress("z_mass",          &z_mass);

    Double_t tempSelection=0, tempSelectionError=0;

    for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++)   // Event loop
    {
        intree->GetEntry(iEntry);

        if(pions1_pt==0 || pions2_pt==0 || neutpions1_pt==0 || neutpions2_pt==0 || nProngTau1!=1 || nProngTau2!=1) continue;

        // Charged pions 4-vectors
        vcpion1.SetPtEtaPhiM(pions1_pt, pions1_eta, pions1_phi, pions1_mass);
        vcpion2.SetPtEtaPhiM(pions2_pt, pions2_eta, pions2_phi, pions2_mass);

        // Neutral pions 4-vectors
        vnpion1.SetPtEtaPhiM(neutpions1_pt, neutpions1_eta, neutpions1_phi, neutpions1_mass);
        vnpion2.SetPtEtaPhiM(neutpions2_pt, neutpions2_eta, neutpions2_phi, neutpions2_mass);

        // Rho 4-vectors
        vrho1=vcpion1+vnpion1;
        vrho2=vcpion2+vnpion2;

        TLorentzVector v1, v2, vZ;

        if(zToLep){
            v1.SetPtEtaPhiM(z1_pt, z1_eta, z1_phi, z1_mass);
            v2.SetPtEtaPhiM(z2_pt, z2_eta, z2_phi, z2_mass);
        }
        else{
            v1.SetPtEtaPhiM(jetz1_pt, jetz1_eta, jetz1_phi, jetz1_mass);
            v2.SetPtEtaPhiM(jetz2_pt, jetz2_eta, jetz2_phi, jetz2_mass);
        }

        TLorentzVector vInit;
        vInit.SetPxPyPzE(0,0,0,240);
        vZ = v1 + v2;

        // If the event has a generated Z use its real 4-vector
        if(z_pt>0.01) vZ.SetPtEtaPhiM(z_pt, z_eta, z_phi, z_mass);

        vHiggs = vInit-vZ;

        // Scale down the backgrounds to simulate a cut on the mass of the Higgs.
        if(signalFlags.at(samp)==0) eventWeight/=2.0;

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
        Double_t phiSol1, phiSol2;
        phiSol1 = getPhi(vcpion1, vnpion1, vcpion2, vnpion2, vTau1Sol1, vTau2Sol1);
        phiSol2 = getPhi(vcpion1, vnpion1, vcpion2, vnpion2, vTau1Sol2, vTau2Sol2);

        // Check whether theta variable is well defined
        if(thetaSol1 != thetaSol1 || thetaSol2 != thetaSol2) continue;

        // Fill histogram for theta variable
        hTheta.at(signalFlags.at(samp))->Fill(thetaSol1);
        hTheta.at(signalFlags.at(samp))->Fill(thetaSol2);

        // Check whether phi variable is well defined
        if(phiSol1 != phiSol1 || phiSol2 != phiSol2) continue;

        // Fill histogram for phi variable
        hPhi.at(signalFlags.at(samp))->Fill(phiSol1);
        hPhi.at(signalFlags.at(samp))->Fill(phiSol2);

        genhistos.at(signalFlags.at(samp))->Fill(vNeutrino1Sol1.Pt());

        vector<double> vars1, vars2;
        vars1.push_back(thetaSol1);
        vars1.push_back(eventWeight/10.0);
        vars2.push_back(thetaSol2);
        vars2.push_back(eventWeight/10.0);

        finalEvents.push_back(vars1);
        finalEvents.push_back(vars2);

        if(signalFlags.at(samp)==0 && iEntry < intree->GetEntries()/10.0){
            datasets.at(0).add(vars1);
            datasets.at(0).add(vars2);
            datasets.at(1).add(vars1);
            datasets.at(1).add(vars2);
            datasets.at(2).add(vars1);
            datasets.at(2).add(vars2);
            datasets.at(3).add(vars1);
            datasets.at(3).add(vars2);
        }
        else if(iEntry < intree->GetEntries()/10.0){
            datasets.at(signalFlags.at(samp)-1).add(vars1);
            datasets.at(signalFlags.at(samp)-1).add(vars2);
        }

        tempSelection+=eventWeight/Double_t(10);
        tempSelectionError++;

    } // end event loop

    TString out = "";
    out += tempSelection;
    out.Resize(9);

    selection.at(samp) += tempSelection;
    if(tempSelectionError > 0) selectionError.at(samp) += sqrtf(tempSelectionError)*tempSelection/tempSelectionError;

    cout << (signalFlags.at(samp) ? "\033[1;32m" : (out == "0       " ?  "\033[1;34m": "\033[1;31m")) << out << "\033[0m" << " events passed all cuts" << endl;

    infile->Close();
}

/*
 * FUNCTION FOR PRINTING AND SAVING THE RESULTS
 */

void saveResults()
{
    cout << endl;

    TH1D *hp1 = new TH1D("hp1", "hp1", 50, -0.002, 0.004);
    TH1D *hp2 = new TH1D("hp2", "hp2", 50, -0.002, 0.004);
    TH1D *hp3 = new TH1D("hp3", "hp3", 50, -0.002, 0.004);
    TH1D *hp4 = new TH1D("hp4", "hp4", 50, -0.002, 0.004);

    vector<TH1D*> hp;
    hp.push_back(hp1); hp.push_back(hp2); hp.push_back(hp3); hp.push_back(hp4);
    vector<TString> hpNames;
    hpNames.push_back("Delta=0 and Delta=0"); hpNames.push_back("Delta=0 and Delta=pi/4"); hpNames.push_back("Delta=0 and Delta=pi/2"); hpNames.push_back("Delta=0 and Delta=3pi/4");

    for(Int_t i=0; i<4; i++){
        Double_t tval = calcTW(datasets.at(0),datasets.at(i));

        Double_t ptval;
        srand(22);
        Int_t pval=0;
        Int_t nperm=100;
        for(Int_t p=0; p<nperm; p++){
            ptval = permCalcTW(datasets.at(0),datasets.at(i));
            hp.at(i)->Fill(ptval);
            if(ptval > tval) pval++;
            if(p!=0) cout << "\e[A";
            cout << "[" << string((p+1)/(nperm/50),'-') << string(50-(p+1)/(nperm/50),' ') << "]  " << 100*(p+1)/nperm << "\% completed.   Permutation " << p+1 << " of " << nperm << endl;
        }
        Double_t p = pval/Double_t(nperm), pError = Sqrt(p*(1-p)/Double_t(nperm));

        cout << "T of " << hpNames.at(i) << ": " << tval << endl;
        cout << "p of " << hpNames.at(i) << ": " << p  << " +- " << pError << endl << endl;

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

    cout << "\n\n\nProcess finished\nPrinting results...\n\n" << "\033[1;34mResults\033[0m\n\n";
    
    for(UInt_t x = 0; x<sampleNames.size();x++){
        
        cout << (signalFlags.at(x) ? "\033[1;32m" : "\033[1;31m") << sampleNames.at(x) << "\033[0m\n";
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
        histos.at(i)->Scale(Double_t(1)/histos.at(i)->Integral());
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

/* 
 * FUNCTION FOR GETTING THE 4-MOMENTUM OF THE FIRST NEUTRINO
 */

TLorentzVector getNeut1(TLorentzVector vZ, TLorentzVector vInit, TLorentzVector vRho1, TLorentzVector vRho2, Int_t sol){

    TLorentzVector vHiggs = vInit - vZ;

    TVector3 v3Higgs;
    v3Higgs.SetXYZ(vHiggs.Px()/vHiggs.E(), vHiggs.Py()/vHiggs.E(), vHiggs.Pz()/vHiggs.E());

    // Boost vectors to Higgs frame
    vHiggs.Boost(-v3Higgs);
    vRho1.Boost(-v3Higgs);
    vRho2.Boost(-v3Higgs);

    Double_t Etau = vHiggs.E()/2, Mtau = 1.77682;

    //cout << vRho1.E() << " " << vRho1.X() << " " << vRho1.Y() << " " << vRho1.Z() << endl;

    Double_t Prho1x = vRho1.X(), Prho1y = vRho1.Y(), Prho1z = vRho1.Z();
    Double_t Prho2x = vRho2.X(), Prho2y = vRho2.Y(), Prho2z = vRho2.Z();

    Double_t Enu1 = Etau - vRho1.E(), Enu2 = Etau - vRho2.E();

    Double_t tempVal = Power(-4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x + 4*Power(Etau,2)*Power(Prho1y,2)*Prho2x - 
            4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x - 4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x + 
            4*Power(Etau,2)*Power(Prho1z,2)*Prho2x - 4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x + 
            8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) + 8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) + 
            4*Power(Prho1y,2)*Power(Prho2x,3) + 4*Power(Prho1z,2)*Power(Prho2x,3) + 
            4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y - 4*Power(Etau,2)*Prho1x*Prho1y*Prho2y + 
            4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y + 
            4*Power(Etau,2)*Prho1y*Prho2x*Prho2y - 4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y - 
            12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y + 4*Power(Prho1y,3)*Prho2x*Prho2y + 
            4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y - 4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y + 
            4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) - 4*Power(Etau,2)*Prho1x*Power(Prho2y,2) + 
            4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) + 4*Power(Prho1x,3)*Power(Prho2y,2) - 
            4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) + 
            4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) + 4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) - 
            4*Prho1x*Prho1y*Power(Prho2y,3) + 4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z - 
            4*Power(Etau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z - 
            4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z + 4*Power(Etau,2)*Prho1z*Prho2x*Prho2z - 
            4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z - 12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z + 
            4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z + 4*Power(Prho1z,3)*Prho2x*Prho2z - 
            4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z - 16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z - 
            4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z + 4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) - 
            4*Power(Etau,2)*Prho1x*Power(Prho2z,2) + 4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) + 
            4*Power(Prho1x,3)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) - 
            4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) + 4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) + 
            4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) - 4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) - 
            4*Prho1x*Prho1z*Power(Prho2z,3),2) - 
          4*(4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
             8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
             4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
             8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
             4*Power(Prho1y,2)*Power(Prho2z,2))*
           (Power(Enu2,4)*Power(Prho1y,2) - 2*Power(Enu2,2)*Power(Etau,2)*Power(Prho1y,2) + 
             Power(Etau,4)*Power(Prho1y,2) + 2*Power(Enu2,2)*Power(Mtau,2)*Power(Prho1y,2) - 
             2*Power(Etau,2)*Power(Mtau,2)*Power(Prho1y,2) + Power(Mtau,4)*Power(Prho1y,2) + 
             Power(Enu2,4)*Power(Prho1z,2) - 2*Power(Enu2,2)*Power(Etau,2)*Power(Prho1z,2) + 
             Power(Etau,4)*Power(Prho1z,2) + 2*Power(Enu2,2)*Power(Mtau,2)*Power(Prho1z,2) - 
             2*Power(Etau,2)*Power(Mtau,2)*Power(Prho1z,2) + Power(Mtau,4)*Power(Prho1z,2) - 
             4*Power(Enu2,2)*Prho1x*Power(Prho1y,2)*Prho2x + 
             4*Power(Etau,2)*Prho1x*Power(Prho1y,2)*Prho2x - 
             4*Power(Mtau,2)*Prho1x*Power(Prho1y,2)*Prho2x - 
             4*Power(Enu2,2)*Prho1x*Power(Prho1z,2)*Prho2x + 
             4*Power(Etau,2)*Prho1x*Power(Prho1z,2)*Prho2x - 
             4*Power(Mtau,2)*Prho1x*Power(Prho1z,2)*Prho2x - 
             2*Power(Enu2,2)*Power(Prho1y,2)*Power(Prho2x,2) + 
             2*Power(Etau,2)*Power(Prho1y,2)*Power(Prho2x,2) - 
             2*Power(Mtau,2)*Power(Prho1y,2)*Power(Prho2x,2) + 
             4*Power(Prho1x,2)*Power(Prho1y,2)*Power(Prho2x,2) - 
             2*Power(Enu2,2)*Power(Prho1z,2)*Power(Prho2x,2) + 
             2*Power(Etau,2)*Power(Prho1z,2)*Power(Prho2x,2) - 
             2*Power(Mtau,2)*Power(Prho1z,2)*Power(Prho2x,2) + 
             4*Power(Prho1x,2)*Power(Prho1z,2)*Power(Prho2x,2) + 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2x,3) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2x,3) + 
             Power(Prho1y,2)*Power(Prho2x,4) + Power(Prho1z,2)*Power(Prho2x,4) + 
             2*Power(Enu1,2)*Power(Enu2,2)*Prho1y*Prho2y - 
             2*Power(Enu1,2)*Power(Etau,2)*Prho1y*Prho2y - 
             2*Power(Enu2,2)*Power(Etau,2)*Prho1y*Prho2y + 2*Power(Etau,4)*Prho1y*Prho2y + 
             2*Power(Enu1,2)*Power(Mtau,2)*Prho1y*Prho2y + 
             2*Power(Enu2,2)*Power(Mtau,2)*Prho1y*Prho2y - 
             4*Power(Etau,2)*Power(Mtau,2)*Prho1y*Prho2y + 2*Power(Mtau,4)*Prho1y*Prho2y + 
             2*Power(Enu2,2)*Power(Prho1x,2)*Prho1y*Prho2y - 
             2*Power(Etau,2)*Power(Prho1x,2)*Prho1y*Prho2y + 
             2*Power(Mtau,2)*Power(Prho1x,2)*Prho1y*Prho2y - 2*Power(Enu2,2)*Power(Prho1y,3)*Prho2y + 
             2*Power(Etau,2)*Power(Prho1y,3)*Prho2y - 2*Power(Mtau,2)*Power(Prho1y,3)*Prho2y - 
             2*Power(Enu2,2)*Prho1y*Power(Prho1z,2)*Prho2y + 
             2*Power(Etau,2)*Prho1y*Power(Prho1z,2)*Prho2y - 
             2*Power(Mtau,2)*Prho1y*Power(Prho1z,2)*Prho2y - 
             4*Power(Enu1,2)*Prho1x*Prho1y*Prho2x*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2x*Prho2y - 
             4*Power(Mtau,2)*Prho1x*Prho1y*Prho2x*Prho2y - 4*Power(Prho1x,3)*Prho1y*Prho2x*Prho2y + 
             4*Prho1x*Power(Prho1y,3)*Prho2x*Prho2y + 4*Prho1x*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y - 
             2*Power(Enu1,2)*Prho1y*Power(Prho2x,2)*Prho2y + 
             2*Power(Etau,2)*Prho1y*Power(Prho2x,2)*Prho2y - 
             2*Power(Mtau,2)*Prho1y*Power(Prho2x,2)*Prho2y - 
             2*Power(Prho1x,2)*Prho1y*Power(Prho2x,2)*Prho2y + 
             2*Power(Prho1y,3)*Power(Prho2x,2)*Prho2y + 
             2*Prho1y*Power(Prho1z,2)*Power(Prho2x,2)*Prho2y + Power(Enu1,4)*Power(Prho2y,2) - 
             2*Power(Enu1,2)*Power(Etau,2)*Power(Prho2y,2) + Power(Etau,4)*Power(Prho2y,2) + 
             2*Power(Enu1,2)*Power(Mtau,2)*Power(Prho2y,2) - 
             2*Power(Etau,2)*Power(Mtau,2)*Power(Prho2y,2) + Power(Mtau,4)*Power(Prho2y,2) + 
             2*Power(Enu1,2)*Power(Prho1x,2)*Power(Prho2y,2) - 
             2*Power(Etau,2)*Power(Prho1x,2)*Power(Prho2y,2) + 
             2*Power(Mtau,2)*Power(Prho1x,2)*Power(Prho2y,2) + Power(Prho1x,4)*Power(Prho2y,2) - 
             2*Power(Enu1,2)*Power(Prho1y,2)*Power(Prho2y,2) - 
             2*Power(Enu2,2)*Power(Prho1y,2)*Power(Prho2y,2) + 
             4*Power(Etau,2)*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Power(Mtau,2)*Power(Prho1y,2)*Power(Prho2y,2) - 
             2*Power(Prho1x,2)*Power(Prho1y,2)*Power(Prho2y,2) + Power(Prho1y,4)*Power(Prho2y,2) - 
             2*Power(Enu1,2)*Power(Prho1z,2)*Power(Prho2y,2) - 
             2*Power(Enu2,2)*Power(Prho1z,2)*Power(Prho2y,2) + 
             2*Power(Prho1x,2)*Power(Prho1z,2)*Power(Prho2y,2) + 
             2*Power(Prho1y,2)*Power(Prho1z,2)*Power(Prho2y,2) + Power(Prho1z,4)*Power(Prho2y,2) + 
             4*Prho1x*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) + 
             4*Prho1x*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 
             2*Power(Prho1y,2)*Power(Prho2x,2)*Power(Prho2y,2) + 
             2*Power(Prho1z,2)*Power(Prho2x,2)*Power(Prho2y,2) - 
             2*Power(Enu1,2)*Prho1y*Power(Prho2y,3) + 2*Power(Etau,2)*Prho1y*Power(Prho2y,3) - 
             2*Power(Mtau,2)*Prho1y*Power(Prho2y,3) - 2*Power(Prho1x,2)*Prho1y*Power(Prho2y,3) + 
             2*Power(Prho1y,3)*Power(Prho2y,3) + 2*Prho1y*Power(Prho1z,2)*Power(Prho2y,3) + 
             Power(Prho1y,2)*Power(Prho2y,4) + Power(Prho1z,2)*Power(Prho2y,4) + 
             2*Power(Enu1,2)*Power(Enu2,2)*Prho1z*Prho2z - 
             2*Power(Enu1,2)*Power(Etau,2)*Prho1z*Prho2z - 
             2*Power(Enu2,2)*Power(Etau,2)*Prho1z*Prho2z + 2*Power(Etau,4)*Prho1z*Prho2z + 
             2*Power(Enu1,2)*Power(Mtau,2)*Prho1z*Prho2z + 
             2*Power(Enu2,2)*Power(Mtau,2)*Prho1z*Prho2z - 
             4*Power(Etau,2)*Power(Mtau,2)*Prho1z*Prho2z + 2*Power(Mtau,4)*Prho1z*Prho2z + 
             2*Power(Enu2,2)*Power(Prho1x,2)*Prho1z*Prho2z - 
             2*Power(Etau,2)*Power(Prho1x,2)*Prho1z*Prho2z + 
             2*Power(Mtau,2)*Power(Prho1x,2)*Prho1z*Prho2z - 
             2*Power(Enu2,2)*Power(Prho1y,2)*Prho1z*Prho2z + 
             2*Power(Etau,2)*Power(Prho1y,2)*Prho1z*Prho2z - 
             2*Power(Mtau,2)*Power(Prho1y,2)*Prho1z*Prho2z - 2*Power(Enu2,2)*Power(Prho1z,3)*Prho2z + 
             2*Power(Etau,2)*Power(Prho1z,3)*Prho2z - 2*Power(Mtau,2)*Power(Prho1z,3)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Prho1z*Prho2x*Prho2z + 
             4*Power(Etau,2)*Prho1x*Prho1z*Prho2x*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2x*Prho2z - 4*Power(Prho1x,3)*Prho1z*Prho2x*Prho2z + 
             4*Prho1x*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z + 4*Prho1x*Power(Prho1z,3)*Prho2x*Prho2z - 
             2*Power(Enu1,2)*Prho1z*Power(Prho2x,2)*Prho2z + 
             2*Power(Etau,2)*Prho1z*Power(Prho2x,2)*Prho2z - 
             2*Power(Mtau,2)*Prho1z*Power(Prho2x,2)*Prho2z - 
             2*Power(Prho1x,2)*Prho1z*Power(Prho2x,2)*Prho2z + 
             2*Power(Prho1y,2)*Prho1z*Power(Prho2x,2)*Prho2z + 
             2*Power(Prho1z,3)*Power(Prho2x,2)*Prho2z + 8*Power(Etau,2)*Prho1y*Prho1z*Prho2y*Prho2z - 
             8*Power(Mtau,2)*Prho1y*Prho1z*Prho2y*Prho2z - 
             8*Power(Prho1x,2)*Prho1y*Prho1z*Prho2y*Prho2z - 
             2*Power(Enu1,2)*Prho1z*Power(Prho2y,2)*Prho2z + 
             2*Power(Etau,2)*Prho1z*Power(Prho2y,2)*Prho2z - 
             2*Power(Mtau,2)*Prho1z*Power(Prho2y,2)*Prho2z - 
             2*Power(Prho1x,2)*Prho1z*Power(Prho2y,2)*Prho2z + 
             2*Power(Prho1y,2)*Prho1z*Power(Prho2y,2)*Prho2z + 
             2*Power(Prho1z,3)*Power(Prho2y,2)*Prho2z + Power(Enu1,4)*Power(Prho2z,2) - 
             2*Power(Enu1,2)*Power(Etau,2)*Power(Prho2z,2) + Power(Etau,4)*Power(Prho2z,2) + 
             2*Power(Enu1,2)*Power(Mtau,2)*Power(Prho2z,2) - 
             2*Power(Etau,2)*Power(Mtau,2)*Power(Prho2z,2) + Power(Mtau,4)*Power(Prho2z,2) + 
             2*Power(Enu1,2)*Power(Prho1x,2)*Power(Prho2z,2) - 
             2*Power(Etau,2)*Power(Prho1x,2)*Power(Prho2z,2) + 
             2*Power(Mtau,2)*Power(Prho1x,2)*Power(Prho2z,2) + Power(Prho1x,4)*Power(Prho2z,2) - 
             2*Power(Enu1,2)*Power(Prho1y,2)*Power(Prho2z,2) - 
             2*Power(Enu2,2)*Power(Prho1y,2)*Power(Prho2z,2) + 
             2*Power(Prho1x,2)*Power(Prho1y,2)*Power(Prho2z,2) + Power(Prho1y,4)*Power(Prho2z,2) - 
             2*Power(Enu1,2)*Power(Prho1z,2)*Power(Prho2z,2) - 
             2*Power(Enu2,2)*Power(Prho1z,2)*Power(Prho2z,2) + 
             4*Power(Etau,2)*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Power(Prho1z,2)*Power(Prho2z,2) - 
             2*Power(Prho1x,2)*Power(Prho1z,2)*Power(Prho2z,2) + 
             2*Power(Prho1y,2)*Power(Prho1z,2)*Power(Prho2z,2) + Power(Prho1z,4)*Power(Prho2z,2) + 
             4*Prho1x*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             2*Power(Prho1y,2)*Power(Prho2x,2)*Power(Prho2z,2) + 
             2*Power(Prho1z,2)*Power(Prho2x,2)*Power(Prho2z,2) - 
             2*Power(Enu1,2)*Prho1y*Prho2y*Power(Prho2z,2) + 
             2*Power(Etau,2)*Prho1y*Prho2y*Power(Prho2z,2) - 
             2*Power(Mtau,2)*Prho1y*Prho2y*Power(Prho2z,2) - 
             2*Power(Prho1x,2)*Prho1y*Prho2y*Power(Prho2z,2) + 
             2*Power(Prho1y,3)*Prho2y*Power(Prho2z,2) + 
             2*Prho1y*Power(Prho1z,2)*Prho2y*Power(Prho2z,2) + 
             2*Power(Prho1y,2)*Power(Prho2y,2)*Power(Prho2z,2) + 
             2*Power(Prho1z,2)*Power(Prho2y,2)*Power(Prho2z,2) - 
             2*Power(Enu1,2)*Prho1z*Power(Prho2z,3) + 2*Power(Etau,2)*Prho1z*Power(Prho2z,3) - 
             2*Power(Mtau,2)*Prho1z*Power(Prho2z,3) - 2*Power(Prho1x,2)*Prho1z*Power(Prho2z,3) + 
             2*Power(Prho1y,2)*Prho1z*Power(Prho2z,3) + 2*Power(Prho1z,3)*Power(Prho2z,3) + 
             Power(Prho1y,2)*Power(Prho2z,4) + Power(Prho1z,2)*Power(Prho2z,4));

    if(tempVal<0){
        //if(sol==1) cout << "Using approximation" << endl;
        tempVal = -tempVal;
    }

    Double_t Px1 = (4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
        4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
        4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
        4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
        8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
        4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
        4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
        4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
        4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
        4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
        4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
        4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
        4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
        4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
        4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
        4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
        4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
        4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
        12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
        4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
        16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
        4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
        4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
        4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
        4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
        4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) - 
        Sqrt(tempVal))/
      (2.*(4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
          8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
          4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
          8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
          4*Power(Prho1y,2)*Power(Prho2z,2)));

    Double_t Py1 = (-(Power(Enu2,2)*Prho1z) + Power(Etau,2)*Prho1z - Power(Mtau,2)*Prho1z + 
        2*Prho1x*Prho1z*Prho2x + Prho1z*Power(Prho2x,2) + 2*Prho1y*Prho1z*Prho2y + 
        Prho1z*Power(Prho2y,2) - Power(Enu1,2)*Prho2z + Power(Etau,2)*Prho2z - Power(Mtau,2)*Prho2z - 
        Power(Prho1x,2)*Prho2z - Power(Prho1y,2)*Prho2z + Power(Prho1z,2)*Prho2z + 
        Prho1z*Power(Prho2z,2) + (Prho1z*Prho2x*
           (4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 
             8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 
             4*Power(Prho1y,2)*Power(Prho2x,3) - 4*Power(Prho1z,2)*Power(Prho2x,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 
             4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 
             12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 4*Power(Prho1y,3)*Prho2x*Prho2y - 
             4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Prho1x,3)*Power(Prho2y,2) + 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 
             4*Prho1x*Prho1y*Power(Prho2y,3) - 4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 
             4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 
             4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 
             4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1z,3)*Prho2x*Prho2z + 
             4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 
             4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Prho1x,3)*Power(Prho2z,2) - 4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 
             4*Prho1x*Prho1z*Power(Prho2z,3) - 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)) - 
        (Prho1x*Prho2z*(4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
             8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
             4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
             4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
             4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
             12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
             16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) - 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)))/(-2*Prho1z*Prho2y + 2*Prho1y*Prho2z);

    Double_t Pz1 = (-(Power(Enu2,2)*Prho1y) + Power(Etau,2)*Prho1y - Power(Mtau,2)*Prho1y + 
        2*Prho1x*Prho1y*Prho2x + Prho1y*Power(Prho2x,2) - Power(Enu1,2)*Prho2y + 
        Power(Etau,2)*Prho2y - Power(Mtau,2)*Prho2y - Power(Prho1x,2)*Prho2y + 
        Power(Prho1y,2)*Prho2y - Power(Prho1z,2)*Prho2y + Prho1y*Power(Prho2y,2) + 
        2*Prho1y*Prho1z*Prho2z + Prho1y*Power(Prho2z,2) + 
        (Prho1y*Prho2x*(4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
             8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
             4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
             4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
             4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
             12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
             16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) - 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)) - 
        (Prho1x*Prho2y*(4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
             8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
             4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
             4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
             4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
             12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
             16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) - 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)))/(2*Prho1z*Prho2y - 2*Prho1y*Prho2z);
        
    Double_t Px2 = (4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 
        4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 
        4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 
        8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 
        4*Power(Prho1y,2)*Power(Prho2x,3) - 4*Power(Prho1z,2)*Power(Prho2x,3) - 
        4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 
        4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 
        4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 
        12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 4*Power(Prho1y,3)*Prho2x*Prho2y - 
        4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 
        4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 
        4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Prho1x,3)*Power(Prho2y,2) + 
        4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 
        4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 
        4*Prho1x*Prho1y*Power(Prho2y,3) - 4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 
        4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 
        4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 
        4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 
        4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1z,3)*Prho2x*Prho2z + 
        4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 
        4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 
        4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 
        4*Power(Prho1x,3)*Power(Prho2z,2) - 4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 
        4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 
        4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 
        4*Prho1x*Prho1z*Power(Prho2z,3) + 
        Sqrt(tempVal))/
      (2.*(4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
          8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
          4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
          8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
          4*Power(Prho1y,2)*Power(Prho2z,2)));

    Double_t Py2 = (-(Power(Enu2,2)*Prho1z) + Power(Etau,2)*Prho1z - Power(Mtau,2)*Prho1z + 
        2*Prho1x*Prho1z*Prho2x + Prho1z*Power(Prho2x,2) + 2*Prho1y*Prho1z*Prho2y + 
        Prho1z*Power(Prho2y,2) - Power(Enu1,2)*Prho2z + Power(Etau,2)*Prho2z - Power(Mtau,2)*Prho2z - 
        Power(Prho1x,2)*Prho2z - Power(Prho1y,2)*Prho2z + Power(Prho1z,2)*Prho2z + 
        Prho1z*Power(Prho2z,2) + (Prho1z*Prho2x*
           (4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 
             8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 
             4*Power(Prho1y,2)*Power(Prho2x,3) - 4*Power(Prho1z,2)*Power(Prho2x,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 
             4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 
             12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 4*Power(Prho1y,3)*Prho2x*Prho2y - 
             4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Prho1x,3)*Power(Prho2y,2) + 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 
             4*Prho1x*Prho1y*Power(Prho2y,3) - 4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 
             4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 
             4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 
             4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1z,3)*Prho2x*Prho2z + 
             4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 
             4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Prho1x,3)*Power(Prho2z,2) - 4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 
             4*Prho1x*Prho1z*Power(Prho2z,3) + 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)) - 
        (Prho1x*Prho2z*(4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
             8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
             4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
             4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
             4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
             12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
             16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) + 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)))/(-2*Prho1z*Prho2y + 2*Prho1y*Prho2z);

    Double_t Pz2 = (-(Power(Enu2,2)*Prho1y) + Power(Etau,2)*Prho1y - Power(Mtau,2)*Prho1y + 
        2*Prho1x*Prho1y*Prho2x + Prho1y*Power(Prho2x,2) - Power(Enu1,2)*Prho2y + 
        Power(Etau,2)*Prho2y - Power(Mtau,2)*Prho2y - Power(Prho1x,2)*Prho2y + 
        Power(Prho1y,2)*Prho2y - Power(Prho1z,2)*Prho2y + Prho1y*Power(Prho2y,2) + 
        2*Prho1y*Prho1z*Prho2z + Prho1y*Power(Prho2z,2) + 
        (Prho1y*Prho2x*(4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
             8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
             4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
             4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
             4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
             12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
             16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) + 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)) - 
        (Prho1x*Prho2y*(4*Power(Enu2,2)*Power(Prho1y,2)*Prho2x - 
             4*Power(Etau,2)*Power(Prho1y,2)*Prho2x + 4*Power(Mtau,2)*Power(Prho1y,2)*Prho2x + 
             4*Power(Enu2,2)*Power(Prho1z,2)*Prho2x - 4*Power(Etau,2)*Power(Prho1z,2)*Prho2x + 
             4*Power(Mtau,2)*Power(Prho1z,2)*Prho2x - 8*Prho1x*Power(Prho1y,2)*Power(Prho2x,2) - 
             8*Prho1x*Power(Prho1z,2)*Power(Prho2x,2) - 4*Power(Prho1y,2)*Power(Prho2x,3) - 
             4*Power(Prho1z,2)*Power(Prho2x,3) - 4*Power(Enu2,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Etau,2)*Prho1x*Prho1y*Prho2y - 4*Power(Mtau,2)*Prho1x*Prho1y*Prho2y + 
             4*Power(Enu1,2)*Prho1y*Prho2x*Prho2y - 4*Power(Etau,2)*Prho1y*Prho2x*Prho2y + 
             4*Power(Mtau,2)*Prho1y*Prho2x*Prho2y + 12*Power(Prho1x,2)*Prho1y*Prho2x*Prho2y - 
             4*Power(Prho1y,3)*Prho2x*Prho2y - 4*Prho1y*Power(Prho1z,2)*Prho2x*Prho2y + 
             4*Prho1x*Prho1y*Power(Prho2x,2)*Prho2y - 4*Power(Enu1,2)*Prho1x*Power(Prho2y,2) + 
             4*Power(Etau,2)*Prho1x*Power(Prho2y,2) - 4*Power(Mtau,2)*Prho1x*Power(Prho2y,2) - 
             4*Power(Prho1x,3)*Power(Prho2y,2) + 4*Prho1x*Power(Prho1y,2)*Power(Prho2y,2) - 
             4*Prho1x*Power(Prho1z,2)*Power(Prho2y,2) - 4*Power(Prho1y,2)*Prho2x*Power(Prho2y,2) - 
             4*Power(Prho1z,2)*Prho2x*Power(Prho2y,2) + 4*Prho1x*Prho1y*Power(Prho2y,3) - 
             4*Power(Enu2,2)*Prho1x*Prho1z*Prho2z + 4*Power(Etau,2)*Prho1x*Prho1z*Prho2z - 
             4*Power(Mtau,2)*Prho1x*Prho1z*Prho2z + 4*Power(Enu1,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Etau,2)*Prho1z*Prho2x*Prho2z + 4*Power(Mtau,2)*Prho1z*Prho2x*Prho2z + 
             12*Power(Prho1x,2)*Prho1z*Prho2x*Prho2z - 4*Power(Prho1y,2)*Prho1z*Prho2x*Prho2z - 
             4*Power(Prho1z,3)*Prho2x*Prho2z + 4*Prho1x*Prho1z*Power(Prho2x,2)*Prho2z + 
             16*Prho1x*Prho1y*Prho1z*Prho2y*Prho2z + 4*Prho1x*Prho1z*Power(Prho2y,2)*Prho2z - 
             4*Power(Enu1,2)*Prho1x*Power(Prho2z,2) + 4*Power(Etau,2)*Prho1x*Power(Prho2z,2) - 
             4*Power(Mtau,2)*Prho1x*Power(Prho2z,2) - 4*Power(Prho1x,3)*Power(Prho2z,2) - 
             4*Prho1x*Power(Prho1y,2)*Power(Prho2z,2) + 4*Prho1x*Power(Prho1z,2)*Power(Prho2z,2) - 
             4*Power(Prho1y,2)*Prho2x*Power(Prho2z,2) - 4*Power(Prho1z,2)*Prho2x*Power(Prho2z,2) + 
             4*Prho1x*Prho1y*Prho2y*Power(Prho2z,2) + 4*Prho1x*Prho1z*Power(Prho2z,3) + 
             Sqrt(tempVal)))/
         (4*Power(Prho1y,2)*Power(Prho2x,2) + 4*Power(Prho1z,2)*Power(Prho2x,2) - 
           8*Prho1x*Prho1y*Prho2x*Prho2y + 4*Power(Prho1x,2)*Power(Prho2y,2) + 
           4*Power(Prho1z,2)*Power(Prho2y,2) - 8*Prho1x*Prho1z*Prho2x*Prho2z - 
           8*Prho1y*Prho1z*Prho2y*Prho2z + 4*Power(Prho1x,2)*Power(Prho2z,2) + 
           4*Power(Prho1y,2)*Power(Prho2z,2)))/(2*Prho1z*Prho2y - 2*Prho1y*Prho2z);

    TLorentzVector vNeutrino1;

    if(sol==1) vNeutrino1.SetPxPyPzE(Px1,Py1,Pz1,Enu1);

    else if(sol==2) vNeutrino1.SetPxPyPzE(Px2,Py2,Pz2,Enu1);

    else {cout << "Invalid solution number." << endl; vNeutrino1.SetPxPyPzE(0,0,0,0);}

    vNeutrino1.Boost(v3Higgs);

    //cout << tempVal << endl;

    return vNeutrino1;
}

/*
 * FUNCTION FOR GETTING THE MOMENTUM OF THE SECOND NEUTRINO
 */

TLorentzVector getNeut2(TLorentzVector vRho1, TLorentzVector vRho2, TLorentzVector vNeutrino1, TLorentzVector vZ, TLorentzVector vInit){
    return vInit - vRho1 - vRho2 - vNeutrino1 - vZ;
}
