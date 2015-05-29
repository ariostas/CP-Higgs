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

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include "Dataset.h"
#include "CalcT.h"

using namespace std;

// Declare functions
void histogram(TH1D*, const TString, TCanvas*, const TString, const TString, const TString);
void histogram(vector<TH1D*>, vector<TString>, TCanvas*, const TString, const TString, const TString);
void saveResults();
void analyze(TString, Double_t, Int_t);
Double_t deltaR(const Float_t, const Float_t, const Float_t, const Float_t);

// Initialize histograms
TH1D *hThetaS0 = new TH1D("hThetaS1", "hThetaS1", 20, -3.14, 3.14);
TH1D *hThetaSpi4 = new TH1D("hThetaS2", "hThetaS2", 20, -3.14, 3.14);
TH1D *hThetaSpi2 = new TH1D("hThetaS3", "hThetaS3", 20, -3.14, 3.14);
TH1D *hThetaS3pi4 = new TH1D("hThetaS4", "hThetaS4", 20, -3.14, 3.14);
TH1D *hThetaB = new TH1D("hThetaB", "hThetaB", 20, -3.14, 3.14);

TH1D *genhistoS1 = new TH1D("histoS1", "histoS1", 50, 50, 200);
TH1D *genhistoS2 = new TH1D("histoS2", "histoS2", 50, 50, 200);
TH1D *genhistoS3 = new TH1D("histoS3", "histoS3", 50, 50, 200);
TH1D *genhistoS4 = new TH1D("histoS4", "histoS4", 50, 50, 200);
TH1D *genhistoB = new TH1D("histoB", "histoB", 50, 50, 200);

// Initialize data sets
Dataset data0(1,10000), datapi4(1,10000), datapi2(1,10000), data3pi4(1,10000);
vector<Dataset> datasets;

// Initialyze storage variables
vector<Double_t> total, selection, kinematicCuts, massCuts;
vector<Double_t> totalError, selectionError, kinematicCutsError, massCutsError;
vector<Int_t> signalFlags;
vector<TString> sampleNames;
vector<TH1D*> hTheta, genhistos;
vector<TString> histogramNames;

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

    genhistos.push_back(genhistoB);
    genhistos.push_back(genhistoS1);
    genhistos.push_back(genhistoS2);
    genhistos.push_back(genhistoS3);
    genhistos.push_back(genhistoS4);

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

    Float_t eventWeight;

    TLorentzVector vtau1, vtau2, vcpion1, vcpion2, vnpion1, vnpion2, vHiggs, vrho1, vrho2, vq1, vq2;
    TVector3 v3Higgs, vE1, vE2, vB1, vB2, v3tau1, v3tau2;
    UInt_t nProngTau1=0, nProngTau2=0, nLeptons=0, zToLep=0;

    Float_t y1, y2;
    Float_t r=0.14;
    Float_t theta;

    TFile* infile = new TFile(inputFileTemp);
    assert(infile);
    TTree* intree = (TTree*) infile->Get("Events");
    assert(intree);

    intree->SetBranchAddress("nLeptons",		&nLeptons);
    intree->SetBranchAddress("eventWeight",		&eventWeight);
    intree->SetBranchAddress("nProngTau1",		&nProngTau1);
    intree->SetBranchAddress("nProngTau2",		&nProngTau2);
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

    Double_t tempSelection=0, tempSelectionError=0;

    for (Long64_t iEntry=0; iEntry<intree->GetEntries(); iEntry++)   // Event loop
    {
        intree->GetEntry(iEntry);

        if(pions1_pt==0 || pions2_pt==0 || neutpions1_pt==0 || neutpions2_pt==0 || nProngTau1!=1 || nProngTau2!=1) continue;
        if(genTau1_pt < 10 || genTau2_pt< 10) continue;

        // Taus 4-vectors
        vtau1.SetPtEtaPhiM(genTau1_pt, genTau1_eta, genTau1_phi, genTau1_mass);
        vtau2.SetPtEtaPhiM(genTau2_pt, genTau2_eta, genTau2_phi, genTau2_mass);

        // Charged pions 4-vectors
        vcpion1.SetPtEtaPhiM(pions1_pt, pions1_eta, pions1_phi, pions1_mass);
        vcpion2.SetPtEtaPhiM(pions2_pt, pions2_eta, pions2_phi, pions2_mass);

        // Neutral pions 4-vectors
        vnpion1.SetPtEtaPhiM(neutpions1_pt, neutpions1_eta, neutpions1_phi, neutpions1_mass);
        vnpion2.SetPtEtaPhiM(neutpions2_pt, neutpions2_eta, neutpions2_phi, neutpions2_mass);

        // Rho 4-vectors
        vrho1=vcpion1+vnpion1;
        vrho2=vcpion2+vnpion2;

        // Higgs 4- and 3-vector
        vHiggs=vtau1+vtau2;
        v3Higgs.SetXYZ(-vHiggs.Px()/vHiggs.E(), -vHiggs.Py()/vHiggs.E(), -vHiggs.Pz()/vHiggs.E());

        vq1=vcpion1-vnpion1;
        vq2=vcpion2-vnpion2;

        y1=(vq1.Dot(vtau1))/(vrho1.Dot(vtau1));
        y2=(vq2.Dot(vtau2))/(vrho2.Dot(vtau2));

        // Boost to Higgs frame
        vcpion1.Boost(v3Higgs);
        vcpion2.Boost(v3Higgs);
        vnpion1.Boost(v3Higgs);
        vnpion2.Boost(v3Higgs);
        vtau1.Boost(v3Higgs);
        vtau2.Boost(v3Higgs);

        // Set up 3-vectors in Higgs frame
        TVector3 v3cpion1, v3cpion2, v3npion1, v3npion2;

        v3tau1.SetXYZ(vtau1.Px()/vtau1.E(), vtau1.Py()/vtau1.E(), vtau1.Pz()/vtau1.E());
        v3tau2.SetXYZ(vtau2.Px()/vtau2.E(), vtau2.Py()/vtau2.E(), vtau2.Pz()/vtau2.E());

        v3cpion1.SetXYZ(vcpion1.Px()/vcpion1.E(), vcpion1.Py()/vcpion1.E(), vcpion1.Pz()/vcpion1.E());
        v3cpion2.SetXYZ(vcpion2.Px()/vcpion2.E(), vcpion2.Py()/vcpion2.E(), vcpion2.Pz()/vcpion2.E());
        v3npion1.SetXYZ(vnpion1.Px()/vnpion1.E(), vnpion1.Py()/vnpion1.E(), vnpion1.Pz()/vnpion1.E());
        v3npion2.SetXYZ(vnpion2.Px()/vnpion2.E(), vnpion2.Py()/vnpion2.E(), vnpion2.Pz()/vnpion2.E());

        TVector3 tempE1, tempE2;

        tempE1=(y1-r)*v3cpion1-(y1+r)*v3npion1;
        tempE2=(y2-r)*v3cpion2-(y2+r)*v3npion2;

        vE1=tempE1-(tempE1.Dot(v3tau1))*v3tau1;
        vE2=tempE2-(tempE2.Dot(v3tau2))*v3tau2;

        vB1=v3tau1.Cross(vE1);
        vB2=v3tau2.Cross(vE2);

        // Compute theta variable
        theta=TMath::Sign(Double_t(1),v3tau1.Dot(vE2.Cross(vE1)))*TMath::ACos(vE1.Dot(vE2)/(vE1.Mag()*vE2.Mag()));

        // Fill histogram for theta variable
        hTheta.at(signalFlags.at(samp))->Fill(theta);

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
        vInit.SetPtEtaPhiE(0,0,0,240);
        vZ = v1 + v2;

        vHiggs = vInit-vZ;

        if(vHiggs.M()<120) continue;

        genhistos.at(signalFlags.at(samp))->Fill(vHiggs.M());

        vector<double> vars;
        vars.push_back(theta);
        //vars.push_back(vHiggs.M());

        if(signalFlags.at(samp)==0 /*&& iEntry < intree->GetEntries()/4*/){
            datasets.at(0).add(vars);
            datasets.at(1).add(vars);
            datasets.at(2).add(vars);
            datasets.at(3).add(vars);
        }
        else /*if(iEntry < intree->GetEntries()/4)*/ datasets.at(signalFlags.at(samp)-1).add(vars);

        tempSelection+=eventWeight;
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
    hpNames.push_back("0 and 0"); hpNames.push_back("0 and pi/4"); hpNames.push_back("0 and pi/2"); hpNames.push_back("0 and 3pi/4");

    for(Int_t i=0; i<4; i++){
        Double_t tval = calcT(datasets.at(0),datasets.at(i));
        cout << "T of " << hpNames.at(i) << ": " << tval << endl;

        Double_t ptval;
        srand(22);
        Int_t pval=0;
        Int_t nperm=200;
        for(Int_t p=0; p<nperm; p++){
            ptval = permCalcT(datasets.at(0),datasets.at(i));
            hp.at(i)->Fill(ptval);
            if(ptval > tval) pval++;
        }
        cout << "p of " << hpNames.at(i) << ": " << pval/Double_t(nperm) << endl << endl;

    }

    TCanvas *c1 = new TCanvas("Histogram", "Histogram", 1000, 600);

    gStyle->SetOptStat(kFALSE);

    histogram(hTheta, histogramNames, c1, "#Theta variable", "Fraction", "Theta.jpg");
    histogram(genhistos, histogramNames, c1, "H mass", "Fraction", "histo.jpg");

    for(Int_t i=0; i<4; i++){
        TString filename="p_"; filename+=i; filename+=".jpg";
        histogram(hp.at(i), hpNames.at(i), c1, "p distribution", "Fraction", filename);
    }

    /*TFile f("histos.root","new");
    hThetaS0->Write();
    hThetaSpi4->Write();
    hThetaSpi2->Write();
    hThetaS3pi4->Write();
    hThetaB->Write();
    f.Close();*/

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

    can->SaveAs(name);
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
