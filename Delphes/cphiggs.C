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
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMath.h"
#include "Math/Vector3D.h"

#include "neutrinos.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
void histogram(vector<TH1D*>, vector<TString>, const TString, const TString, const TString);
void saveResults();
void analyze(TString, Double_t, Int_t);
Double_t getTheta(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
TLorentzVector getNeut1(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, Int_t);
TLorentzVector getNeut2(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);
Double_t getLikelihood(TH1D*, TH1D*);
Double_t computeUncertainty(TH1D*);

// Initialize storage variables
vector<Double_t> selLep, selHad, massLep, massHad, kinLep, kinHad, sel, mass, kin;
vector<Double_t> selLepErr, selHadErr, massLepErr, massHadErr, kinLepErr, kinHadErr, selErr, massErr, kinErr;
vector<Int_t> signalFlags;
vector<Double_t> crossSections;
vector<TString> sampleNames, sampleFiles;
vector<TH1D*> hTheta, hThetaLep, hThetaHad, genHist;
TF1 *fitFunc = new TF1("theta shape", "[0] + [1]*TMath::Cos(x-2.0*[2])");
TF1 *fitLike = new TF1("likelihood fit", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x");

/*
 * MAIN FUNCTION
 */

void cphiggs(TString inputFile = "xsec.txt"){
    
    cout << "\n\nStarting process...\n\n";
    
    ifstream ifs(inputFile); if(!ifs.is_open()){cout << "Error. File " << inputFile << " not found. Exiting...\n"; return;}
    
    Int_t nSamples = 0;
    TString sampleName, sampleFile, crossSection, signalFlag;
    
    while(ifs >> sampleName >> sampleFile >> crossSection >> signalFlag){
        
        if(sampleName.Contains("#")) continue;
        
        sampleNames.push_back(sampleName);
        sampleFiles.push_back(sampleFile);
        crossSections.push_back(atof(string(crossSection).c_str()));
        signalFlags.push_back(atof(string(signalFlag).c_str()));
        
    }

    for(UInt_t x = 0; x < sampleNames.size(); x++){

        selLep.push_back(0);        selLepErr.push_back(0);
        selHad.push_back(0);        selHadErr.push_back(0);
        sel.push_back(0);           selErr.push_back(0);
        massLep.push_back(0);       massLepErr.push_back(0);
        massHad.push_back(0);       massHadErr.push_back(0);
        mass.push_back(0);          massErr.push_back(0);
        kinLep.push_back(0);        kinLepErr.push_back(0);
        kinHad.push_back(0);        kinHadErr.push_back(0);
        kin.push_back(0);           kinErr.push_back(0);

        hTheta.push_back(new TH1D(sampleNames[x]+"Lep+Had", sampleNames[x]+"Lep+Had", 20, -3.1416, 3.1416));
        hThetaLep.push_back(new TH1D(sampleNames[x]+"Lep", sampleNames[x]+"Lep", 20, -3.1416, 3.1416));
        hThetaHad.push_back(new TH1D(sampleNames[x]+"Had", sampleNames[x]+"Had", 20, -3.1416, 3.1416));

        genHist.push_back(new TH1D(sampleNames[x]+"genHist", sampleNames[x]+"genHist", 100, 0, 250));

        analyze(sampleFiles[x], crossSections[x], x);

    }
    
    // Save results
    saveResults();
    
}

void analyze(TString inputfile, Double_t xsec, Int_t samp){
    
    TString inputFile = inputfile;

    //Event info
    Double_t eventWeight;
    Int_t NLeptons, NJets, NTauJets, hasH, sameCharge;

    //Charged pions
    Double_t CPion1_Pt, CPion2_Pt;
    Double_t CPion1_Eta, CPion2_Eta;
    Double_t CPion1_Phi, CPion2_Phi;
    Double_t CPion1_Mass, CPion2_Mass;

    //Neutral pions
    Double_t NPion11_Pt, NPion12_Pt;
    Double_t NPion11_Eta, NPion12_Eta;
    Double_t NPion11_Phi, NPion12_Phi;
    Double_t NPion11_Mass, NPion12_Mass;
    Double_t NPion21_Pt, NPion22_Pt;
    Double_t NPion21_Eta, NPion22_Eta;
    Double_t NPion21_Phi, NPion22_Phi;
    Double_t NPion21_Mass, NPion22_Mass;

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
    intree->SetBranchAddress("CPion1_Pt",       &CPion1_Pt);
    intree->SetBranchAddress("CPion1_Eta",      &CPion1_Eta);
    intree->SetBranchAddress("CPion1_Phi",      &CPion1_Phi);
    intree->SetBranchAddress("CPion1_Mass",     &CPion1_Mass);
    intree->SetBranchAddress("CPion2_Pt",       &CPion2_Pt);
    intree->SetBranchAddress("CPion2_Eta",      &CPion2_Eta);
    intree->SetBranchAddress("CPion2_Phi",      &CPion2_Phi);
    intree->SetBranchAddress("CPion2_Mass",     &CPion2_Mass);
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

    Double_t tempSelLep = 0, tempSelHad = 0, tempMassLep = 0, tempMassHad = 0, tempKinLep = 0, tempKinHad = 0;
    Double_t tempSelLepErr = 0, tempSelHadErr = 0, tempMassLepErr = 0, tempMassHadErr = 0, tempKinLepErr = 0, tempKinHadErr = 0;

    for (Int_t iEntry=0; iEntry<numberOfEntries; iEntry++) {
        intree->GetEntry(iEntry);

        // Ignore background events containing a Higgs
        if(hasH == 1 && signalFlags.at(samp) > 0) continue;

        TLorentzVector ZP4, HP4, CPion1P4, CPion2P4, NPion11P4, NPion12P4, NPion21P4, NPion22P4, NPion1P4, NPion2P4, InitP4;

        InitP4.SetPxPyPzE(0, 0, 0, 240);

        CPion1P4.SetPtEtaPhiM(CPion1_Pt, CPion1_Eta, CPion1_Phi, CPion1_Mass);
        CPion2P4.SetPtEtaPhiM(CPion2_Pt, CPion2_Eta, CPion2_Phi, CPion2_Mass);

        NPion11P4.SetPtEtaPhiM(NPion11_Pt, NPion11_Eta, NPion11_Phi, NPion11_Mass);
        NPion12P4.SetPtEtaPhiM(NPion12_Pt, NPion12_Eta, NPion12_Phi, NPion12_Mass);
        NPion21P4.SetPtEtaPhiM(NPion21_Pt, NPion21_Eta, NPion21_Phi, NPion21_Mass);
        NPion22P4.SetPtEtaPhiM(NPion22_Pt, NPion22_Eta, NPion22_Phi, NPion22_Mass);

        NPion1P4 = NPion11P4 + NPion12P4;
        NPion2P4 = NPion21P4 + NPion22P4;

        Int_t ZFromLep;
        if(ZLepton1_Pt != 0. && ZLepton2_Pt != 0.){
            TLorentzVector Temp1P4, Temp2P4;
            Temp1P4.SetPtEtaPhiM(ZLepton1_Pt, ZLepton1_Eta, ZLepton1_Phi, ZLepton1_Mass);
            Temp2P4.SetPtEtaPhiM(ZLepton2_Pt, ZLepton2_Eta, ZLepton2_Phi, ZLepton2_Mass);
            ZP4 = Temp1P4 + Temp2P4;
            ZFromLep = 1;
        }
        else if(ZReco_Pt > 0.1 && ZReco_Mass > 0.1){
            ZP4.SetPtEtaPhiM(ZReco_Pt, ZReco_Eta, ZReco_Phi, ZReco_Mass);
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

        genHist[samp]->Fill(ZP4.M());

        HP4 = InitP4 - ZP4;

        // Rho and neutrino 4-vectors
        TLorentzVector Rho1P4, Rho2P4, Neutrino1Sol1P4, Neutrino1Sol2P4, Neutrino2Sol1P4, Neutrino2Sol2P4, Tau1Sol1P4, Tau1Sol2P4, Tau2Sol1P4, Tau2Sol2P4;

        Rho1P4 = CPion1P4 + NPion1P4;
        Rho2P4 = CPion2P4 + NPion2P4;
        Neutrino1Sol1P4 = getNeut1(ZP4, InitP4, Rho1P4, Rho2P4, 1);
        Neutrino2Sol1P4 = getNeut2(Rho1P4, Rho2P4, Neutrino1Sol1P4, ZP4, InitP4);
        Neutrino1Sol2P4 = getNeut1(ZP4, InitP4, Rho1P4, Rho2P4, 2);
        Neutrino2Sol2P4 = getNeut2(Rho1P4, Rho2P4, Neutrino1Sol2P4, ZP4, InitP4);
        Tau1Sol1P4 = Rho1P4 + Neutrino1Sol1P4;
        Tau1Sol2P4 = Rho1P4 + Neutrino1Sol2P4;
        Tau2Sol1P4 = Rho2P4 + Neutrino2Sol1P4;
        Tau2Sol2P4 = Rho2P4 + Neutrino2Sol2P4;

        if(ZFromLep == 1){

            if(sameCharge != 0 || NLeptons != 2) continue;

            tempSelLep += eventWeight;
            tempSelLepErr++;

            if(fabs(ZP4.M() - 91.2) > 5.) continue;
            if(fabs(HP4.M() - 125) > 5.) continue;

            tempMassLep += eventWeight;
            tempMassLepErr++;
            
            if(fabs(ZP4.P()-51.6) > 5.) continue;
            if(Neutrino1Sol1P4.M() < -5 || Neutrino2Sol1P4.M() < -5 || Neutrino1Sol2P4.M() < -5 || Neutrino2Sol2P4.M() < -5) continue;

            tempKinLep += eventWeight;
            tempKinLepErr++;

        }

        else if(ZFromLep == -1){

            if(NLeptons != 0) continue;

            tempSelHad += eventWeight;
            tempSelHadErr ++;

            if(fabs(ZP4.M() - 91.2) > 10.) continue;
            if(fabs(HP4.M() - 125) > 10.) continue;

            tempMassHad += eventWeight;
            tempMassHadErr++;

            if(fabs(ZP4.P()-51.6) > 5.) continue;
            if(Neutrino1Sol1P4.M() < -5 || Neutrino2Sol1P4.M() < -5 || Neutrino1Sol2P4.M() < -5 || Neutrino2Sol2P4.M() < -5) continue;

            tempKinHad += eventWeight;
            tempKinHadErr++;

        }

        // Compute theta variable
        Double_t thetaSol1, thetaSol2;
        thetaSol1 = getTheta(CPion1P4, NPion1P4, CPion2P4, NPion2P4, Tau1Sol1P4, Tau2Sol1P4);
        thetaSol2 = getTheta(CPion1P4, NPion1P4, CPion2P4, NPion2P4, Tau1Sol2P4, Tau2Sol2P4);

        hTheta.at(samp)->Fill(thetaSol1, eventWeight/2.);
        hTheta.at(samp)->Fill(thetaSol2, eventWeight/2.);

        (ZFromLep == 1 ? hThetaLep : hThetaHad).at(samp)->Fill(thetaSol1, eventWeight/2.);
        (ZFromLep == 1 ? hThetaLep : hThetaHad).at(samp)->Fill(thetaSol2, eventWeight/2.);

    } // end event loop

    selLep.at(samp) += tempSelLep;
    selHad.at(samp) += tempSelHad;
    sel.at(samp) += tempSelLep + tempSelHad;
    if(tempSelLep > 0) selLepErr.at(samp) += tempSelLep/sqrtf(tempSelLepErr);
    if(tempSelHad > 0) selHadErr.at(samp) += tempSelHad/sqrtf(tempSelHadErr);
    if(tempSelLep + tempSelHad > 0) selErr.at(samp) += (tempSelLep + tempSelHad)/sqrtf(tempSelLepErr + tempSelHadErr);
    massLep.at(samp) += tempMassLep;
    massHad.at(samp) += tempMassHad;
    mass.at(samp) += tempMassLep + tempMassHad;
    if(tempMassLep > 0) massLepErr.at(samp) += tempMassLep/sqrtf(tempMassLepErr);
    if(tempMassHad > 0) massHadErr.at(samp) += tempMassHad/sqrtf(tempMassHadErr);
    if(tempMassLep + tempMassHad > 0) massErr.at(samp) += (tempMassLep + tempMassHad)/sqrtf(tempMassLepErr + tempMassHadErr);
    kinLep.at(samp) += tempKinLep;
    kinHad.at(samp) += tempKinHad;
    kin.at(samp) += tempKinLep + tempKinHad;
    if(tempKinLep > 0) kinLepErr.at(samp) += tempKinLep/sqrtf(tempKinLepErr);
    if(tempKinHad > 0) kinHadErr.at(samp) += tempKinHad/sqrtf(tempKinHadErr);
    if(tempKinLep + tempKinHad > 0) kinErr.at(samp) += (tempKinLep + tempKinHad)/sqrtf(tempKinLepErr + tempKinHadErr);

    TString out = TString::Format("%4.2f", tempKinLep + tempKinHad);
    out.Resize(9);

    cout << "\e[A";
    cout << inputfile << (signalFlags.at(samp)<=0 ? "\033[1;32m" : (out == "0       " ?  "\033[1;34m": "\033[1;31m")) << out << "\033[0m" << " events passed" << endl;

}

void saveResults()
{
    cout << endl << endl;   

    gStyle->SetOptStat(kFALSE);

    cout << "\n\n\nProcess finished\n";

    cout << "\n\n\033[1;34mEvent yields\033[0m\n\n";
    
    for(UInt_t x = 0; x<sampleNames.size();x++){
        
        cout << (signalFlags.at(x)<=0 ? "\033[1;32m" : "\033[1;31m") << sampleNames.at(x) << "\033[0m\n";
        cout << "Leptonic events:" << endl;
        cout << "  Events after selection: " << selLep.at(x) << " +- " << selLepErr.at(x) << endl;
        cout << "  Events after mass cuts: " << massLep.at(x) << " +- " << massLepErr.at(x) << endl;
        cout << "  Events after kine cuts: " << kinLep.at(x) << " +- " << kinLepErr.at(x) << endl << endl;
        cout << "Hadronic events:" << endl;
        cout << "  Events after selection: " << selHad.at(x) << " +- " << selHadErr.at(x) << endl;
        cout << "  Events after mass cuts: " << massHad.at(x) << " +- " << massHadErr.at(x) << endl;
        cout << "  Events after kine cuts: " << kinHad.at(x) << " +- " << kinHadErr.at(x) << endl << endl;
        cout << "All events:" << endl;
        cout << "  Events after selection: " << sel.at(x) << " +- " << selErr.at(x) << endl;
        cout << "  Events after mass cuts: " << mass.at(x) << " +- " << massErr.at(x) << endl;
        cout << "  Events after kine cuts: " << kin.at(x) << " +- " << kinErr.at(x) << endl << endl;
        
    }

    Double_t uncLep = computeUncertainty(hThetaLep.at(0));
    Double_t uncHad = computeUncertainty(hThetaHad.at(0));
    Double_t unc = computeUncertainty(hTheta.at(0));

    cout << "\n\n\033[1;34mCP measurement uncertainties\033[0m\n\n";

    cout << "Uncertainty from leptonic events: " << uncLep << " degrees" << endl;
    cout << "Uncertainty from hadronic events: " << uncHad << " degrees" << endl;
    cout << "Uncertainty from all events: " << unc << " degrees" << endl;

    cout << "\n\n\033[1;34mSaving histograms\033[0m\n\n";

    histogram(hTheta, sampleNames, "#Theta variable", "Normalized yield", "Theta");
    histogram(hThetaLep, sampleNames, "#Theta variable (leptonic events)", "Normalized yield", "ThetaLep");
    histogram(hThetaHad, sampleNames, "#Theta variable (hadronic events)", "Normalized yield", "ThetaHad");
    histogram(genHist, sampleNames, "Z mass", "Normalized yield", "Zmass");

    cout << "Done\nExiting...\n\n\n";

}

Double_t computeUncertainty(TH1D *thetaHisto){

    fitFunc->SetParameter(0, thetaHisto->GetMaximum() - (thetaHisto->GetMaximum() - thetaHisto->GetMinimum())/2.);
    fitFunc->SetParameter(1, -(thetaHisto->GetMaximum() - thetaHisto->GetMinimum())/2.);
    fitFunc->SetParameter(2, 0);

    thetaHisto->Fit(fitFunc, "ME");

    fitFunc->SetParameter(2, 0);

    TH1D *thetaScalar = new TH1D("theta scalar", "theta scalar", 20, -3.1416, 3.1416);

    for(Int_t x = 0; x < 20; x++){
        thetaScalar->Fill(3.1416*2.0/20.0*(0.5 + x - 10), fitFunc->Eval(3.1416*2.0/20.0*(0.5 + x - 10)));
    }

    Int_t nsteps = 21;
    Int_t middle = (nsteps-1)/2;
    Double_t phase[nsteps], loglike[nsteps];
    for(Int_t x = 0; x < nsteps; x++){
        fitFunc->SetParameter(2, (x - middle)*3.1416/2.0/Double_t(middle));
        TH1D *thetaTemp = new TH1D("theta temp", "theta temp", 20, -3.1416, 3.1416);
        for(Int_t y = 0; y < 20; y++){
            thetaTemp->Fill(3.1416*2.0/20.0*(0.5 + y - 10), fitFunc->Eval(3.1416*2.0/20.0*(0.5 + y - 10)));
        }
        phase[x] = (x - middle)*3.1416/2.0/Double_t(middle);
        loglike[x] = -Log(getLikelihood(thetaScalar, thetaTemp));
        delete thetaTemp;
    }
    delete thetaScalar;

    TGraph *loglikeGraph = new TGraph(nsteps, phase, loglike);
    loglikeGraph->Fit(fitLike, "ME", "", 0, 3.1416/2.0);
    Double_t uncert = fitLike->GetX(0.5, 0, 3.1416);

    delete loglikeGraph;

    return uncert*180./3.1415;
}

Double_t getLikelihood(TH1D *histo1, TH1D *histo2){

    Double_t Likelihood=1.;

    for(Int_t x = 1; x <= 20; x++){

        Likelihood*=PoissonI(histo1->GetBinContent(x), histo2->GetBinContent(x));
        Likelihood/=PoissonI(histo1->GetBinContent(x), histo1->GetBinContent(x));
    }

    return Likelihood;
}

/*
 * FUNCTION FOR SAVING MULTIPLE HISTOGRAMS
 */

void histogram(vector<TH1D*> histos, vector<TString> histNames, const TString xTitle, const TString yTitle, const TString name){

    TCanvas *can = new TCanvas(name, name, 1600, 900);   
    
    if(histos.size()!=histNames.size()){ cout << "Number of histograms and names don't match." << endl; return;}

    Double_t max=0, min=1;

    for(UInt_t i=0; i<histos.size(); i++){
        Double_t integral = histos.at(i)->Integral();
        if(integral != 0.){
            // histos.at(i)->Scale(Double_t(1)/integral);
            if(histos.at(i)->GetMaximum()>max) max=histos.at(i)->GetMaximum();
            if(histos.at(i)->GetMinimum()<min) min=histos.at(i)->GetMinimum();
        }
    }

    max*=1.25;
    min*=0.75;
    //if(name == "Theta"){max=.1; min=0;}

    vector<Int_t> colors, lines;
    colors.push_back(kRed); colors.push_back(kBlue); colors.push_back(kBlack); colors.push_back(kGreen+2); colors.push_back(kRed); colors.push_back(kBlack); 
    lines.push_back(1); lines.push_back(1); lines.push_back(2); lines.push_back(7); lines.push_back(9);

    for(UInt_t i=0; i<histos.size(); i++){
        if(histNames.at(i) == "Backgrounds") continue;
        histos.at(i)->SetMaximum(max);
        histos.at(i)->SetMinimum(min);
        histos.at(i)->SetLineWidth(4);
        histos.at(i)->SetLineColor(colors.at(i%6));
        histos.at(i)->SetLineStyle(lines.at(i%5));
        if(i==0){
            histos.at(i)->Draw("][ histo");
        }
        else histos.at(i)->Draw("same ][ histo");
        
    }

    histos.at(0)->GetXaxis()->SetTitle(xTitle);
    histos.at(0)->GetXaxis()->CenterTitle();
    histos.at(0)->GetXaxis()->SetTitleSize(0.055);
    histos.at(0)->GetXaxis()->SetTitleOffset(0.87);
    histos.at(0)->GetXaxis()->SetLabelOffset(0.010);
    histos.at(0)->GetXaxis()->SetLabelSize(0.05);
    histos.at(0)->GetYaxis()->SetTitle(yTitle);
    histos.at(0)->GetYaxis()->CenterTitle();
    histos.at(0)->GetYaxis()->SetTitleSize(0.055);
    histos.at(0)->GetYaxis()->SetTitleOffset(0.95);
    histos.at(0)->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    histos.at(0)->SetTitle("");
    can->Update();

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.70,0.675,0.90,0.875);
    //leg->SetTextFont(72);
    leg->SetTextSize(0.04);
    leg->SetFillColor(kWhite);
    for(UInt_t i=0; i<histos.size(); i++){
        if(histNames.at(i) == "Backgrounds") continue;
        leg->AddEntry(histos.at(i),histNames.at(i),"l");
    }
    leg->Draw("same");

    TString filename=name; filename+=".png";

    can->SaveAs(filename);
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
