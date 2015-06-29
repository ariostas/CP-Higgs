#include "TMath.h"

using namespace TMath;

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
