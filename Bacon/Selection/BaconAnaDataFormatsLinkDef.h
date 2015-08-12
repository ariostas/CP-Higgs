#ifndef BACONANA_DATAFORMATS_LINKDEF_H
#define BACONANA_DATAFORMATS_LINKDEF_H
#include "TGenEventInfo.hh"
#include "TGenParticle.hh"
#include "TGenJet.hh"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace baconhep;

#pragma link C++ class baconhep::TGenEventInfo+;
#pragma link C++ class baconhep::TGenParticle+;
#pragma link C++ class baconhep::TGenJet+;
#endif
