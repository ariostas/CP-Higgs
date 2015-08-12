{
	if (!TClass::GetDict("TGenEventInfo")) {
		gROOT->ProcessLine(".L TGenEventInfo.cc+");
	}
	
	if (!TClass::GetDict("TGenJet")) {
		gROOT->ProcessLine(".L TGenJet.cc+");
	} 

	if (!TClass::GetDict("TGenParticle")) {
		gROOT->ProcessLine(".L TGenParticle.cc+");
	} 

}

