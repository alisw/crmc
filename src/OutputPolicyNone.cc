#include <OutputPolicyNone.h>

#include <CRMCoptions.h>
#include <CRMCinterface.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>

#include <CRMCconfig.h> //cmake generated


using namespace std;


OutputPolicyNone::OutputPolicyNone()
{
}


void
OutputPolicyNone::InitOutput(const CRMCoptions& cfg)
{
}


void
OutputPolicyNone::FillEvent(const CRMCoptions& cfg, const int nEvent)
{
  //#ifdef HEPMC_HAS_CROSS_SECTION
  // set cross section information for this event
  //None::GenCrossSection theCrossSection;
  //theCrossSection.set_cross_section(double(gCRMC_data.sigineaa)*1e9); //required in pB
  //fEvtNone->set_cross_section(theCrossSection);
  //#endif

}


void
OutputPolicyNone::CloseOutput(const CRMCoptions& cfg)
{
  if(cfg.fCSMode) PrintCrossSections(cfg);
}

void
OutputPolicyNone::PrintCrossSections(const CRMCoptions& cfg)
{
    cout << "\n          >> Cross Sections <<\n\n"

         << "  sqrt(s/GeV**2)=" << cfg.fSqrts << "\n\n"

         << "  Total Cross Section (mb):    " << gCRMC_data.sigtot << "\n"
         << "  Elastic Cross Section (mb):  " << gCRMC_data.sigela << "\n"
         << "  Inel. Cross Section (mb) :   " << gCRMC_data.sigine << "\n" ;
    if(cfg.fProjectileId>1 || cfg.fTargetId>1)
    cout << "  Inel. AA Cross Section (mb): " << gCRMC_data.sigineaa << "\n"
         << "  Elastic AA Cross Section (mb): " << gCRMC_data.sigelaaa << "\n"
         << "  Total AA Cross Section (mb): " << gCRMC_data.sigtotaa << "\n" ;
}
