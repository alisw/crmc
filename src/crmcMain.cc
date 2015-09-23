#include <CRMCoptions.h>
#include <CRMC.h>
#ifdef WITH_ROOT
#include <OutputPolicyROOT.h>
#endif
#include <OutputPolicyHepMC.h>
#include <OutputPolicyLHE.h>
#include <OutputPolicyNone.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;



int
main(int argc, char **argv)
{

  const CRMCoptions cfg(argc, argv);
  if (cfg.OptionsError()){
    cout << "\nConfiguration Error\n" << endl;
    return 2;
  }

  /*
  ofstream out("out.dat");

  CRMCoptions cfgCpy = cfg;

  for (int i=0; i<50; ++i) {
    const double labE = std::pow(10., 1 + 0.5 * i); // [GeV]

    cfgCpy.SetProjectileMomentum(labE);
    cfgCpy.SetTargetMomentum(0.);

    CRMC2NONE crmc(cfgCpy);
    crmc.init();

    double xsigtot=0,xsigine=0,xsigela=0,xsigdd=0,xsigsd=0,
      xsloela=0,xsigtotaa=0,xsigineaa=0,xsigelaaa=0;
    crmc.GetInterface().crmc_xsection(xsigtot,xsigine,xsigela,xsigdd,xsigsd,
				      xsloela,xsigtotaa,xsigineaa,xsigelaaa);

    out << labE << " " << xsigtot<< " " <<xsigine<< " " <<xsigela<< " " <<xsigdd<< " " <<xsigsd<< " " <<
      xsloela<< " " <<xsigtotaa<< " " <<xsigineaa<< " " <<xsigelaaa << endl;

  }

  out.close();
  return 0;
  */


  switch(cfg.GetOutputMode()) {

#ifdef WITH_ROOT
  case CRMCoptions::eROOT:
    {
      CRMC2ROOT crmc(cfg);
      crmc.init();
      crmc.run();
      crmc.finish();
    }
    break;
#endif

  case CRMCoptions::eHepMC:
  case CRMCoptions::eHepMCGZ:
    {
      CRMC2HepMC crmc(cfg);
      crmc.init();
      crmc.run();
      crmc.finish();
    }
    break;

  case CRMCoptions::eLHE:
  case CRMCoptions::eLHEGZ:
    {
      CRMC2LHE crmc(cfg);
      crmc.init();
      crmc.run();
      crmc.finish();
    }
    break;

  case CRMCoptions::eNone:
    {
      CRMC2NONE crmc(cfg);
      crmc.init();
      crmc.run();
      crmc.finish();
    }
    break;

  }

  // default options
}
