#ifndef __CRMC_H
#define __CRMC_H

#include <CRMCinterface.h>
//#include <CRMCfilter.h>

// //////////
// //////////

class CRMCoptions;

/*
class VCRMC {

  VCRMC();
 public:

  static VCRMC* Create(const CRMCoptions& cfg);

  bool init();
  bool run();
  bool finish();

  void CleanVector();

};
*/

template<class OutputPolicy>
class CRMC : public OutputPolicy {

  CRMC();

 public:
  CRMC(const CRMCoptions& cfg);

  bool init();
  bool run();
  bool finish();

  CRMCinterface& GetInterface() { return fInterface; }

 private:

  const CRMCoptions& fCfg;
  CRMCinterface fInterface;

  //CRMCfilter fFilter;

};



#ifdef WITH_ROOT
#include <OutputPolicyROOT.h>
typedef CRMC<OutputPolicyROOT> CRMC2ROOT;
#endif

#include <OutputPolicyHepMC.h>
typedef CRMC<OutputPolicyHepMC> CRMC2HepMC;

#include <OutputPolicyLHE.h>
typedef CRMC<OutputPolicyLHE> CRMC2LHE;

#include <OutputPolicyNone.h>
typedef CRMC<OutputPolicyNone> CRMC2NONE;


#endif
