#ifndef _OutputPolicyHepMC_h_
#define _OutputPolicyHepMC_h_

class CRMCoptions;

namespace HepMC {
  class IO_HEPEVT;
  class IO_GenEvent;
  //class GenEvent;
}

#include <boost/iostreams/filtering_stream.hpp>


class OutputPolicyHepMC {

 public:
  OutputPolicyHepMC();

  void InitOutput(const CRMCoptions& cfg);
  void FillEvent(const CRMCoptions& cfg, const int nEvent);
  void CloseOutput(const CRMCoptions& cfg);

 private:

  void PrintTestEvent(const CRMCoptions& cfg);
  boost::iostreams::filtering_ostream *fOut;
  HepMC::IO_HEPEVT* hepevtio;
  HepMC::IO_GenEvent* ascii_out;
  //  HepMC::GenEvent* fEvtHepMC;

 protected:
// --------------  test observables
  double TotalEnergy ;
  double Multiplicity ;
  double MeanPseudorapidity ;
  double PlateauHeight ;
  double MeanPt ;
  double ErrorTotalEnergy ;
  double ErrorMultiplicity ;
  double ErrorMeanPseudorapidity ;
  double ErrorPlateauHeight ;
  double ErrorMeanPt ;
};


#endif

