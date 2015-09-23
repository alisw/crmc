#ifndef _OutputPolicyROOT_h_
#define _OutputPolicyROOT_h_

class TTree;
class TFile;

class CRMCoptions;

class OutputPolicyROOT {

 public:
  OutputPolicyROOT();
  
  void InitOutput(const CRMCoptions& cfg);
  void FillEvent(const CRMCoptions& cfg,const int nEvent);
  void CloseOutput(const CRMCoptions& cfg);

 protected:

  double fSigmaPairTot;
  double fSigmaPairInel;
  double fSigmaPairEl;
  double fSigmaTot;
  double fSigmaInel;
  double fSigmaEl;

  TFile* fFile;
  TTree* fHead;
  TTree* fParticle;
};


#endif
