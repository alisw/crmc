#ifndef _OutputPolicyNone_h_
#define _OutputPolicyNone_h_

class CRMCoptions;

class OutputPolicyNone {

 public:
  OutputPolicyNone();

  void InitOutput(const CRMCoptions& cfg);
  void FillEvent(const CRMCoptions& cfg, const int nEvent);
  void CloseOutput(const CRMCoptions& cfg);

 private:

  void PrintTestEvent(const CRMCoptions& cfg) {};
  void PrintCrossSections(const CRMCoptions& cfg);

};


#endif
