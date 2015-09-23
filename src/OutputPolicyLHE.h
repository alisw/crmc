#ifndef _OutputPolicyLHE_h_
#define _OutputPolicyLHE_h_


class CRMCoptions;

class OutputPolicyLHE {

 public:
  OutputPolicyLHE();
  
  void InitOutput(const CRMCoptions& cfg);
  void FillEvent(const CRMCoptions& cfg,const int nEvent);
  void CloseOutput(const CRMCoptions& cfg);

};


#endif
