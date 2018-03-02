#ifndef _CRMCoptions_h_
#define _CRMCoptions_h_

#include <string>

class CRMCoptions {

  template<typename> friend class CRMC;
  friend class OutputPolicyROOT;
  friend class OutputPolicyHepMC;
  friend class OutputPolicyLHE;
  friend class OutputPolicyNone;

  CRMCoptions();

 public:

  enum EOutputMode {
    eHepMC,
    eHepMCGZ,
    eLHE,
    eLHEGZ,
    eROOT,
    eNone,
  };

  CRMCoptions(int argc, char** argv);

  bool OptionsError() const { return fError; }
  void DumpConfig() const;
  EOutputMode GetOutputMode() const { return fOutputMode; }
  std::string GetOutputTypeEnding() const;
  std::string GetOutputFileName() const;

  std::string ParticleName(const int pid) const;

  //std::string GetFilter() const { return fFilter; }

  void SetProjectileMomentum(const double p) { fProjectileMomentum = p; }
  void SetTargetMomentum(const double p) { fTargetMomentum = p; }

 protected:

  bool fError;
  EOutputMode fOutputMode;

  // real data members

  int fNCollision;
  int fSeed;
  int fProjectileId;
  int fTargetId;
  int fHEModel;
  int fTypout;
  double fProjectileMomentum;
  double fTargetMomentum;
  double fSqrts;
  std::string fParamFileName;
  std::string fOutputFileName;

  bool fProduceTables;
  bool fSeedProvided;
  //std::string fFilter;
  bool fTest;
  bool fCSMode;

 private:

  void CheckEnvironment();
  void ParseOptions(int argc, char** argv);

};

#endif
