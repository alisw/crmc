#include <OutputPolicyROOT.h>

#include <CRMCoptions.h>
#include <CRMCinterface.h>

#include <TTree.h>
#include <TFile.h>

#include <iostream>

using namespace std;


OutputPolicyROOT::OutputPolicyROOT()
  :  fSigmaPairTot(0),
     fSigmaPairInel(0),
     fSigmaPairEl(0),
     fSigmaTot(0),
     fSigmaInel(0),
     fSigmaEl(0),
     fFile(0),
     fHead(0),
     fParticle(0)
{
}


void 
OutputPolicyROOT::InitOutput(const CRMCoptions& cfg) 
{

  fFile = new TFile(cfg.GetOutputFileName().c_str(), "RECREATE");

  // file header
  
  fHead = new TTree("Header", "run header");

  fHead->Branch("Seed", const_cast<int*>(&cfg.fSeed), "Seed/I");     // random seedj
  fHead->Branch("ProjectileId", const_cast<int*>(&cfg.fProjectileId), "Projectile/I");   // beam/projectile Id
  fHead->Branch("ProjectileMomentum", const_cast<double*>(&cfg.fProjectileMomentum), "ProjectileMomentum/D");             // beam momentum
  fHead->Branch("TargetMomentum", const_cast<double*>(&cfg.fTargetMomentum),	"TargetMomentum/D");             // target momentum
  fHead->Branch("TargetId", const_cast<int*>(&cfg.fTargetId), "TargetId/I");                 // target Id
  fHead->Branch("HEModel", const_cast<int*>(&cfg.fHEModel), "HEModel/I");   // HE model flag
  fHead->Branch("sigmaPairTot", const_cast<double*>(&fSigmaPairTot), "sigmaPairTot/D");   // projectile-Nucleon tot sigma
  fHead->Branch("sigmaPairInel", const_cast<double*>(&fSigmaPairInel), "sigmaPairInel/D");   // projectile-Nucleon inel sigma
  fHead->Branch("sigmaPairEl", const_cast<double*>(&fSigmaPairEl), "sigmaPairEl/D");   // projectile-Nucleon el sigma
  fHead->Branch("sigmaTot", const_cast<double*>(&fSigmaTot), "sigmaTot/D");   // overal tot sigma
  fHead->Branch("sigmaInel", const_cast<double*>(&fSigmaInel), "sigmaInel/D");   // overal inel sigma
  fHead->Branch("sigmaEl", const_cast<double*>(&fSigmaEl), "sigmaEl/D");   // overal el sigma
  
  
  
  // particle list
  fParticle = new TTree("Particle","particles produced");
  
  fParticle->Branch("nPart", &gCRMC_data.fNParticles, "nPart/I");
  fParticle->Branch("ImpactParameter", &gCRMC_data.fImpactParameter, "ImpactParameter/D");
  fParticle->Branch("pdgid", gCRMC_data.fPartId, "pdgid[nPart]/I");
  fParticle->Branch("status", gCRMC_data.fPartStatus, "status[nPart]/I");
  fParticle->Branch("px", gCRMC_data.fPartPx, "px[nPart]/D");
  fParticle->Branch("py", gCRMC_data.fPartPy, "py[nPart]/D");
  fParticle->Branch("pz", gCRMC_data.fPartPz, "pz[nPart]/D");
  fParticle->Branch("E", gCRMC_data.fPartEnergy, "E[nPart]/D");
  fParticle->Branch("m", gCRMC_data.fPartMass, "m[nPart]/D");
}


void
OutputPolicyROOT::FillEvent(const CRMCoptions& cfg,const int nEvent)
{
  if (fSigmaPairTot==0) { 
    fSigmaPairTot = gCRMC_data.sigtot;
    fSigmaPairInel = gCRMC_data.sigine;
    fSigmaPairEl = gCRMC_data.sigela;
    fSigmaTot = gCRMC_data.sigtotaa;
    fSigmaInel = gCRMC_data.sigineaa;    
    fSigmaEl = gCRMC_data.sigelaaa;
    fHead->Fill(); // do only once
  }
  
  fParticle->Fill();
}


void
OutputPolicyROOT::CloseOutput(const CRMCoptions& cfg)
{      
  fFile->cd();
  fFile->Write();
  fFile->Close();
}
