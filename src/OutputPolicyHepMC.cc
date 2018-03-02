#include <OutputPolicyHepMC.h>

#include <CRMCoptions.h>
#include <CRMCinterface.h>

// TODO check for has hepmc
#include <HepMC/GenEvent.h>
#include <HepMC/HeavyIon.h>
#include <HepMC/HEPEVT_Wrapper.h>
#include <HepMC/IO_GenEvent.h>
#include <HepMC/IO_HEPEVT.h>
#include <HepMC/PdfInfo.h>

#ifdef HEPMC_HAS_CROSS_SECTION
#include <HepMC/GenCrossSection.h>
#endif

#ifdef HEPMC_HAS_UNITS
#include <HepMC/Units.h>
#endif

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <CRMCconfig.h> //cmake generated

namespace io = boost::iostreams;

using namespace std;


OutputPolicyHepMC::OutputPolicyHepMC()
{
}


void
OutputPolicyHepMC::InitOutput(const CRMCoptions& cfg)
{
  boost::filesystem::path oldFile(cfg.GetOutputFileName());
  if(!boost::filesystem::is_other(cfg.GetOutputFileName())) //protect fifo file
    boost::filesystem::remove(oldFile); //before liboost v1.44 truncate does not seem to work properly in boost

  //io::filtering_ostream out; //top to bottom order
  fOut = new io::filtering_ostream();
  if (cfg.fOutputMode==CRMCoptions::eHepMCGZ)
    fOut->push(io::gzip_compressor(io::zlib::best_compression));
  fOut->push(io::file_descriptor_sink(cfg.GetOutputFileName()), ios_base::trunc);

  // Instantiate an IO strategy for reading from HEPEVT.
  hepevtio = new HepMC::IO_HEPEVT();
  // Instantiate an IO strategy to write the data to file
  ascii_out = new HepMC::IO_GenEvent(*fOut);

  //We need to explicitly pass this information to the
  //  HEPEVT_Wrapper.
  HepMC::HEPEVT_Wrapper::set_max_number_entries(gCRMC_data.fMaxParticles); //as used in crmc-aaa.f!!!
  HepMC::HEPEVT_Wrapper::set_sizeof_real(8); //as used in crmc-aaa.f!!!
  if (cfg.fTest){
    TotalEnergy=0. ;
    Multiplicity=0. ;
    MeanPseudorapidity=0. ;
    PlateauHeight=0. ;
    MeanPt=0. ;
    ErrorTotalEnergy=0. ;
    ErrorMultiplicity=0. ;
    ErrorMeanPseudorapidity=0. ;
    ErrorPlateauHeight=0. ;
    ErrorMeanPt=0. ;
  }
}


void
OutputPolicyHepMC::FillEvent(const CRMCoptions& cfg, const int nEvent)
{
 //void AddEvent const int nEvent, HepMC::IO_HEPEVT& hepevtio, HepMC::IO_GenEvent& ascii_out

#ifdef HEPMC_HAS_UNITS
 HepMC::GenEvent* fEvtHepMC = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
#else
 HepMC::GenEvent* fEvtHepMC = new HepMC::GenEvent();
#endif

 const bool res = hepevtio->fill_next_event(fEvtHepMC); // here hepevt_ a COMMON fortran block gets read in
 if (! res) {
   // delete fEvtHepMC;
   throw std::runtime_error("!!!Could not read next event");
 }

 //   if(!HepMC::HEPEVT_Wrapper::check_hepevt_consistency(std::cout))
 //    HepMC::HEPEVT_Wrapper::print_hepevt(std::cout);

 if (cfg.fTest){
   // Test mode : compute directly some observables
   if(!HepMC::HEPEVT_Wrapper::check_hepevt_consistency(std::cout))
     HepMC::HEPEVT_Wrapper::print_hepevt(std::cout);
   double Mul2=0.;
   double Pla2=0.;
   double egy2=0.;
   for (HepMC::GenEvent::particle_const_iterator p = fEvtHepMC->particles_begin(); p != fEvtHepMC->particles_end(); ++p)
     {
       if( (*p)->status()==1){
         ++Mul2;
         double eta=(*p)->momentum().pseudoRapidity();
         double pt=(*p)->momentum().perp();
         if(fabs(eta)<10000){
           MeanPseudorapidity=MeanPseudorapidity+eta;
           ErrorMeanPseudorapidity=ErrorMeanPseudorapidity+eta*eta;}
         MeanPt=MeanPt+pt;
         ErrorMeanPt=ErrorMeanPt+pt*pt;
         if( fabs(eta) < 0.5 )++Pla2;
         egy2=egy2+(*p)->momentum().e();
       }
     }

   TotalEnergy=TotalEnergy+egy2;
   ErrorTotalEnergy=ErrorTotalEnergy+egy2*egy2;
   Multiplicity=Multiplicity+Mul2;
   ErrorMultiplicity=ErrorMultiplicity+Mul2*Mul2;
   PlateauHeight=PlateauHeight+Pla2;
   ErrorPlateauHeight=ErrorPlateauHeight+Pla2*Pla2;

 }else{

#ifdef HEPMC_HAS_CROSS_SECTION
  // set cross section information for this event
  HepMC::GenCrossSection theCrossSection;
  if(cfg.fProjectileId>1 || cfg.fTargetId>1){
  theCrossSection.set_cross_section(double(gCRMC_data.sigineaa)*1e9); //required in pB
  }else{
  theCrossSection.set_cross_section(double(gCRMC_data.sigine)*1e9); //required in pB
  }
  fEvtHepMC->set_cross_section(theCrossSection);
#endif

  // provide optional pdf set id numbers for CMSSW to work
  // flavour of partons and stuff. hope it's optional
  HepMC::PdfInfo pdf(0, 0, 0, 0, 0, 0, 0);
  fEvtHepMC->set_pdf_info(pdf);

  //Setting heavy ion infromation
  //      int   Ncoll_hard          // Number of hard scatterings
  //      int   Npart_proj          // Number of projectile participants
  //      int   Npart_targ          // Number of target participants
  //      int   Ncoll               // Number of NN (nucleon-nucleon) collisions
  //      int   spectator_neutrons           // Number of spectator neutrons
  //      int   spectator_protons            // Number of spectator protons
  //      int   N_Nwounded_collisions        // Number of N-Nwounded collisions (here Glauber number of participants with at least 1 interaction)
  //      int   Nwounded_N_collisions        // Number of Nwounded-N collisons (here Glauber number of participants with at least 2 interaction2)
  //      int   Nwounded_Nwounded_collisions // Number of Nwounded-Nwounded collisions (here GLauber number of collisions)
  //      float impact_parameter        // Impact Parameter(fm) of collision
  //      float event_plane_angle       // Azimuthal angle of event plane
  //      float eccentricity            // eccentricity of participating nucleons
  //                                        in the transverse plane
  //                                        (as in phobos nucl-ex/0510031)
  //      float sigma_inel_NN           // nucleon-nucleon inelastic
  //                                        (including diffractive) cross-section

  HepMC::HeavyIon ion(gCRMC_data.kohevt,
		      gCRMC_data.npjevt,
		      gCRMC_data.ntgevt,
		      gCRMC_data.kolevt,
		      gCRMC_data.npnevt + gCRMC_data.ntnevt,
		      gCRMC_data.nppevt + gCRMC_data.ntpevt,
		      gCRMC_data.ng1evt,
		      gCRMC_data.ng2evt,
		      gCRMC_data.nglevt,
		      gCRMC_data.bimevt,
		      gCRMC_data.phievt,
		      gCRMC_data.fglevt,  //defined only if phimin=phimax=0.
		      gCRMC_data.sigine*1e9); //required in pB
  fEvtHepMC->set_heavy_ion(ion);

  // add some information to the event
  fEvtHepMC->set_event_number(nEvent);

  //an integer ID uniquely specifying the signal process (i.e. MSUB in Pythia)
  int sig_id = -1;
  switch (gCRMC_data.typevt) // if negative typevt mini plasma was created by event (except -4)
    {
    case  0: break; //unknown for qgsjetII
    case  1: sig_id = 101; break;
    case -1: sig_id = 101; break;
    case  2: sig_id = 105; break;
    case -2: sig_id = 105; break;
    case  3: sig_id = 106; break;
    case -3: sig_id = 106; break;
    case  4: sig_id = 103; break;
    case -4: sig_id = 104; break;
    default: cerr << "Signal ID not recognised for setting HEPEVT" << endl;
    }
  fEvtHepMC->set_signal_process_id(sig_id);

  if(fEvtHepMC->vertices_begin()!=fEvtHepMC->vertices_end())
    fEvtHepMC->set_signal_process_vertex(*(fEvtHepMC->vertices_begin()));

  //DEBUG OUTPUT
  /*  for (HepMC::GenEvent::particle_const_iterator p = fEvtHepMC->particles_begin(); p != fEvtHepMC->particles_end(); ++p)
     {
       (*p)->print();
     }
  for (HepMC::GenEvent::vertex_const_iterator v = fEvtHepMC->vertices_begin(); v != fEvtHepMC->vertices_end(); ++v)
     {
       (*v)->print();
     }
  */

  // write the event out to the ascii file
  (*ascii_out) << fEvtHepMC;
 }
  // we also need to delete the created event from memory
 delete fEvtHepMC;
}


void
OutputPolicyHepMC::CloseOutput(const CRMCoptions& cfg)
{
  if(cfg.fTest) PrintTestEvent(cfg);
  //fOut->close();
  delete hepevtio;
  delete ascii_out;
  delete fOut;
}

void
OutputPolicyHepMC::PrintTestEvent(const CRMCoptions& cfg)
{
  cout.setf(ios::showpoint);
  cout.setf(ios::fixed);
  cout.precision(3);

  if(Multiplicity > 0){

    ErrorMeanPseudorapidity=sqrt(max(0.,ErrorMeanPseudorapidity
                                     -MeanPseudorapidity*MeanPseudorapidity/Multiplicity)
                                 /(max(Multiplicity-1.,Multiplicity)*Multiplicity));
    MeanPseudorapidity=MeanPseudorapidity/Multiplicity;
    ErrorMeanPt=sqrt(max(0.,ErrorMeanPt
                                     -MeanPt*MeanPt/Multiplicity)
                     /(max(Multiplicity-1.,Multiplicity)*Multiplicity));
    MeanPt=MeanPt/Multiplicity;
    ErrorMultiplicity=sqrt(max(0.,ErrorMultiplicity
                                     -Multiplicity*Multiplicity/cfg.fNCollision)
                           /(max(cfg.fNCollision-1.,double(cfg.fNCollision))*cfg.fNCollision));
    Multiplicity=Multiplicity/cfg.fNCollision;
    ErrorPlateauHeight=sqrt(max(0.,ErrorPlateauHeight
                                     -PlateauHeight*PlateauHeight/cfg.fNCollision)
                            /(max(cfg.fNCollision-1.,double(cfg.fNCollision))*cfg.fNCollision));
    PlateauHeight=PlateauHeight/cfg.fNCollision;
    ErrorTotalEnergy=sqrt(max(0.,ErrorTotalEnergy
                                     -TotalEnergy*TotalEnergy/cfg.fNCollision)
                            /(max(cfg.fNCollision-1.,double(cfg.fNCollision))*cfg.fNCollision));
    TotalEnergy=TotalEnergy/cfg.fNCollision;


    cout << "\n          >> Test output <<\n\n"
         << "  Total Cross Section (mb):    " << gCRMC_data.sigtot << "\n"
         << "  Elastic Cross Section (mb):  " << gCRMC_data.sigela << "\n"
         << "  Inel. Cross Section (mb) :   " << gCRMC_data.sigine << "\n" ;
    if(cfg.fProjectileId>1 || cfg.fTargetId>1)
    cout << "  Inel. AA Cross Section (mb): " << gCRMC_data.sigineaa << "\n"
         << "  Elastic AA Cross Section (mb): " << gCRMC_data.sigelaaa << "\n"
         << "  Total AA Cross Section (mb): " << gCRMC_data.sigtotaa << "\n" ;

    cout << endl;
    cout << "  Energy (GeV):                " << TotalEnergy
         <<                           " +/- " << ErrorTotalEnergy << "\n"
         << "  Multiplicity:                " << Multiplicity
         <<                           " +/- " << ErrorMultiplicity << "\n"
         << "  PlateauHeight:               " << PlateauHeight
         <<                           " +/- " << ErrorPlateauHeight << "\n"
         << "  MeanPseudorapidity:          " << MeanPseudorapidity
         <<                           " +/- " << ErrorMeanPseudorapidity << "\n"
         << "  MeanPt (GeV/c):              " << MeanPt
         <<                           " +/- " << ErrorMeanPt << "\n";
    cout << endl;

  }else{

    cout << "Error during test : no particles !" << endl;

  }

}
