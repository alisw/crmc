#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdlib.h>


#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>


#include "HepMC/IO_HEPEVT.h"
#include "HepMC/GenParticle.h"
#include <HepMC/HeavyIon.h>
#include "HepMC/HEPEVT_Wrapper.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"
#include "HepMC/GenEvent.h"


#include <TFile.h>
#include <TH1D.h>


#include "analysis.h"


namespace io = boost::iostreams;
using namespace std;




int main (int argc, char **argv)
{
  //-------------------SET UP DATA
  TFile* theOutFile;
  string outFileName ("new_histogram_file.root");
  cout << " ! Opening output file: " << outFileName << endl;
  theOutFile = new TFile (outFileName.c_str(),"RECREATE");

  vector<string> filesModel1;
  filesModel1.push_back ("files/test.hepmc"); //add your files here with additional push_back()
  DataManager data;
  data.SetFiles (filesModel1); //for more models, loop over models and call SetFiles each time
  theOutFile->mkdir ("model1");

  //------------------SET UP HISTOGRAMS
  TH1D* exampleHist = new TH1D ("dNdeta",";#eta;dN/d#eta",21,-10,10);


  //-------------------EVENT LOOP
  int nEvts = 0;
  while (data.GetNextEvent ())
    {
      ++nEvts;
      //HepMC::HeavyIon* heavyIonInfo = NULL; //helpful to get cross section: heavyIonInfo->sigma_inel_NN ()
      //heavyIonInfo = data.evt->heavy_ion ();

      //-------------------PARTICLE LOOP
      HepMC::GenEvent::particle_const_iterator par = data.evt->particles_begin ();
      for (; par != data.evt->particles_end (); ++par)
	{
	  HepMC::GenParticle* p = (*par);

	  if (p->status () != 1) continue; //get final state particles. status == 2 are decayed particles, status == 4 is beam particles

	  //HepMC::GenVertex* parent_vertex = p->production_vertex();
	  //const int id = p->pdg_id ();
	  const double eta = p->momentum ().eta ();
	  //const double pt = p->momentum ().perp ();
	  //const double e = p->momentum ().e ();

          //for more advance paramters see HepMC documentation or load #include <TParticle.h> and fill object (see analysis.h)

          //-------------------EVENT SELECTION
          //-------------------FILL HISTOGRAMS WITH PER PARTICLE VARIABLES
          exampleHist->Fill (eta);
	}//PARTICLE LOOP
      //-------------------FILL HISTOGRAMS WITH PER EVENT VARIABLES
    }//EVENT LOOP
    //---------------FINALISE HISTOGRAMS
  exampleHist->Scale (1. / nEvts, "width");

  //----------------Closing Files
  std::cout << " ! Writing output file: " << outFileName << std::endl;
  theOutFile->Write();
  delete exampleHist;
  exampleHist = 0;
  std::cout << " ! Closing output file: " << outFileName << std::endl;
  theOutFile->Close();
  delete theOutFile;
  theOutFile = 0;
  return 0;
}
