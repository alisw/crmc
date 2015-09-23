#ifndef __analysis__h__
#define __analysis__h__

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <exception>

#include "HepMC/IO_GenEvent.h"

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>

// #include <TParticle.h> //nice class for more information about particles
// #include <TParticlePDG.h>

// TParticle* CopyHepMC2ROOT(HepMC::GenParticle* a);

class DataManager
{
private:
  std::vector<std::string> filelist;
  std::vector<std::string>::iterator current_file;

  boost::iostreams::filtering_istream in;
  HepMC::IO_GenEvent* ascii_in;
  int count;


  bool GetNextFile()
  {
    if(++current_file == filelist.end())
      return 0;
    return ReadFile();
  }

  bool ReadFile()
  {
    std::string filename = *current_file;
    bool use_compression = false;
    if((filename).find(".gz") != std::string::npos) use_compression = true;
    std::cout << "DataMangager::Opening file " << (filename) << " with" << (use_compression?" ":"out ") << "compression. " << std::endl;

    in.reset();
    if(use_compression)
      in.push (boost::iostreams::gzip_decompressor ());
    in.push (boost::iostreams::file_descriptor_source (filename));
    if(ascii_in) delete ascii_in;
    ascii_in = new HepMC::IO_GenEvent(in);
    return true;
  }

public:
  HepMC::GenEvent* evt;
  DataManager() : ascii_in(0), count(0), evt(0)
  {

  }
  ~DataManager()
  {
    if(ascii_in) delete ascii_in;
    std::cout << "DataManger::Closing after reading " << count << " events." << std::endl;
  }

  void SetFiles(const std::vector<std::string>& files)
  {
    filelist = files; //copy
    //open first file
    current_file = filelist.begin();
    ReadFile();
  }


  bool GetNextEvent()
  {
    delete evt;
    evt = 0;
    try
      {
        (*ascii_in) >> evt;
      }
    catch (std::exception& e)
      {
        evt = NULL;
        std::cout << "Event could not be read: " << e.what() << std::endl;
      }
    if (!evt)
      {
	if (!GetNextFile()) return false;
	else // new file selected
	  {
            delete evt;
	    evt = NULL;
	    try
              {
                (*ascii_in) >> evt;
              }
            catch(std::exception& e)
              {
                evt = NULL;
                std::cout << "Event could not be read again: " << e.what() << std::endl;
              }
	    if (!evt)  return false;
	  }
      }

    if(count % 1000 == 0)
      std::cout << "Reading event: " << count << std::endl;
    count++;
    return true;
  }

};

// TParticle* CopyHepMC2ROOT(HepMC::GenParticle* a)
// {
//   TParticle* b = new TParticle;
//   if(!b)
//     return 0;
//   b->SetPdgCode(a->pdg_id ());

//   b->SetStatusCode(a->status ());

//   b->SetMomentum(a->momentum ().px (),
// 		 a->momentum ().py (),
// 		 a->momentum ().pz (),
// 		 a->momentum ().e ()
// 		 );

//   b->SetProductionVertex(a->production_vertex ()->position ().x (),
// 			 a->production_vertex ()->position ().y (),
// 			 a->production_vertex ()->position ().z (),
// 			 a->production_vertex ()->position ().t ()
// 			 );
//   return b;
// }

#endif
