#include <CRMCinterface.h>
#include <CRMCconfig.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
using namespace std;

CRMCdata gCRMC_data;

CRMCinterface::CRMCinterface() :
  crmc_generate(NULL),
  crmc_set(NULL),
  crmc_init(NULL),
  crmc_xsection(NULL),
  fLibrary(NULL)
{
}

CRMCinterface::~CRMCinterface()
{
  if (fLibrary)
    {
      dlclose(fLibrary);
      fLibrary = NULL;
    }
}

bool CRMCinterface::init(int HEmodel)
{
#ifdef __CRMCSTATIC__
  crmc_generate = &crmc_f_;
  crmc_set = &crmc_set_f_;
  crmc_init = &crmc_init_f_;
  crmc_xsection = &crmc_xsection_f_;
#else
  ostringstream libname;
  if (!fLibrary)
    {
      //rpath is being checked for non absolute paths
      libname << "lib";
      switch (HEmodel)
        {
        case 0: libname << "Epos"; break;
        case 1: libname << "Epos"; break;
        case 2: libname << "Qgsjet01"; break;
        case 3: libname << "Gheisha"; break;
        case 4: libname << "Pythia"; break;
        case 5: libname << "Hijing"; break;
        case 6: libname << "Sibyll"; break;
#ifdef __QGSJETII04__
        case 7: libname << "QgsjetII04"; break;
#endif
        case 8: libname << "Phojet"; break;
#ifdef __QGSJETII03__
        case 11: libname << "QgsjetII03"; break;
#endif
        case 12: libname << "Dpmjet"; break;
        default: libname << "UnknownModel"; break;
        }
      libname << ".so";
      fLibrary = dlopen(libname.str().c_str(), RTLD_NOW);
      cout << "Opening: " << libname.str() << endl;
      if (!fLibrary )
        {
          ostringstream errMsg;
          errMsg << "\n cannot open shared library " << libname.str() << "\'\n\n"
                 << " Dynamic-link error:\n \"" << dlerror() << "\"\n";

          cerr << errMsg.str() << endl;
          exit(1);
        }
    }

  crmc_generate = (void(*)( const int&, const int&, int&, double&, int&, double&,
                      double&, double&, double&, double&, int&)) dlsym(fLibrary, "crmc_f_");
  if(crmc_generate == NULL)
    {
      ostringstream errMsg;
      errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

      cerr << errMsg.str() << endl;
      exit(1);
    }

  crmc_set = (void(*)( const int&, const int&, const double&, const double&,
                           const int&, const int&, const int&, const int&,
                           const int&, const char*)) dlsym(fLibrary, "crmc_set_f_");
  if(crmc_set == NULL)
    {
      ostringstream errMsg;
      errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

      cerr << errMsg.str() << endl;
      exit(1);
    }

  crmc_init = (void(*)(const char*,const int&)) dlsym(fLibrary, "crmc_init_f_");
  if(crmc_init == NULL)
    {
      ostringstream errMsg;
      errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

      cerr << errMsg.str() << endl;
      exit(1);
    }
  
  crmc_xsection = (void(*)( double&, double&, double&, double&, double&,
                              double&, double&, double&, double&)) dlsym(fLibrary, "crmc_xsection_f_");
  if(crmc_xsection == NULL)
    {
      ostringstream errMsg;
      errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

      cerr << errMsg.str() << endl;
      exit(1);
    }

  //common blocks from library are not used. they come from DummyHepEvt library
  //grabbed in header with extern "C"
  // cevt_  = (cevt*)  dlsym(fLibrary, "cevt_");
  // if(cevt_ == NULL)
  //   {
  //     ostringstream errMsg;
  //     errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

  //     cerr << errMsg.str() << endl;
  //     exit(1);
  //   }
  // c2evt_ = (c2evt*) dlsym(fLibrary, "c2evt_");
  // if(c2evt_ == NULL)
  //   {
  //     ostringstream errMsg;
  //     errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

  //     cerr << errMsg.str() << endl;
  //     exit(1);
  //   }
  // hadr5_ = (hadr5*) dlsym(fLibrary, "hadr5_");
  // if(hadr5_ == NULL)
  //   {
  //     ostringstream errMsg;
  //     errMsg << " dlsym error:\n \"" << dlerror() << "\"\n";

  //     cerr << errMsg.str() << endl;
  //     exit(1);
  //   }

#endif
  return 1;
}
