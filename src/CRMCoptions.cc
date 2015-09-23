#include <CRMCoptions.h>
#include <CRMCconfig.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

CRMCoptions::CRMCoptions(int argc, char** argv)
  : fError(false),
    fOutputMode(eHepMC),
    fNCollision(500),
    fSeed(0),
    fProjectileId(1),
    fTargetId(1),
    fHEModel(0),
    fTypout(0),
    fProjectileMomentum(3500),
    fTargetMomentum(-3500),
    fParamFileName("crmc.param "),
    fOutputFileName(""),
    fProduceTables(false),
    fSeedProvided(false),
    fTest(false),
    fCSMode(false)
{
  CheckEnvironment();
  ParseOptions(argc, argv);
}


void
CRMCoptions::CheckEnvironment()
{
  const char* crmcRoot = getenv("CONEX_ROOT");
  if (crmcRoot == 0)
    {
      crmcRoot = getenv("PWD");
    }

  /*
  const fs::path basepath(argv[0]);
  iUtl::Config::Init(basepath.parent_path() / "../..", argc, argv);

    ptime now = second_clock::local_time();
    date today = now.date();
    const date ignoreUntil = today - days(deltaDay) - months(deltaMonth) - years(deltaYear);
    fIgnoreUntilDay = ignoreUntil.day();
    fIgnoreUntilMonth = ignoreUntil.month();
    fIgnoreUntilYear = ignoreUntil.year();


  */


}


void
CRMCoptions::ParseOptions(int argc, char** argv)
{
  po::options_description desc("Options of CRMC");

  ostringstream model_desc;
  model_desc << "model [0=EPOS_LHC, 1=EPOS_1.99"
#ifdef __QGSJET01__
	     << ", 2=QGSJET01"
#endif
#ifdef __GHEISHA__
	     << ", 3=Gheisha"
#endif
#ifdef __PYTHIA__
	     << ", 4=Pythia_6.115"
#endif
#ifdef __HIJING__
	     << ", 5=Hijing_1.38"
#endif
#ifdef __SIBYLL__
	     << ", 6=Sibyll_2.1"
#endif
#ifdef __QGSJETII04__
	     << ", 7=QGSJETII-04"
#endif
#ifdef __PHOJET__
	     << ", 8=Phojet"
#endif
#ifdef __QGSJETII03__
	     << ", 11=QGSJETII-03"
#endif
#ifdef __DPMJET__
	     << ", 12=DPMJet 3.0-6"
#endif
	     << "]";

  desc.add_options()
    ("help,h", "description of options")
    ("version,v", "show version and exits")
    ("output,o", po::value<string>(), "output mode: hepmc (default), hepmcgz, root, lhe, lhegz")
    //    ("verbosity,v", po::value<int>(), "verbosity (default: 0)")
    ("seed,s", po::value<int>(), "random seed between 0 and 1e9 (default: random)")
    ("number,n", po::value<int>(), "number of collisions")
    ("model,m", po::value<int>(), model_desc.str().c_str())
    ("projectile-momentum,p", po::value<double>(), "momentum/(GeV/c)")
    ("target-momentum,P", po::value<double>(), "momentum/(GeV/c)")
    ("projectile-id,i", po::value<int>(), "PDG or Z*10000+A*10")
    ("target-id,I", po::value<int>(), "PDG or Z*10000+A*10")
    ("config,c", po::value<string>(), "config file")
    ("out,f", po::value<string>(), "output file name (auto if none provided)")
    ("produce-tables,t", po::value<bool>()->implicit_value(1), "create tables if none are found")
    //("filter,F", po::value<string>(), "specify file with filter commands")
    ("test,T", po::value<bool>()->implicit_value(1), "test mode")
    ("cross-section,x", po::value<bool>()->implicit_value(1), "calculate and print cross section only")
    ;

  po::variables_map opt;
  po::store(po::parse_command_line(argc, argv, desc), opt);
  po::notify(opt);

  if (opt.count("version")) {
    cout << "crmc v" << CRMC_VERSION_MAJOR << "." << CRMC_VERSION_MINOR << endl;
    exit(0);
  }

  // check if user needs help-printout
  if (opt.count("help")) {
    cout << endl << desc << endl;
    exit(0);
  }

  if(opt.count("output")) {

    const string om = opt["output"].as<string>();
    if (om=="hepmc") fOutputMode = eHepMC;
    else if (om=="hepmcgz") fOutputMode = eHepMCGZ;
    else if (om=="lhe") fOutputMode = eLHE;
    else if (om=="lhegz") fOutputMode = eLHEGZ;
    else if (om=="root") {
#ifdef WITH_ROOT
      fOutputMode = eROOT;
#else
      cerr << " Compile with ROOT first " << endl;
      exit(1);
#endif
    } else {
      cerr << " Wrong output type: " << om << endl;
      cerr << endl << desc << endl;
      exit(1);
    }
  }

  if (fOutputMode == eLHE || fOutputMode == eLHEGZ)
    fTypout = 1;

#ifdef WITH_ROOT
  if (fOutputMode == eROOT)
    fTypout = -1;
#endif

  // parameter readout
  if (opt.count("seed")) {
    fSeed = opt["seed"].as<int>();
    if (fSeed < 0) {
      cerr << " Seed is negative: " << fSeed << endl;
      cerr << endl << desc << endl;
      exit(1);
    }
    if (fSeed > 1e9) {
      cerr << " Seed too large (>1e9): " << fSeed << endl;
      cerr << endl << desc << endl;
      exit(1);
    }
  }

  if (opt.count("number"))
    fNCollision = opt["number"].as<int>();

  if (opt.count("model"))
    fHEModel = opt["model"].as<int>();

  if (opt.count("projectile-momentum"))
    fProjectileMomentum = opt["projectile-momentum"].as<double>();

  if (opt.count("target-momentum"))
    fTargetMomentum = opt["target-momentum"].as<double>();

  if (opt.count("projectile-id"))
    fProjectileId = opt["projectile-id"].as<int>();

  if (opt.count("target-id"))
    fTargetId = opt["target-id"].as<int>();

  if (opt.count("config")) {
    fParamFileName = opt["config"].as<string>();
    fParamFileName += ' ';      //space needed at the end for Fortran command "index". "trim" needs f90 therefore not used here
}

  if (opt.count("out"))
    fOutputFileName = opt["out"].as<string>();

  if (opt.count("produce-tables")) {
    fProduceTables = opt["produce-tables"].as<bool>();
  }

  //if (opt.count("filter")) {
  //fFilter = opt["filter"].as<string>();
  //}

  if (opt.count("test"))
    {
      fOutputMode = eHepMC;
      fTest = opt["test"].as<bool>();
    }

  if (opt.count("cross-section"))
    {
      fCSMode = opt["cross-section"].as<bool>();
      fOutputMode = eNone;
      fNCollision = 1;
    }

  // check if random seed was provided, otherwise generate one
  fSeedProvided = fSeed;
  if (!fSeedProvided) {
    if (fTest)
      fSeed = 123; //if test and no seed provided fix to get same results
    else {
      ifstream urandom("/dev/urandom", ios::in|ios::binary);
      urandom.read((char*)&fSeed, sizeof(fSeed)/sizeof(char));
      urandom.close();
      fSeed = abs(fSeed) % 999999999;
    }
  }

  DumpConfig();

}

string
CRMCoptions::ParticleName(const int pid)
const
{
  if (pid < 10000) {

    switch (pid) {
    case 120  : return "pi"; break;
    case -120 : return "antipi"; break;
    case 1    : return "p";  break;
    case -1   : return "antip";  break;
    case 12   : return "C";   break;
    case 208  : return "Pb";  break;
    default:
      {
	ostringstream ss;
	ss << "pdg" << pid;
	return ss.str();
	break;
      }
    }
  }
  const int Z =  pid/10000;
  const int A = (pid%10000) / 10;

  switch(Z) {
  case 1: return "p"; break;
  case 2: return "He"; break;
  case 6: return "C"; break;
  case 7: return "N"; break;
  case 8: return "O"; break;
  case 26: return "Fe"; break;
  case 82: return "Pb"; break;
  }
  ostringstream ss;
  ss << "A" << A << "Z" << Z;
  return ss.str();
}

void
CRMCoptions::DumpConfig() const
{
  //if particles have same momentum the minus sign was probably forgotten
  if(fTargetMomentum == fProjectileMomentum)
    {
      cerr << " Beam particles at rest in the centre-of-mass frame. Did you forget a (-) sign for the target momentum? " << endl;
      cerr << "          exit ..." << endl;
      exit(1);
    }

 if ( fTest ) {
  cout << "\n          >> crmc Test mode <<\n\n" << endl;
 }else{
  cout << "\n          >> crmc <<\n\n" ;
 }
  cout << "  seed:                       " << fSeed << (fSeedProvided ? " (provided by user)" : " (automatic)") << "\n"
       << "  projectile id:              " << fProjectileId;
  if (fProjectileId/10000>0) {
    const int Z =  fProjectileId/10000;
    const int A = (fProjectileId%10000) / 10;
    cout << " (A=" << A << ", Z=" << Z << ")";
  }
  cout << "\n"
       << "  projectile momentum:        " << fProjectileMomentum << "\n"
       << "  target id:                  " << fTargetId;
  if (fTargetId/10000>0) {
    const int Z =  fTargetId/10000;
    const int A = (fTargetId%10000) / 10;
    cout << " (A=" << A << ", Z=" << Z <<")";
  }
  cout << "\n"
       << "  target momentum:            " << fTargetMomentum << "\n\n";

  cout << "  number of collisions:       " << fNCollision << "\n"
       << "  parameter file name:        " << fParamFileName << "\n";
  if (!fTest && !fCSMode)
    cout << "  output file name:           " << GetOutputFileName() << "\n";
  cout << "  HE model:                   " << fHEModel;

  switch(fHEModel) {
  case 0: cout << " (EPOS-LHC) \n"; break;
  case 1: cout << " (EPOS 1.99) \n"; break;
  case 2: cout << " (QGSJET01) \n"; break;
  case 3: cout << " (Gheisha)\n "; break;
  case 4: cout << " (Pythia)\n "; break;
  case 5: cout << " (Hijing)\n "; break;
  case 6: cout << " (Sibyll 2.1)\n "; break;
  case 7: cout << " (QGSJETII-04) \n"; break;
  case 8: cout << " (Phojet) \n"; break;
  case 11: cout << " (QGSJETII-03) \n"; break;
  case 12: cout << " (DPMJet 3.0-6) \n"; break;
  default:
    cerr << " (unknown model) \n";
    exit(1);
  }
  cout << endl;


  cout.setf(ios::showpoint);
  cout.setf(ios::fixed);
  cout.precision(3);

}


string
CRMCoptions::GetOutputTypeEnding() const
{
  switch (fOutputMode) {
  case eHepMC:
    return ".hepmc";
    break;
  case eHepMCGZ:
    return ".hepmc.gz";
    break;
  case eLHE:
    return ".lhe";
    break;
  case eLHEGZ:
    return ".lhe.gz";
    break;
#ifdef WITH_ROOT
  case eROOT:
    return ".root";
    break;
#endif
  }
  return ".unknown";
}

string
CRMCoptions::GetOutputFileName() const
{
  // open output file and connect tree

  const char* crmcOutDir = getenv("CRMC_OUT");
  if (crmcOutDir == 0 ) {
    crmcOutDir = getenv("PWD");
  }

  //check compatibility of path with lhe files
  if ( GetOutputTypeEnding().find("lhe") != string::npos && string(crmcOutDir).find(".lhe") != string::npos ) {
    cerr << " path error - path contains '.lhe' : " << crmcOutDir << endl;
    cerr << "          exit ..." << endl;
    exit(1);
  }

  //if file name provided
  if (fOutputFileName != "") {
  //check compatibility of file name with lhe files
  if ( GetOutputTypeEnding().find("lhe") != string::npos && fOutputFileName.find(".lhe") == string::npos ) {
    cerr << " file name error - extension is not '.lhe' : " << fOutputFileName << endl;
    cerr << "          exit ..." << endl;
    exit(1);
  }{
    return fOutputFileName;
  }}

  //create file name based on options
  ostringstream crmcFileName;
  crmcFileName << crmcOutDir << "/" << "crmc_";
  switch (fHEModel) {
  case 0: crmcFileName << "eposlhc";   break;
  case 1: crmcFileName << "epos199";      break;
  case 2: crmcFileName << "qgsjet";    break;
  case 3: crmcFileName << "gheisha";   break;
  case 4: crmcFileName << "pythia";    break;
  case 5: crmcFileName << "hijing";    break;
  case 6: crmcFileName << "sibyll";    break;
  case 7: crmcFileName << "qgsjetII04";  break;
  case 8: crmcFileName << "phojet";    break;
  case 11:crmcFileName << "qgsjetII03";  break;
  case 12: crmcFileName << "dpmjet";    break;
  default:
    cerr << " crmcOut: error - unknown model " << fHEModel << endl;
    cerr << "          exit ..." << endl;
    exit(1);
    break;
  }

  crmcFileName << "_" << fSeed << "_";

  crmcFileName << ParticleName(fProjectileId) << "_";
  crmcFileName << ParticleName(fTargetId);

  crmcFileName << "_" << fProjectileMomentum
	       << GetOutputTypeEnding();

return crmcFileName.str();
}
