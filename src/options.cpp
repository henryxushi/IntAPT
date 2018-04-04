/*
 *  common.h
 */
#include "options.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <libgen.h>
#include <iostream>
#include <sstream>

string OPT_bamfile = "NaN";
string OPT_InstFile = "NaN";
string OPT_outdir = "NaN";
string OPT_readtype = "NaN";
int OPT_Ncores = 1;
double OPT_conf = 0.5;
bool OPT_help = false;
bool OPT_InstOnly = false;
string OPT_cempath = "";
string OPT_samtoolspath = "";
bool OPT_output_all = false;

void Options::assign(string s1, string s2, string s3, int s4, double s5, string s6, bool s7, string s8, bool s10, bool s11)
{
  bamfile = s1;
  outdir = s2;
  readtype = s3;
  N_cores = s4;
  conf = s5;
  InstFile = s6;
  InstOnly = s7;
  cempath = s8;
  samtoolspath = s8;
  output_all = s10;
  use_inst = s11;
}

int Options::parse_options(int argc, char* argv[]) {

  int option_index = 0;
  int next_option;

  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
    case -1:     /* Done with options. */
      break;
    case 'b':
      OPT_bamfile = optarg;
      break;
    case 'o':
      OPT_outdir = optarg;
      break;
    case 'r':
      OPT_readtype = optarg;
      break;
    case 'h':
      OPT_help = true;
      break;
    case 'p':
      OPT_Ncores = atoi(optarg);
      break;
    case 'c':
      OPT_conf = atof(optarg);
      break;
    case OPT_INSTFILE:
      OPT_InstFile = optarg;
      break;
    case OPT_CEMPATH:
      OPT_cempath = optarg;
      OPT_samtoolspath = optarg;
      if (!OPT_cempath.empty())
      {
        if (OPT_cempath[OPT_cempath.size()-1] != '/')
        {
          stringstream temp;
          temp << OPT_cempath << "/";
          OPT_cempath = temp.str();
          OPT_samtoolspath = temp.str();
        }
      }
      break;
    case OPT_OUTPUT_ALL:
      OPT_output_all = true;
      break;
    case OPT_INSTONLY:
      OPT_InstOnly = true;
      break;
    default:
      std::cout << usage();
      exit(1);
    }
  } while(next_option != -1);

  if (OPT_help) {
    std::cout << usage() ;
    exit(1);
  }
  
  //check input argument
  bool use_inst_bool = false;
  bool baminvalid = OPT_bamfile == "" || OPT_bamfile == "NaN";
  bool instfileinvalid = OPT_InstFile == "" || OPT_InstFile == "NaN";
  if (OPT_InstOnly)
  {
    if (baminvalid)
    {
      std::cerr << "You want to generate Inst file, please check the input bam file" << std::endl;
      std::cout << usage();
      exit(1);
    }
  }
  if (baminvalid and instfileinvalid) {
  	std::cerr << "Please check bam or inst file option\n";
    std::cout << usage() ;
    exit(1);
  }
  else if (!baminvalid and !instfileinvalid) {
    std::cout << "Found both bamlist and inst file, use inst file...\n";
    use_inst_bool = true;
  }
  else if (baminvalid and !instfileinvalid)
  {
    std::cout << "Use inst file...\n";
    use_inst_bool = true;
  }
  else if (!baminvalid and instfileinvalid)
  {
    std::cout << "Use BAM file...\n";
  }
  else
  {
    std::cerr << "Please check the input bam file" << std::endl;
    std::cout << usage();
    exit(1);
  }
  
  if (OPT_readtype != "p" && OPT_readtype != "s")
  {
  	std::cerr << "Please check read type option\n";
  	std::cout << usage() ;
    exit(1);
  }
  
  if (OPT_outdir == "" || OPT_outdir == "NaN") {
    std::cerr << "Please check output file directory option\n";
    std::cout << usage() ;
    exit(1);
  }

  if (OPT_Ncores < 0)
  {
  	std::cerr << "Please check threads option\n";
  	std::cout << usage() ;
    exit(1);
  }
  
  if (OPT_conf < 0 || OPT_conf > 1)
  {
  	std::cerr << "Please check confidence option\n";
  	std::cout << usage() ;
    exit(1);
  }

  assign(OPT_bamfile, OPT_outdir, OPT_readtype, OPT_Ncores, OPT_conf, OPT_InstFile, OPT_InstOnly, OPT_cempath, OPT_output_all, use_inst_bool);

  return 0;
}


string Options::usage () {

  std::stringstream usage_info;
  usage_info
    << std::endl
    << "===============================================================================" << std::endl
    << " Usage: IntAPT [--bam/b] <filename>  [opts] " << std::endl
    << "===============================================================================" << std::endl
    << " **Required :" << std::endl
    << " --bam/-b <string>         " << ": start from the list of the names of bam file " << std::endl
    << " --inst <string>           " << ": start from the processed inst file by IntAPT " << std::endl
    << " --out-dir/-o <string>     " << ": the directory stored all output files" << std::endl
    << " --readtype/-r <p/s>       " << ": the type of reads paired-end(p) or single-end(s)" << std::endl;

  usage_info
    << std::endl
    << " ** Prerequisite program (not needed to set if the tools in system path) :" << std::endl
    << " --IntAPTpath <string>     " << ": the path to the tools folder in the package e.g. $PATHTOINTAPTpath/tools." << std::endl;

  usage_info
    << std::endl
    << " ** Optional :" << std::endl
    << " --threads/-p <int>        " << ": the number of threads to be used, default: 1." << std::endl
    << " --outputall               " << ": output all isoforms enumberated from the graph" << std::endl
    << " --instonly                " << ": only output inst file (not compatible with --inst), default: false. " << std::endl
    << " --conf/-c <double>        " << ": the confidence level for isoform identification, default: 0.5." << std::endl
    << " --help/-h                 " << ": display the help information."<< std::endl
    << std::endl;
    
  usage_info
    << "================================================================================" << std::endl
    << std::endl;

  return (usage_info.str());
}

