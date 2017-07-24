/////////////////////////////////////////////////////////////////////////
// 
// author: Alejandro Hernandez Wences
// date: 21 Feb 2013
// 
// Merge two assemblies
//
// see merge.cc and merge.hh
/////////////////////////////////////////////////////////////////////////

#include "merge.hh"
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string>
#include <map>
#include <vector>
#include <getopt.h>

using namespace std;

//-------------------------------------- Globals
//--Required options
char * OUT;
char * MUMCOORDS;
char * DELTAFILE;
char * RCE;
char * RNS;
char * QCE;
char * QNS;
//-- Metassembly options
bool PRINT_MUMS=true;
int MIN_MIS_OR=13;
double MIN_ZSTAT=3;
double MIN_DIFF_Z=2;
int MIN_COV=15;
int MIN_SUP=0;
string SCF_PREFIX="metaScf";
//-- N-insert options
int MIN_N_INS=50;
double MIN_N_PER=0.65;
//-- Linking options
int MIN_LINK_COV=20;
int MIN_ALN_LEN=60;
int MIN_BASES_ALN=60;
int MAX_EDGE_UNALN=200;

//-------------------------------------- Function declarations

ostream& printParameters(ostream& out);
bool overlap(int a, int b, int c, int d);
void ParseArgs (int argc, char ** argv);
void PrintHelp (const char * s);
void PrintUsage (const char * s);

int main(int argc, char** argv){

  ParseArgs(argc,argv);
  printParameters(cout); 
  ofstream OUTPA (string(OUT).append(".parameters").c_str());
  printParameters(OUTPA);

  std::cout<<"\n\n---- Loading data"<<std::endl
           <<"---------------------------------"<<std::endl;

  Metassembly Meta(string(OUT),
                   SCF_PREFIX,
                   DELTAFILE,
                   MUMCOORDS,
                   MIN_ALN_LEN,
                   MIN_BASES_ALN,
                   MAX_EDGE_UNALN,
                   RCE,
                   QCE,
                   MIN_MIS_OR,
                   MIN_SUP,
                   RNS,
                   QNS,
                   MIN_COV,
                   MIN_ZSTAT,
                   MIN_DIFF_Z,
                   MIN_LINK_COV,
                   MIN_N_INS,
                   MIN_N_PER);

  if(PRINT_MUMS)
    Meta.mum.printMums(string(OUT));

  int metaScfN=0;
  cout<<"\n\n---- Computing Metassembly" <<std::endl
      <<"---------------------------------"<<endl;

  map<string,vector<int> >::iterator refMUMaln;
  for(refMUMaln=Meta.mum.refScfs.begin();
      refMUMaln != Meta.mum.refScfs.end();
      refMUMaln ++ 
     ){

    string refScf=(*refMUMaln).first;
    vector<int>& refAlns=(*refMUMaln).second;
    metaScfN+=1;
    string metaScf("metaScf_");
    std::ostringstream c;
    c << metaScfN;
    metaScf.append(c.str());
    Meta.Scf_Meta[refScf]=metaScf;
    
    //Analyse possible left and right links
    vector<int>::iterator leftAln=Meta.mum.findRefEdgeAln(refScf,"start");
    vector<int>::iterator rightAln=Meta.mum.findRefEdgeAln(refScf,"end");
    vector<int>::iterator beginScfAnalysis;
    vector<int>::iterator endScfAnalysis;
    bool leftLink=false;
    bool rightLink=false;
    if(leftAln != refAlns.end()){ // A non-spurious leftAln was found
      string scfLeftLink=Meta.findScfLink(leftAln,"start");
      if(scfLeftLink.empty())
        beginScfAnalysis=refAlns.begin();
      else{
        beginScfAnalysis=leftAln;
        leftLink=true;
      }
    }else
      beginScfAnalysis=refAlns.begin();

    if(rightAln != refAlns.end()){
      string scfRightLink=Meta.findScfLink(rightAln,"end");
      if(scfRightLink.empty())
        endScfAnalysis=refAlns.end()-1;
      else{
        endScfAnalysis=rightAln;
        rightLink=true;
      }
    }else
      endScfAnalysis=refAlns.end()-1;

    //Analyse and merge refScf with the corresponding secondary scaffolds
    vector<int>::iterator it;
    it=beginScfAnalysis;
    while(it <= endScfAnalysis){
      it=Meta.analyseAln(refScf,
                         refAlns,
                         it,
                         metaScf,
                         leftLink,
                         beginScfAnalysis,
                         rightLink,
                         endScfAnalysis
                        );
    }
  }

  Meta.printLinks();

  Meta.outp=string(OUT).append(".noLink");
  Meta.printMetassem();
  Meta.printStats();

  cout<<"\n\n---- Joining Links "<<std::endl
      <<"---------------------------------"<<std::endl;
  Meta.joinLinks();
  Meta.outp=string(OUT);
  cout<<"\n\n---- Printing Metassembly "<<std::endl
      <<"---------------------------------"<<std::endl;
  Meta.printMetassem();
  cout<<"\n\n---- Printing Stats "<<std::endl
      <<"---------------------------------"<<std::endl;
  Meta.printStats();

}

ostream& printParameters(ostream& OUTF){

  int i;

  OUTF << "\n\n-----asseMerge Parameters:"<<std::endl
       << "---------------------------------"<<std::endl;
  OUTF  << "-R "<<string(RCE)<<"\n"
        << "-Q "<<string(QCE)<<"\n"
        << "-r "<<string(RNS)<<"\n"
        << "-q "<<string(QNS)<<"\n"
        << "-D "<<string(DELTAFILE)<<"\n"
        << "-M "<<MUMCOORDS<<"\n"
        << "-O "<<string(OUT)<<"\n"
        << "Optional:\n"
        << "\n  Signal:\n"
        << "  -z "<<MIN_ZSTAT<<"\n"
        << "  -d "<<MIN_DIFF_Z<<"\n"
        << "  -i "<<MIN_MIS_OR<<"\n"
        << "  -s "<<MIN_SUP<<"\n"
        << "  -c "<<MIN_COV<<"\n"
        << "\n  Linking:\n"
        << "  -e "<<MAX_EDGE_UNALN<<"\n"
        << "  -l "<<MIN_ALN_LEN<<"\n"
        << "  -a "<<MIN_BASES_ALN<<"\n"
        << "  -L "<<MIN_LINK_COV<<"\n"
        << "\n  N-insert:\n"
        << "  -t "<<MIN_N_INS<<"\n"
        << "  -p "<<MIN_N_PER<<"\n"
        << "\n Other:"<<std::endl;
  if(!PRINT_MUMS)
    OUTF<< "  -C Do not print .coords file ordered by primary and secondary assembly sequence names\n";
  OUTF  << "  -x "<<SCF_PREFIX<<std::endl
      <<std::endl;

}


void ParseArgs(int argc, char**argv)
{

  int ch, errflg = 0;
  optarg = NULL;

  while ( !errflg && ((ch = getopt (argc, argv, "hCR:Q:r:q:M:D:O:e:l:a:z:d:i:s:t:p:c:x:L:")) != EOF) )
    switch (ch)
      {
      case 'h':
        PrintHelp (argv[0]);
        exit(1);
        break;

      case 'D':
        DELTAFILE=optarg;
        break;

      case 'R':
        RCE=optarg;
        break;

      case 'Q':
        QCE=optarg;
        break;

      case 'r':
        RNS=optarg;
        break;

      case 'q':
        QNS=optarg;
        break;

      case 'M':
        MUMCOORDS = optarg;
        break;

      case 'O':
        OUT = optarg;
        break;

      case 'e':
        MAX_EDGE_UNALN=atoi(optarg);
        break;

      case 'l':
        MIN_ALN_LEN=atoi(optarg);
        break;

      case 'a':
        MIN_BASES_ALN=atoi(optarg);
        break;

      case 'z':
        MIN_ZSTAT=atof(optarg);
        break;

      case 'd':
        MIN_DIFF_Z=atof(optarg);
        break;

      case 'i':
        MIN_MIS_OR=atoi(optarg);
        break;

      case 's':
        MIN_SUP=atoi(optarg);
        break;

      case 't':
        MIN_N_INS=atoi(optarg);
        break;

      case 'p':
        MIN_N_PER=atof(optarg);
        break;

      case 'C':
        PRINT_MUMS=false;
        break;

      case 'c':
        MIN_COV=atoi(optarg);
        if(MIN_COV<=0){
          std::cerr<<"c must be >0"<<endl;
          exit(1);
        }
        break;

      case 'L':
        MIN_LINK_COV=atoi(optarg);
        if(MIN_LINK_COV<0){
          std::cerr<<"L must be > 0"<<endl;
          exit(1);
        }
        break;

      case 'x':
        SCF_PREFIX=string(optarg);
        break;

      default:
        cout<<"Error parsing arguments"<<std::endl;
        errflg ++;
      }

  if (errflg > 0 || !RCE || !QCE || 
      !RNS || !QNS || !MUMCOORDS || !DELTAFILE)
  {
    PrintUsage (argv[0]);
    std::cerr << "Try '" << argv[0] << " -h' for more information.\n";
    exit (1);
  }

  if(!OUT){
    char x[3];
    OUT=x; 
    OUT[0]='o';
    OUT[1]='u';
    OUT[2]='t';
  }
}

void PrintHelp (const char * s){

  PrintUsage (s);
  std::cerr
    << "Pair-wise merge."
    << "\n\n"
    << "Required:" << std::endl
    << "-R -Q    Paths for .mateAn files. -R for primary genome. -Q for secondary genome\n"
    << "-r -q    Paths for .Ns files. -r for primary genome. -q for secondary genome\n"
    << "         The .Ns file should be a TAB delimited table with the following columns:\n"
    << "         1)chr_name, 2)start_pos, 3)end_pos. The start position should be 0-based\n"
    << "         while the end position should be 1-based. The ranges specify gap sequences\n"
    << "         present in the genome\n"
    << "-D       Path to the output of:\n"
    << "         delta-filter -1 in.maxmatch.delta\n"
    << "         Where \'in.maxmatch.delta\' is the delta file generated by aligning the two input\n"
    << "         assemblies using nucmer --maxmatch\n" 
    << "-M       Path to the output of:\n"
    << "         delta-filter -1 in.maxmatch.delta > show-coords -H -c -l -r -T\n\n"
    << "Optional:\n"
    << "\n  Output:\n"
    << "  -O       Output prefix. (Default: out )\n"
    << "  -C       Do not print .coords file ordered by primary assembly and"
    << "           secondary assembly sequence names\n"
    << "           (Default: true)\n"
    << "  -x       Prefix of scaffold names for the new metassembly. Names will be:\n"
    << "           prefix_0, prefix_1, and so on.\n"
    << "\n  CE statistic and coverage:\n"
    << "  -z       Threshold for the absolute value of the CE-statistic for considering local deviations from the expected\n" 
    << "           mean insert length as significant (Default: "<<MIN_ZSTAT<<")\n"
    << "  -d       Only make changes to the \"primary\" assembly if the difference in the CE-statistic values\n"
    << "           of the primary and the secondary assemblies is at least d. (Default: "<<MIN_DIFF_Z<<")\n"
    << "  -i       For inversion events, only consider those that have at least i misoriented mate pairs\n"
    << "           (Default: "<<MIN_MIS_OR<<")\n"
    << "  -s       Only consider those CE-statistic deviations that are supported by at least s mate pairs.\n"
    << "           (Default: "<<MIN_SUP<<")\n"
    << "  -c       Minimum coverage for considering a region as correct. This value is used along with the\n" 
    << "           CE-statistic value in order to find errors in the primary assembly (Default:"<<MIN_COV<<")\n"
    << "\n  Scaffolding:\n"
    << "  Whenever the edge of two primary scaffolds align contiguously within the same secondary scaffold\n"
    << "  there is evidence that two primary scaffolds can be linked into a single sequence:\n"
    << "  -e       When looking for alignments at the edge of primary scaffolds in order to find\n"
    << "           links, only consider the first non-spurious alignment in the range\n"
    << "           [edge,e]. Alignments that do not fall in this range will be ignored (Default: "<<MAX_EDGE_UNALN<<")\n"
    << "  -l       When linking, only use alignments with length>=l (Default: "<<MIN_ALN_LEN<<")\n"
    << "  -a       When linking, only use alignments with at least \'a\' perfectly matched bases\n"
    << "           (Default: "<<MIN_BASES_ALN<<")\n"
    << "  -L       When linking, only link scaffolds when the secondary scaffold has coverage higher than\n"
    << "           L (Default: "<<MIN_LINK_COV<<")\n"
    << "\n  Gap closure:\n"
    << "  Blocks of unaligned sequence in the primary assembly containing at least -t Ns in a row,\n"
    << "  and such that -p percent of the bases are Ns, are considered gap sequences. These gap\n"
    << "  sequences may be closed using the corresponding secondary sequence if it is not itself a\n" 
    << "  a gap sequence. By default, secondary sequences with at least one row of Ns equal or larger\n"
    << "  than 10 and such that 10\% of the total range are Ns are considered gap sequences. For primary\n"
    << "  sequences the following parameters apply:\n"
    << "  -t       Only close primary gap sequences with at least one row of Ns \n"
    << "           with length >= t (Default:"<<MIN_N_INS<<")\n"
    << "  -p       Only close primary gap sequences such that Ncount/length(block)>=p\n"
    << "           (Default:"<<MIN_N_PER<<")\n"
    << std::endl;
}

void PrintUsage (const char * s)
{
  std::cerr
    << "\nUSAGE: " << s << "  -R in.mateAn.primary -Q in.mateAn.secondary -r in.Ns.primary -q in.Ns.secondary -D in.1delta -M in.nuc.1coords -O out.prefix \n\n";
  return;
}
