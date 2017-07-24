///////////////////////////////////////////////////////////////////////////////
////// 
////// author: Alejandro Hernandez Wences
////// date: 21 Feb 2013
////// 
////// Compute CE-statistic from mate pair alignments
//////
////// see CEstat.hh
/////////////////////////////////////////////////////////////////////////////////

#include "CEstat.hh"
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <list>
#include <stdlib.h>
#include <getopt.h>

using std::string;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::vector;
using std::map;
using std::ostream;
using std::min;
using std::max;
using std::list;

//-------------------------------------- Globals
//-------------------Options
int OPT_edgeDist=0;
int OPT_minCov=1;
int OPT_MIN_MAPQ=20;
bool OPT_SORT_BDP=true;
bool OPT_SAM=true;
bool OPT_CTGLEN=false;
bool OPT_FIT_MEAN_ZSTAT=true;
const char* OPT_MATES_FILE;
const char* OPT_CTGLEN_FILE;
const char* OPT_OUT_PREFIX;
double OPT_mu0=0;
double OPT_sd0=0;
double OPT_A=0;
double OPT_B=0;
double OPT_MINZ_FILTER=6;
double OPT_ZSTAT_CUTOFF=3;
double OPT_FIT_MU=0.05;
int OPT_MIN_N_STAT=30;

//-------------------Variables

//-------------------------------------- Function Defs
bool compareInsertEnds(Insert A, Insert B);

//-------------------------------------- Classes

struct RangeInsertSummary
{
  //current range info
  std::string chr;
  std::list<Insert> endStack;
  mateAnRow mateRow;
  
  //Parameters
  double mu0;
  double sd0;
  double zstat_threshold;
  
  mateAn * mateInf;
  std::map< std::string, int >ctglens;

  RangeInsertSummary();
  void clear();
  void conclude(int newStart);
  void loadInsert(Insert& insInf);
  void addMateAnRow(int endRange);

};

void RangeInsertSummary::clear()
{

  endStack.clear();
  mateRow.N=0;
  mateRow.goodOrientCount=0;
  mateRow.start=1;
  mateRow.sum=0;
  mateRow.end=1;
  mateRow.zneg=0;
  mateRow.zpos=0;
  chr.clear();

}

void RangeInsertSummary::addMateAnRow(int endRange)
{
  
  if(mateRow.start < OPT_edgeDist || endRange > ctglens[chr] - OPT_edgeDist)
    return;

  if(mateRow.N < OPT_minCov)
    return;

  mateRow.end=endRange;  

  mateInf -> addMateAnRow(chr, mateRow);

}

void RangeInsertSummary::conclude(int newStart)
{

  std::list<Insert>::reverse_iterator leftend;
  leftend=endStack.rbegin();
  int leftendpos;

  while(endStack.size() > 0 && leftend -> end < newStart){
    addMateAnRow(leftend -> end);
    mateRow.start=leftend -> end + 1;
    leftendpos= leftend -> end;
    while(endStack.size() > 0 && leftend -> end == leftendpos){
      //Update goodOrient and badOrient
      if(leftend -> propOrient)
        mateRow.goodOrientCount--;

      //Update zneg or zpos
      if(leftend -> zneg(mu0,sd0,zstat_threshold))
        mateRow.zneg--;
      else if(leftend -> zpos(mu0,sd0,zstat_threshold))
        mateRow.zpos--;
      
      mateRow.sum -= leftend -> size;
      mateRow.N--;
      endStack.pop_back();
      leftend=endStack.rbegin();
    }
  }
  if(mateRow.start < newStart)
    addMateAnRow(newStart -1);
}

void RangeInsertSummary::loadInsert(Insert& insInf)
{

  mateRow.start=insInf.start;
  
  if(insInf.propOrient)
    mateRow.goodOrientCount++;

  if(insInf.zneg(mu0,sd0,zstat_threshold))
    mateRow.zneg++;
  else if(insInf.zpos(mu0,sd0,zstat_threshold))
    mateRow.zpos++;

  mateRow.sum += insInf.size;
  mateRow.N++;
  endStack.push_front(insInf);
  endStack.sort(compareInsertEnds);
//  sort(endStack.begin(), endStack.end(), compareInsertEnds);
}

RangeInsertSummary::RangeInsertSummary ()
{

  clear();
  ctglens.clear();
  mateInf = NULL;

}

//-------------------------------------- Functions
bool compareInsertEnds(Insert A, Insert B)
{

  if(A.chr > B.chr)
    return(true);
  else if(A.chr == B.chr)
    if(A.end > B.end)
      return(true);

  return(false);

}

bool loadSeqLens(ifstream& CTG, map<string, int>& ctglens)
{

  string row,seqName;
  int seqLength;
  int found;
  while(CTG.good() && CTG.peek() != EOF){
    CTG >> seqName;
    getline(CTG,row);
    found=row.find_last_of("\t");
    if(found != string::npos){
      seqLength=atol(row.substr(found+1).c_str());
      ctglens[seqName]=seqLength;
    }else{
      return(false);
    }
    if(CTG.bad())
      return(false);
  }
  return(true);
}

void computeMu0Sd0(double A, double B , string in_bedpe, RangeInsertSummary& currMates)
{

  BDPfile BDP (in_bedpe);
  BDProw bdpr;

  if(A == B){
    cerr << "ERROR: A and B cannot be equal" << std::endl;
    exit(EXIT_FAILURE);
  }
  //Infer mu0 and sd0
  double mu0=0;
  double sd0=0;
  double N0=0;
  while(BDP.bdpin.good() && BDP.bdpin.peek() != EOF){
    BDP.nextBDP(bdpr);
    if(! bdpr.filterOut(OPT_MIN_MAPQ)){
      Insert ins=bdpr;
      if( ins.size >= A && ins.size <= B ){
        mu0 += ins.size;
        sd0 += pow(ins.size,2);
        N0++;
      } 
    }
  }
  currMates.mu0 = mu0/N0;
  currMates.sd0 = sd0/N0 - pow(currMates.mu0,2);
  currMates.sd0 = sqrt(currMates.sd0);

  BDP.bdpin.close();

}

void computeSd0(double mu0, double A, double B, string in_bedpe, RangeInsertSummary& currMates)
{

  BDPfile BDP (in_bedpe);
  BDProw bdpr;
  if(A == B){
    cerr << "ERROR: A and B cannot be equal" << std::endl;
    exit(EXIT_FAILURE);
  }
  //Infer mu0 and sd0
  double sd0=0;
  double N0=0;
  while(BDP.bdpin.good() && BDP.bdpin.peek() != EOF){
    BDP.nextBDP(bdpr);
    if(! bdpr.filterOut(OPT_MIN_MAPQ)){
      Insert ins=bdpr;
      if( ins.size >= A && ins.size <= B ){
        sd0 += pow(ins.size,2);
        N0++;
      } 
    }
  }
  currMates.sd0 = sd0/N0 - pow(currMates.mu0,2);
  currMates.sd0 = sqrt(currMates.sd0);

}

void computeCEstats(BDPfile& BDP, RangeInsertSummary& currMates, mateAn& mates)
{

  currMates.clear();
  BDProw bdpr;

  while(BDP.bdpin.good() && BDP.bdpin.peek() != EOF){
    BDP.nextBDP(bdpr);
    if(bdpr.filterOut(OPT_MIN_MAPQ))
      continue;

    Insert ins(bdpr);
    if(ins.filterOut(mates.mu0,mates.sd0,OPT_MINZ_FILTER))
      continue;

    if(currMates.chr.empty())
      currMates.chr=ins.chr;

    if(currMates.chr == ins.chr){
      if(ins.start < currMates.mateRow.start){
        cerr << "\n\n\nThe bedpe file must be ordered by chromosome and insert start position:\n"
             << currMates.mateRow.start << "\t>\t" << ins.start<<"\n\n";
        bdpr.print(cerr);
        exit(EXIT_FAILURE);
      }
      currMates.conclude(ins.start);
      currMates.loadInsert(ins);
    } else {
      currMates.conclude(currMates.ctglens[currMates.chr]+1);
      currMates.clear();
      currMates.chr=ins.chr;
      currMates.conclude(ins.start);
      currMates.loadInsert(ins);
    }
  }
  currMates.conclude(currMates.ctglens[currMates.chr]+1);

}

void ParseArgs(int argc, char** argv);
void PrintHelp(const char * s);
void PrintUsage(const char * s);
void printSetParameters(ostream& OUT);

//-------------------------------------- Main
int main (int argc, char ** argv)
{

  ParseArgs(argc, argv);

  RangeInsertSummary currMates;
  string in_bedpe;

  ofstream REPORT (string(OPT_OUT_PREFIX).append(".report").c_str());
  printSetParameters(REPORT);
  printSetParameters(cout);

  if(OPT_SAM){

    SAMfile SAM (OPT_MATES_FILE);

    if(!OPT_CTGLEN){
      if(! SAM.loadSeqLens(currMates.ctglens)){ //sam file has no headers and ctglens file was not provided
        std::cerr << "\n\nError:\n"
                  << "Your sam file contains no headers and you didn't\n"
                  << "provide any file with contig sizes."
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      ifstream CTG (OPT_CTGLEN_FILE);
      if(!loadSeqLens(CTG,currMates.ctglens)){
        std::cerr << "ERROR:\nBe sure that your contig size file has a TAB delimited\n"
                  << "table with the following columns:\n"
                  << "contig_name contig_size"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    
    string out_bedpe=string(OPT_OUT_PREFIX).append(".sort.bedpe");
    ofstream BDP_OUT(out_bedpe.c_str());
    if(OPT_SORT_BDP){
      cout << "---------- Writing Sort bedpe\n"
           << "----------------------"<<std::endl<<std::endl;
      REPORT << "---------- Writing Sort bedpe\n"
           << "----------------------"<<std::endl<<std::endl;
      SAM.writeBDPsort(BDP_OUT, OPT_MIN_MAPQ);
      SAM.printReadReport(cout);
      SAM.printReadReport(REPORT);
      in_bedpe=out_bedpe;
    }else{
      cout << "---------- Writing bedpe\n"
           << "----------------------"<<std::endl<<std::endl;
      REPORT << "---------- Writing bedpe\n"
           << "----------------------"<<std::endl<<std::endl;
      SAM.writeBDP(BDP_OUT, OPT_MIN_MAPQ);
      SAM.printReadReport(cout);
      SAM.printReadReport(REPORT);
      exit(EXIT_SUCCESS);
    }
  }

  if(currMates.ctglens.empty())
    if(!OPT_CTGLEN){
      cerr << "ERROR:\nYou must provide a file with the following columns:\n"
           << "contig_name contig_size\n"
           << std::endl;
      exit(EXIT_FAILURE);
    } else {
      ifstream CTG (OPT_CTGLEN_FILE);
      if(! loadSeqLens(CTG,currMates.ctglens)){
        std::cerr << "ERROR:\nBe sure that your contig size file has a TAB delimited\n"
                  << "table with the following columns:\n"
                  << "contig_name contig_size"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }

  if(in_bedpe.empty())
    in_bedpe=string(OPT_MATES_FILE);

  if(OPT_A > 0 && OPT_B > 0){
    cout << "---------- Infering mu0 and sd0 with A: "<<min(OPT_A,OPT_B)
                                          <<" and B: "<<max(OPT_A,OPT_B)<<"\n"
         << "----------------------"<<std::endl<<std::endl;
    REPORT << "---------- Infering mu0 and sd0 with A: "<<min(OPT_A,OPT_B)
                                            <<" and B: "<<max(OPT_A,OPT_B)<<"\n"
           << "----------------------"<<std::endl<<std::endl;
    computeMu0Sd0(min(OPT_A,OPT_B),max(OPT_A,OPT_B),in_bedpe,currMates);
  } else if (OPT_mu0 >= 0 && OPT_sd0 >= 0 ){
    currMates.mu0 = OPT_mu0;
    currMates.sd0 = OPT_sd0;
  } else {
    cerr<< "ERROR: You must provide either A and B, or mu0 and sd0\n";
    exit(EXIT_FAILURE);
  }

  cout<< "---------- Computing CE-stat:" <<std::endl
      << "----------------------\n\n" <<std::endl;
  REPORT<< "---------- Computing CE-stat:" <<std::endl
        << "----------------------\n\n" <<std::endl;

  BDPfile BDP (in_bedpe);
  mateAn mates (currMates.mu0,currMates.sd0);
  mates.MIN_N_STAT=OPT_MIN_N_STAT;
  currMates.mateInf = &mates;
  currMates.zstat_threshold = OPT_ZSTAT_CUTOFF;

  computeCEstats(BDP, currMates, mates);

  BDP.printReadReport(REPORT);
  BDP.printReadReport(cout);
  mates.printStats(cout);
  mates.printStats(REPORT);

  if(OPT_FIT_MEAN_ZSTAT){
    bool corrected=false;
    cout<< "---------- Fitting Mu0:" <<std::endl
      << "----------------------\n\n" <<std::endl;
  REPORT<< "---------- Fitting Mu0:" <<std::endl
        << "----------------------\n\n" <<std::endl;
    while(fabs(mates.meanZstat()) > OPT_FIT_MU){
      mates.mu0=mates.suggestMu0();
      mates.updateStats();
    mates.printStats(cout);
    mates.printStats(REPORT);
      corrected=true;
    }

    cout<< "---------- Computing CE-stat:" <<std::endl
      << "----------------------\n\n" <<std::endl;
  REPORT<< "---------- Computing CE-stat:" <<std::endl
        << "----------------------\n\n" <<std::endl;

    if(corrected){
      BDPfile BDP(in_bedpe); 
      currMates.mu0=mates.mu0;
      mates.clear();
      computeCEstats(BDP,currMates, mates);
      BDP.printReadReport(REPORT);
      BDP.printReadReport(cout);
      mates.printStats(cout);
      mates.printStats(REPORT);
    }
  }

  cout<< "---------- Writing to file:"<<string(OPT_OUT_PREFIX).append(".mateAn") <<std::endl
      << "----------------------\n\n" <<std::endl;
  REPORT<< "---------- Writing to file:"<<string(OPT_OUT_PREFIX).append(".mateAn") <<std::endl
        << "----------------------\n\n" <<std::endl;

  ofstream MATE (string(OPT_OUT_PREFIX).append(".mateAn").c_str());
  mates.print(MATE);

}

void PrintUsage(const char * s)
{

  cerr << "\nUSAGE: " << s << " [options] -p <prefix> (-m <float> -s <float> || -A <float> -B <float>) (-s <file> [-l <file>] || -b <file> -l <file>)\n\n";

}

void PrintHelp(const char * s)
{

  PrintUsage(s);

  cerr << "Print the compression-expansion (CE) statistic value at the beginning and end\n"
       << "of each insert across each contig. The CE-statistic quantifies how compressed\n"
       << "or expanded are the set of mate pairs spanning a given position in comparison\n" 
       << "to the expected insert size. Formally, it computes a z-test to detect differences\n"
       << "between the local mean insert length and the global (expected) mean insert length"
       << "\n\n"
       << "-------Required: \n\n"
       << "--Mapped mate-pair information:\n"
       << "Mapped mate-pairs can be provided either through a SAM file:\n" 
       << "  -a           SAM file\n"
       << "or a bedpe file and a sequence-lengths file:\n"
       << "  -b           -b BEDPE file with one range for each mapped insert. This file can be generated\n"
       << "               by mateAn with the -u option\n"
       << "  -l           file with scf_name/scf_length TAB delimited table\n"
       << "If the SAM file lacks \'@SQ  SN:scf_name  LN:scf_length\' headers the -l option must be specified\n"
       << "If a SAM file is provided, a sorted bedpe file is written by default unless you specify\n"
       << "the -u option, in which case you will have to sort the resulting bedpe file by chromosome name and\n"
       << "insert start position using your own tools and re-run mateAn using the -b and -l options\n"
       << "(not recommended).\n\n"
       << "--Null hypothesis mean insert length (mu0) and standard deviation (sd0):\n"
       << "mu0 and sd0 can either be inferred from the mapped mate-pairs:\n"
       << "  -A -B        mu0 and sd0 will be inferred using inserts with length within the range [A,B]\n"
       << "               Setting appropriate values for A and B is recommended to reduce the effect\n"
       << "               of misaligned mate-pairs\n"
       << "or can be provided directly:\n"
       << "  -m   Null hypothesis mean insert length (mu0)\n"
       << "  -s   Null hypothesis insert length standard deviation (sd0)\n"
       << "If -m and -s are provided then -A and -B will be ignored. Otherwise, the null hypothesis mean insert length\n"
       << "and standard deviation will be computed to be used for the CE-statistic computation."
       << "\n\n"
       << "  -p   output prefix\n\n"
       << "-------Optional:\n\n"
       << "  -h   print this message\n"
       << "  -e   only print the CE-statistic for chromosome positions that are at least <int> bases\n"
       << "       apart from the edges (Default: "<<OPT_edgeDist<<")\n"
       << "  -z   zstat threshold for considering deviations from the mean (mu0) as significant.(Default: "<<OPT_ZSTAT_CUTOFF<<" sd0s)\n"    
       << "  -q   minimum mapping quality (Only if a SAM file is provided).(Default: "<<OPT_MIN_MAPQ<<")\n"
       << "  -N   only print the CE-statistic for positions with coverage >= <int>. (Default: "<<OPT_minCov<<")\n"
       << "  -f   filter out inserts whose lengths deviate by at least \'f\' sd0 standard deviations from the expected\n"
       << "       mean insert length (mu0). These inserts won't be used in the CE-statistic computation. Setting an \n"
       << "       appropriate value for -f is recommended to reduce the effect of misaligned mate-pairs on the CE-statistic\n"
       << "       (Default: "<<OPT_MINZ_FILTER<<")\n"
       << "  -c   adjust mu0 and sd0 (and recompute the CE-statistic) until |mean CE-stat| <= <double>. (Default: "<<OPT_FIT_MU<<")\n"
       << "  -o   do not adjust mu0 and sd0, compute the CE-statistic only once.\n"
       << "  -n   when computing mu0 and sd0, only consider positions with coverage >= <int>\n"
       << "       (Default: "<<OPT_MIN_N_STAT<<")\n";
} 

void ParseArgs (int argc, char ** argv)
{

  int ch, errflg =0;
  optarg = NULL;

  if(argc <= 1){
    PrintHelp(argv[0]);
    exit(EXIT_FAILURE);
  } 

  while(!errflg && (( ch = getopt (argc, argv, "hm:a:e:z:q:s:b:l:N:p:f:n:A:B:c:ou")) != EOF) )
    switch(ch){
      
      case 'h':
        PrintHelp(argv[0]);
        exit(EXIT_FAILURE);
        break;

      case 'm':
        OPT_mu0=atof(optarg);
        break;

      case 's':
        OPT_sd0=atof(optarg);
        break;

      case 'e':
        OPT_edgeDist=atoi(optarg);
        break;

      case 'z':
        OPT_ZSTAT_CUTOFF=atof(optarg);
        break;

      case 'q':
        OPT_MIN_MAPQ=atoi(optarg);
        break;

      case 'a':
        OPT_SAM=true;
        OPT_MATES_FILE=optarg;
        break;

      case 'b':
        OPT_SAM=false;
        OPT_MATES_FILE=optarg;
        break;

      case 'l':
        OPT_CTGLEN=true;
        OPT_CTGLEN_FILE=optarg;
        break;

      case 'N':
        OPT_minCov=atoi(optarg);
        break;

      case 'p':
        OPT_OUT_PREFIX=optarg;
        break;

      case 'f':
        OPT_MINZ_FILTER=atof(optarg);
        break;

      case 'n':
        OPT_MIN_N_STAT=atoi(optarg);
        break;

      case 'A':
        OPT_A=atof(optarg);
        break;

      case 'B':
        OPT_B=atof(optarg);
        break;

      case 'c':
        OPT_FIT_MEAN_ZSTAT=true;
        OPT_FIT_MU=atof(optarg);
        break;

      case 'o':
        OPT_FIT_MEAN_ZSTAT=false;
        break;

      case 'u':
        OPT_SORT_BDP=false;
        break;

      default:
        cerr << "ERROR" << std::endl;
        errflg++;
    }
  
  if(errflg > 0                     ||
     OPT_SORT_BDP && (string(OPT_MATES_FILE).empty() || 
     string(OPT_OUT_PREFIX).empty() || 
     !( (OPT_A && OPT_B) || (OPT_mu0 && OPT_sd0) ) )){
    PrintUsage(argv[0]);
    cerr << "Try " << argv[0] << " -h for more information\n";
    exit(EXIT_FAILURE);
  }

  if(OPT_SORT_BDP && !OPT_SAM && !OPT_CTGLEN){
    PrintUsage(argv[0]);
    cerr << "You must provide a file with sequence lengths in a\n"
         << "TAB delimited table: seq_name length\n\n";
    cerr << "Try " << argv[0] << " -h for more information\n";
    exit(EXIT_FAILURE);
  }

  if(OPT_SORT_BDP && OPT_A && OPT_B && OPT_mu0 && OPT_sd0 ){
    cerr <<"ERROR:\n" 
        <<"You should provide either A and B, or mu0 and sd0\n";
    exit(EXIT_FAILURE);
  }
}

void printSetParameters(ostream& out)
{
  
  out << "---------- Parameters: \n"
      << "----------------------\n";
  if(OPT_SAM){
    out << "-a : " << OPT_MATES_FILE << "\n";
    if(OPT_CTGLEN)
      out << "-l : " << OPT_CTGLEN_FILE << "\n";
  }else{
    out << "-b : " << OPT_MATES_FILE << "\n"
        << "-l : " << OPT_CTGLEN_FILE << "\n";
  }
  if(OPT_A && OPT_B){
    out << "-A : " << OPT_A << "\n"
        << "-B : " << OPT_B << "\n";
  }
  if(OPT_mu0 && OPT_sd0){
    out << "-m : " << OPT_mu0 << "\n"
        << "-s : " << OPT_sd0 << "\n";
  }
  
  out << "-p : " << OPT_OUT_PREFIX << "\n"
      << "-e : " << OPT_edgeDist << "\n"
      << "-z : " << OPT_ZSTAT_CUTOFF << "\n"
      << "-q : " << OPT_MIN_MAPQ << "\n"
      << "-N : " << OPT_minCov << "\n"
      << "-f : " << OPT_MINZ_FILTER << "\n";
  
  if(OPT_FIT_MEAN_ZSTAT)
    out << "-c : " << OPT_FIT_MU << "\n";
  else
    out << "-o : (Do not fit CE-stat mean and sd)\n";

  out << "-n : " << OPT_MIN_N_STAT << "\n";

}             

