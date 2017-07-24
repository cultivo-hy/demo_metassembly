///////////////////////////////////////////////////////////////////////////
// 
// author: Alejandro Hernandez Wences
// date: 21 Feb 2013
// 
// Class declarations for computing CEstatistic from a sam file
//
// see CEstat.cc
///////////////////////////////////////////////////////////////////////////


#ifndef __CESTAT_HH
#define __CESTAT_HH

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <map>
#include <vector>

// Variables
const int MAX_HEADER_SIZE=50;

// Classes

struct BDPfile;
struct SAMfile;
struct BDProw;

struct Insert
{

  std::string chr;
  int start,end,size;
  char strand1,strand2;
  bool propOrient;
  BDPfile * bdpfile;
  
  Insert(){};
  Insert(BDProw& bdp ); 
  bool filterOut (double mu0, double sd0, double zstat_threshold) const;
  bool zpos (double mu0, double sd0, double zstat_threshold) const; // true if zvalue > zstat_threshold
  bool zneg (double mu0, double sd0, double zstat_threshold) const; // true if zvalue < zstat_threshold 
  double zvalue (double mu0, double sd0) const;
  std::ostream& print (double mu0, double sd0, std::ostream& out) const;

};

struct BDProw
{

  std::string chr1,chr2,qname;
  int mapq;
  int start1,start2, end1, end2;
  char strand1,strand2;
  BDPfile * bdpfile;
  SAMfile * samfile;
 
  BDProw(){}; 
  BDProw(SAMfile& sam);
  BDProw(BDPfile& bdp);
  bool filterOut (int minMapQ) const;
  std::ostream& print (std::ostream& out);

};

struct BDPfile
{

  bool printFiltered;
  std::string bdp_file;
  std::ifstream bdpin;
  std::vector< BDProw > mpairs;
  int filtFz,filtMapq, filtJmpChrom, totalReads;
  
  BDPfile(){};
  BDPfile(std::string file);
  BDPfile(char * file);
  std::ostream& print (std::ostream& out) ;
  void nextBDP (BDProw& bdp);
  void loadBDPs();
  void sort();
  bool check_sort () ;
  std::ostream& printReadReport(std::ostream& out);
};

struct SAMrow
{

  std::string qname,rname, seq, cigar, mrname, rqual, opt;
  int flag,mapq, inslen;
  int start, mstart;

};

struct SAMfile
{

  bool printFiltered;
  std::string sam_file;
  int filtSingleReads, filtMapq, filtJmpChrom, totalReads; 
  std::ifstream samin;

  SAMfile(){};
  SAMfile(std::string file);
  SAMfile(char * file);
  SAMrow nextSAM();
  BDProw nextBDP();
  std::ostream& writeBDP(std::ostream& out, int minMapQ);
  std::ostream& writeBDPsort(std::ostream& out, int minMapQ);
  bool loadSeqLens(std::map<std::string, int >& ctglens);
  std::ostream& printReadReport(std::ostream& out);

};

struct mateAnRow
{

  int N,zneg,zpos,goodOrientCount;
  int sum,start,end;

  mateAnRow(){};
  int getN();
  int getZneg();
  int getZpos();
  int getZneut();
  int getGoodOrient();
  int getBadOrient();
  int getSum();
  int getStart();
  int getEnd();
  double zvalue(double mu0, double sd0);
  double newZstat(double mu0, double sd0, int delLength, int insLength);
  double newDeltaMean(double mu0, int delLength, int insLength);
  std::ostream& print(std::ostream& out, std::string chr, double mu0, double sd0);
  
};

struct mateAn
{

  std::map< std::string, std::vector<mateAnRow> > mateAns;

  //Parameters
  int MIN_N_STAT; // Only use regions with coverage >= MIN_N_STAT when approximating
                   // zstat distribution to standard normal
  int MIN_MIS_OR;
  int ZSTAT_THRE;
  int MIN_SUP;
  int MIN_COV;

  //Stats only considering mateAnRows with N>=MIN_N_STAT
  int totalRanges;
  int totalBases;
  double sumZstat;
  double sumZstatSquare;
  double minZstat;
  double maxZstat;
  double sumCov;
  double sumCovSquare;
  int minCov;
  int maxCov;
  
  //null hypothesis mu and sd
  double mu0,sd0;  

  void addMateAnRow(std::string chr, mateAnRow& mateRow);
  double meanZstat();
  double sdZstat();
  double meanCov();
  double sdCov();
  double suggestMu0();
  bool isRangeGoodCov(std::string chr, int start, int end, int minCov);
  std::ostream& print(std::ostream& out);  
  int search(std::string chr, int pos);
  char zsignal(std::string scf, int pos, char type);
  double zvalue(std::string scf, int pos);
  void updateStats(double mu0, double sd0);
  void updateStats();
  std::ostream& printStats(std::ostream& out);
  void readMateAn(std::ifstream& in, std::string prefix);  
  void clear();

  std::vector<mateAnRow>& operator[](std::string i){return mateAns[i];}

  mateAn(){};
  mateAn(std::ifstream& in, std::string prefix);
  mateAn(std::ifstream& in);
  mateAn(std::string file, std::string prefix);
  mateAn(double mu, double sd);
};

#endif // #ifndef __CESTAT__HH
