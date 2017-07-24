/////////////////////////////////////////////////////////////////////////
// 
// author: Alejandro Hernandez Wences
// date: 21 Feb 2013
// 
// Class declarations for merging assemblies
//
// see merge.cc
/////////////////////////////////////////////////////////////////////////

#ifndef __MERGE_ASSEM_HH
#define __MERGE_ASSEM_HH

#include "delta.hh"
#include "tigrinc.hh"
#include "CEstat.hh"
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <map>
#include <vector>

struct MetaRow;
struct Metassembly;

struct AlignStats
  //-- DeltaFile alignment statistics data auxiliary structure
{
  int sQ, eQ, sR, eR;
  std::vector<int> Delta;
};

struct DeltaFile
{
  std::string DeltaName;

  int query2ref(std::string qscf, int qpos, std::string rscf, int rposA, int rposB); 
  std::vector<int> query2ref(std::string qscf, int qpos, std::string rscf);
  std::map< std::string, std::vector<int> >query2ref(std::string qscf, int qpos);
  int ref2query(std::string rscf, int rpos, std::string qscf, int qposA, int qposB);
  std::vector<int> ref2query(std::string rscf, int rpos, std::string qscf);
  std::map< std::string, std::vector<int> > ref2query(std::string rscf, int rpos);

  DeltaFile(std::string inDeltaName): DeltaName(inDeltaName){};
  DeltaFile(){};
  private:

  std::vector<int> transCoords(std::vector<AlignStats> Aligns, 
                               std::string Id1, int pos1, 
                               std::string Id2, 
                               bool isReferenceCoord);
  int transCoords(std::vector<AlignStats> Aligns, 
                               std::string Id1, int pos1, 
                               std::string Id2, int posBs, int posBe, 
                               bool isReferenceCoord);
  void parseDeltaR(std::map< std::string, std::vector<AlignStats> > & Aligns, std::string IdR);
  void parseDeltaQ(std::map< std::string, std::vector<AlignStats> > & Aligns, std::string IdR);
  
};

// mumRow -> CoordsRow
struct CoordsRow
{
  int start1,end1,start2,end2,alnsize1,alnsize2;
  double ident,cover1,cover2;
  int scfLen1,scfLen2;
  std::string scf1,scf2;

  std::ostream& print(std::ostream& out);

  CoordsRow(std::ifstream& MUM);
};

struct ordCoords
{
// Auxiliary structure to sort .coords file by reference/query coords.
// See MUMfile::loadFile
  std::string scf;
  int start,end;
  int line;
};

bool compareOrdCoords(ordCoords A, ordCoords B);

struct MUMfile
{
  std::vector<CoordsRow> coords;
  std::map<std::string,std::vector<int> > queryScfs;
  std::map<std::string,std::vector<int> > refScfs;
  std::ifstream coordsin;
  int minOverAlnLength;
  int minBasesAligned;
  int maxEdgeNotAligned;
  DeltaFile delta; 
 
  std::ifstream& loadFile();
  void printMums(std::string out_pre);
  std::vector<int>::iterator findQueEdgeAln(std::string scf, std::string edge);
  std::vector<int>::iterator findRefEdgeAln(std::string scf, std::string edge);
  std::vector<int>::iterator adjacentAlnQue(std::vector<int>&alns, std::vector<int>::iterator pos, char direction);
  std::vector<int>::iterator adjacentAlnRef(std::vector<int>&alns, std::vector<int>::iterator pos, char direction);
  bool QscfBefore(std::string scf, std::vector<int>& alns, std::vector<int>::iterator last);
  std::vector<int>::iterator QscfAfter(std::string scf, std::vector<int>& alns, std::vector<int>::iterator begin);
  bool analyseInterInsert(std::vector<int>&alns, std::vector<int>::iterator A, std::vector<int>::iterator B, mateAn& ce);
  MUMfile(std::string in_coords, std::string inDelta, int iminOver, int iminBas, int imaxEdge);
  MUMfile(){};
};

struct Nsfile
{

  std::map< std::string, std::vector< std::vector<int> > > Ndata;
  std::string ns_file;
  std::ifstream nsin;

  std::ifstream& loadFile(std::ifstream& in, std::string prefix);
  bool isNinsert(std::string scf, int pos1, int pos2, int minNinsert, double minNinsPer);
  int search(std::string scf, int pos1, int pos2);
  std::vector< std::vector<int> >& operator[](std::string i){return Ndata[i];}

  Nsfile(std::string filename, std::string prefix);
  Nsfile(){};
};

struct MetaLink
{

  std::string refScf1,refScf2,queScf,edgeScf1,edgeScf2;
  int refStart1, refEnd1, refStart2, refEnd2, queStart1, queEnd1, queStart2, queEnd2, aln1,aln2;
  char queOrient1, queOrient2;

  std::ostream& print(std::ostream& out);
};

struct metaLNcoords
{

  MetaLink mln;
  char orientation;
  std::string edge; //either "start" or "end"
  std::string desc;
  int delta;

  void print(std::string metaScf, //metassembly scf
             int mstart, //metassembly start
             int mend,   //metassembly end
             MetaRow* mrow,
             Metassembly * meta,
             std::ofstream& out
             );
  metaLNcoords(MetaLink imln, char ior, std::string iedge, int idelta): desc("")
  {
    mln=imln;
    orientation=ior;
    edge=iedge;
    delta=idelta;
  }
};

struct metaNAcoords
{
  std::string rscf,qscf;
  int rnstart,rnend,qostart,qoend;
  std::string qoside; //left or righti
  char qor;
  void print(std::string metaScf, int mstart, int mend, 
             MetaRow* mrow, Metassembly* meta, std::ofstream& out);
  metaNAcoords(std::string irscf, int irnstart, int irnend, //reference N row
               std::string iqscf, int iqostart, int iqoend,       //query overhang
               std::string iqoside, char iqor
              )
  {
    rscf=irscf;
    qscf=iqscf;
    rnstart=irnstart;
    rnend=irnend;
    qostart=iqostart;
    qoend=iqoend;
    qoside=iqoside;
    qor=iqor;
  }
};

struct metaNBcoords
{
  std::string rscf,qscf;
  int qostart,qoend,rnstart,rnend,ristart,riend,rbkp;
  std::string qoedge;
  char qor;
  void print(std::string metaScf, int mstart, int mend,
             MetaRow* mrow, Metassembly* meta, std::ofstream& out);
  metaNBcoords(std::string iqscf,int iqostart, int iqoend, //query overhang
               std::string irscf,int irnstart, int irnend, //reference N row 
               int iristart, int iriend, //reference unmapped insert
               int irbkp,                  //reference breakpoint
               std::string iqoedge, char iqor
              )
  {
    qscf=iqscf;
    rscf=irscf;
    qostart=iqostart;
    qoend=iqoend;
    rnstart=irnstart;
    rnend=irnend;
    ristart=iristart;
    riend=iriend;
    rbkp=irbkp;
    qoedge=iqoedge;
    qor=iqor;
  }
};

struct metaNCcoords
{
  std::string qscf1,qscf2,rscf;
  int ristart,riend,rlnstart,rlnend,rrnstart,rrnend;
  int qrostart,qroend,qlostart,qloend;
  char qror,qlor; 
  void print(std::string metaScf, int mstart, int mend,
             MetaRow* mrow, Metassembly* meta, std::ofstream& out);
  metaNCcoords(std::string iqscf1, int iqlostart, int iqloend, char iqlor, //left query scf overhang
               std::string irscf, int irlnstart, int irlnend, //left reference N row
               int iristart, int iriend,   //reference unmapped insert
               int irrnstart, int irrnend,  //right reference N row
               std::string iqscf2, int iqrostart, int iqroend, char iqror //right query scf overhang
              )
  {
    qscf1=iqscf1;
    qscf2=iqscf2;
    rscf=irscf;
    ristart=iristart;
    riend=iriend;
    rlnstart=irlnstart;
    rlnend=irlnend;
    rrnstart=irrnstart;
    rrnend=irrnend;
    qrostart=iqrostart;
    qroend=iqroend;
    qlostart=iqlostart;
    qloend=iqloend;
    qror=iqror;
    qlor=iqlor;
  }
};

struct metaN1coords
{
  std::string qscf,rscf;
  int qistart,qiend,ristart,riend,rnstart,rnend; 
  char qor;
  void print(std::string metaScf, int mstart, int mend,
             MetaRow* mrow, Metassembly* meta, std::ofstream& out);
  metaN1coords(std::string iqscf,int iqistart, int iqiend, char iqor,//query insert 
               std::string irscf,int iristart, int iriend, //reference insert
               int irnstart, int irnend  //reference N row
              )
  {
    qscf=iqscf;
    rscf=irscf;
    qistart=iqistart;
    qiend=iqiend;
    qor=iqor;
    ristart=iristart;
    riend=iriend;
    rnstart=irnstart;
    rnend=irnend;
  }
};

struct metaQIcoords
{
  std::string rscf,qscf;
  int qistart,qiend,ristart,riend;
  char qor;
  void print(std::string metaScf, int mstart, int mend,
             MetaRow* mrow, Metassembly* meta, std::ofstream& out);
  metaQIcoords(std::string iqscf,int iqistart, int iqiend, char iqor,//query insert
               std::string irscf,int iristart, int iriend  //reference insert
              )
  {
    qscf=iqscf;
    rscf=irscf;
    qistart=iqistart;
    qiend=iqiend;
    ristart=iristart; 
    riend=iriend;
  }
};

struct metaBKPcoords
{
  std::string rscf,qscf;
  int rbkp,qostart,qoend;
  char qor;
  std::string qoside;
  void print(std::string metaScf, int mstart, int mend,
             MetaRow* mrow, Metassembly* meta, std::ofstream& out);
  metaBKPcoords(std::string irscf, int irbkp, 
                std::string iqscf,int iqostart, int iqoend,
                std::string iqoside, char iqor
               )
  {
    rscf=irscf;
    qscf=iqscf;
    rbkp=irbkp;
    qostart=iqostart;
    qoend=iqoend;
    qor=iqor;
    qoside=iqoside;
  }
};

struct MetaRow
{

  std::string scf;
  int start, end, file, aln;
  char orientation, sigRef1, sigRef2, sigQue1, sigQue2;
  bool isInv;
  metaLNcoords* mlnc;
  metaNAcoords* mNAc;
  metaNBcoords* mNBc;
  metaNCcoords* mNCc;
  metaN1coords* mN1c;
  metaQIcoords* mQIc;
  metaBKPcoords* mBKPc;
  void changeOrientation();
};

struct name_len
{ 
//auxiliary structure to sort scaffold name-length pairs by length
  std::string name;
  int length;
};

struct Metassembly
{
  std::string outp;
  std::string mscfPrefix;
  std::map<std::string, std::vector<MetaRow> > metassem;
  std::map<std::string, int> metaLengths;
  std::map<std::string, std::string> Scf_Meta;
  std::map<std::string, std::map<int, std::map<int,char> > > inversions;
  std::map<std::string, MetaLink>scfsLinks;

  MUMfile mum;

  mateAn Qce;
  mateAn Rce;

  Nsfile QNs;
  Nsfile RNs;

  int ScfCount;
  int countInversions;
  int countNevents;
  int countNret;
  int countNA;
  int countNB;
  int countNC;
  int countLN;
  int countBKP;
  int countQI;
  int countNrows;

  std::vector<metaLNcoords*> LNs;
  std::vector<metaNAcoords*> NAs;
  std::vector<metaNBcoords*> NBs;
  std::vector<metaNCcoords*> NCs;
  std::vector<metaN1coords*> N1s;
  std::vector<metaQIcoords*> QIs;
  std::vector<metaBKPcoords*> BKPs;

  int minLinkCov;
  int minCov;
  int minNinsert;
  double minNinsPer;
  int QmaxNinsert; //Only close gaps with sequences with Ncount < QmaxNinsert
  double QmaxNinsPer;
  double minDiffZ;
  double minZstat;
  double minLenSDs;
  double maxQinsLength;
  double maxRinsLength;

  MetaRow* addMetaRow(std::string mscf,
                  std::string edge,
                  int file,
                  std::string seqScf,
                  int seqStart,
                  int seqEnd,
                  char orientation,
                  char sigRef1,
                  char sigRef2,
                  char sigQue1,
                  char sigQue2,
                  int aln);
  std::string findScfLink(std::vector<int>::iterator refAln, std::string edgeScf1);
  void joinLinks();
  void printStats();
  void printLinks();
  void printMetassem();
  std::vector<int>::iterator analyseAln(std::string refScf,
                                        std::vector<int>& rScfAlns,
                                        std::vector<int>::iterator curScfAln,
                                        std::string metaScf,
                                        bool leftLink,
                                        std::vector<int>::iterator leftAln,
                                        bool rightLink,
                                        std::vector<int>::iterator rightAln
                                       );
  Metassembly(std::string ioutp,
              std::string scf_prefix,
              char * inDelta,
              char * mumcoords,
              int iminOver,
              int iminBas,
              int imaxEdge,
              char * Rcefile,
              char * Qcefile,
              int minMisor,
              int minSup,
              char * RNsfile,
              char * QNsfile,
              int iminCov,
              double iminZstat,
              double iminDiffZ,
              int iminLinkCov,
              int iminNinsert,
              double iminNinsPer);
  Metassembly(){};

private:
  std::ofstream& printMetaRow(std::ofstream& OUTF, std::string metScf, int mstart, int mend,
                              MetaRow& mrow);
  char getRefOrientation(std::string scf, int start, int end);
  int lastPosAdded2Meta(std::string scf, std::string metaScf);  
  void addLink(int alnScf1, std::string edgeScf1, int alnScf2, std::string edgeScf2);
  std::string findEdgeInMeta(std::string metaScf, 
                             std::string refScf, 
                             int position, 
                             MetaLink ln, 
                             std::string prefEdge);  
  void joinMetassem(MetaLink& ln, std::string edgeMeta1, int startScf1copy, std::string edgeMeta2, int startScf2copy, bool addedQueInter);
  void addLNcoord(MetaRow* meta, MetaLink imln, std::string iedge, int idelta);
  void addNAcoord(MetaRow* meta, std::string irscf, int irnstart, int irnend, 
                                 std::string iqscf, int iqostart, int iqoend,
                                 std::string iqoside, char iqor);
  void addNBcoord(MetaRow* meta, std::string iqscf, int iqostart, int iqoend, 
                                 std::string irscf, int irnstart, int irnend, 
                                 int iristart, int iriend, 
                                 int irbkp,
                                 std::string iqoedge, char iqor);
  void addNCcoord(MetaRow* meta, std::string iqscf1, int iqlostart, int iqloend, char iqlor,
                                 std::string irscf, int irlnstart, int irlnend,
                                 int iristart, int iriend,   
                                 int irrnstart, int irrnend,  
                                 std::string iqscf2, int iqrostart, int iqroend, char iqror );
  void addN1coord(MetaRow* meta, std::string iqscf, int iqistart, int iqiend, char iqor,
                                 std::string irscf, int iristart, int iriend, 
                                 int irnstart, int irnend);
  void addQIcoord(MetaRow* meta, std::string iqscf, int iqistart, int iqiend, char iqor, 
                                 std::string irscf, int iristart, int iriend);
  void addBKPcoord(MetaRow* meta, std::string irscf, int irbkp, 
                                  std::string iqscf, int iqostart, int iqoend,
                                   std::string iqoside, char iqor);
};

#endif //#ifndef __MERGE_ASSEM_HH
