///////////////////////////////////////////////////////////////////////////////
//// 
//// author: Alejandro Hernandez Wences
//// date: 21 Feb 2013
//// 
//// Source for functions in CEstat.hh
////
//// see CEstat.hh
///////////////////////////////////////////////////////////////////////////////

#include "CEstat.hh"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>

bool compareBDProws(BDProw A, BDProw B)
{

  Insert Ai(A);
  Insert Bi(B);

  if(Ai.chr < Bi.chr)
    return(true);
  else if (Ai.chr == Bi.chr && Ai.start < Bi.start)
    return(true);
  else 
    return(false);

}

//===================================================== Insert ===============//
//-----------------------------------filterOut--------------------------------//
bool Insert::filterOut(double m0, double sd0, double zstat_threshold) const
{

  if(fabs(zvalue(m0,sd0)) > zstat_threshold ){
    if(bdpfile != NULL)
      bdpfile->filtFz++;
    return(true);
  }else
    return(false);

}

//-----------------------------------print------------------------------------//
std::ostream& Insert::print(double m0, double sd0, std::ostream& out) const
{

  out << chr << "\t"
      << start << "\t"
      << end << "\t"
      << size << "\t"
      << strand1 << "\t"
      << strand2 << "\t"
      << zvalue(m0,sd0) << std::endl;

  return(out);
}

//-----------------------------------zvalue-----------------------------------//
double Insert::zvalue(double mu0, double sd0) const
{

  return((double(size)-mu0)/sd0);

}

//-----------------------------------zpos-------------------------------------//
bool Insert::zpos(double mu0, double sd0, double zstat_threshold) const
{

  return(zvalue(mu0,sd0) > zstat_threshold);

}

//-----------------------------------zneg-------------------------------------//
bool Insert::zneg(double mu0, double sd0, double zstat_threshold) const
{
 
  return(zvalue(mu0,sd0) < -1*zstat_threshold);

}

//-----------------------------------Insert constructor-----------------------//
Insert::Insert(BDProw& bdp)
{

  if(bdp.chr1 != bdp.chr2){
    std::cerr<< "ERROR: Cannot instantiate Insert from bdp, bdp.chr1 != bdp.chr2"<<std::endl;
    exit(EXIT_FAILURE);
  }
  chr=bdp.chr1;
  start=std::min(bdp.end1,bdp.end2)+1; //one past the end of downstream read
  end=std::max(bdp.start1,bdp.start2); //start of upstream read
  size=end-start;
  if(size < 0){
    end=start;
    size=0;
  }
  strand1=bdp.strand1;
  strand2=bdp.strand2;
  propOrient = strand1 != strand2;
  if(bdp.bdpfile == NULL && bdp.samfile == NULL){
    std::cerr<< "ERROR: Cannot instantiate Insert from bdp, bdpfile && samfile are NULL"<<std::endl;
    exit(EXIT_FAILURE);
  }
  bdpfile=bdp.bdpfile;
 
}

//===================================================== BDProw ===============//
//----------------------------------- filterOut ------------------------------//
bool BDProw::filterOut (int minMapQ) const
{

  if(chr1 != chr2){
    if(bdpfile != NULL)
      bdpfile->filtJmpChrom++;
    else
      samfile->filtJmpChrom++;
    return(true);
  }  

  if(mapq < minMapQ){
    if(bdpfile != NULL)
      bdpfile->filtMapq++;
    else
      samfile->filtMapq++;
    return(true);
  }
  
  return(false);

}
 
//----------------------------------- print ----------------------------------//
std::ostream& BDProw::print (std::ostream& out)
{

  out << chr1 << "\t"
      << start1 - 1 << "\t" //bedpe start is 0-based
      << end1 << "\t" 
      << chr2 << "\t"
      << start2 - 1 << "\t" //bedpe start is 0-based
      << end2 << "\t"
      << qname << "\t"
      << mapq << "\t"
      << strand1 << "\t"
      << strand2 << std::endl;

  return(out);
}

//----------------------------------- BDProw(SAMfile) -------------------------//
BDProw::BDProw(SAMfile& sam){

  BDProw bdp=sam.nextBDP();
  chr1=bdp.chr1;
  chr2=bdp.chr2;
  qname=bdp.qname;
  mapq=bdp.mapq;
  start1=bdp.start1;
  start2=bdp.start2;
  end1=bdp.end1;
  end2=bdp.end2;
  strand1=bdp.strand1;
  strand2=bdp.strand2;
  
  bdpfile=bdp.bdpfile;
  samfile=bdp.samfile;  

}

//----------------------------------- BDProw(BDPfile) ------------------------//
BDProw::BDProw(BDPfile& bdpf)
{

  BDProw bdp;
  bdpf.nextBDP(bdp);
  chr1=bdp.chr1;
  chr2=bdp.chr2;
  qname=bdp.qname;
  mapq=bdp.mapq;
  start1=bdp.start1;
  start2=bdp.start2;
  end1=bdp.end1;
  end2=bdp.end2;
  strand1=bdp.strand1;
  strand2=bdp.strand2;

  bdpfile=bdp.bdpfile;
  samfile=bdp.samfile;

}

//===================================================== BDPfile ==============//
//----------------------------------- print ----------------------------------//
std::ostream& BDPfile::print(std::ostream& out)
{

  std::vector< BDProw >::iterator it;
  for(it=mpairs.begin(); it<mpairs.end(); it++){
    it->print(out);
    if(out.bad())
      return(out);
  }
  return(out);

}


//----------------------------------- nextBDP --------------------------------//
void BDPfile::nextBDP(BDProw& bdp)
{

  if(!bdpin.eof() && 
     !bdpin.fail() && 
     !bdpin.bad() &&
     bdpin.peek() != EOF){
    bdpin >> bdp.chr1
          >> bdp.start1
          >> bdp.end1
          >> bdp.chr2
          >> bdp.start2
          >> bdp.end2
          >> bdp.qname
          >> bdp.mapq
          >> bdp.strand1
          >> bdp.strand2;
  }

  if(bdpin.fail() || bdpin.bad())
  {
    std::cerr<<"ERROR: reading BDProw from "
             <<bdp_file<<"\n"
             <<"be sure that your bedpe file is tab delimited\n"
             <<"with the following columns:\n"
             <<"chr1 start1 end1 chr2 start2 end2 read_name min_mapQ strand1 strand2"
             <<std::endl;

    exit(EXIT_FAILURE);

  }

  bdpin.ignore(1000,'\n');

  bdp.start1 ++; // bedpe start is 0-based, turn it to 1-based
  bdp.start2 ++; // for handling

  bdp.samfile=NULL;
  bdp.bdpfile=this;
  totalReads++;

}

//----------------------------------- loadBDPs -------------------------------//
void BDPfile::loadBDPs()
{
  
  BDProw bdp;
  while(bdpin.good() && bdpin.peek() != EOF){
    nextBDP(bdp);
    mpairs.push_back(bdp);
  }

}

//----------------------------------- check_sort -----------------------------//
bool BDPfile::check_sort ()
{

  std::vector<BDProw>::iterator it;

  for(it=mpairs.begin(); it<mpairs.end()-2; it++)
    if(!compareBDProws(*it, *(it+1)))
      return(false);

  return(true);

}

//----------------------------------- sort -----------------------------------//
void BDPfile::sort()
{

  if(!check_sort())
    std::sort(mpairs.begin(),mpairs.end(),compareBDProws);

}

//----------------------------------- printReadReport ------------------------//
std::ostream& BDPfile::printReadReport(std::ostream& out)
{

  out << "\n\n----- Read Report for: "<<bdp_file<<std::endl
      << "---------------------------------" << std::endl
      << "Total read pairs:                  "<< totalReads << std::endl
      << "Zstat filtered:                    "<< filtFz << std::endl
      << "Low mapping quality:               "<< filtMapq << std::endl
      << "Map to different contigs:          "<< filtJmpChrom << std::endl
      << "Final read pairs:                  "<< totalReads - filtFz - filtMapq - filtJmpChrom
      << std::endl <<std::endl; 

}

//----------------------------------- BDPfile(char* file) --------------------//
BDPfile::BDPfile (char * bdpf)
{

  bdp_file=std::string(bdpf);
  bdpin.open(bdpf);

  printFiltered=false;
  filtFz=0;
  filtMapq=0;
  filtJmpChrom=0;
  totalReads=0;

}

//----------------------------------- BDPfile(string file) -------------------//
BDPfile::BDPfile(std::string bdpf)
{

  bdp_file=bdpf;
  bdpin.open(bdpf.c_str());

  printFiltered=false;
  filtFz=0;
  filtMapq=0;
  filtJmpChrom=0;
  totalReads=0;

}

//===================================================== SAMfile ==============//
//----------------------------------- nextSAM --------------------------------//
SAMrow SAMfile::nextSAM()
{

  SAMrow sam;

  if(!samin.eof() && 
     !samin.fail() && 
     !samin.bad() &&
     samin.peek() != EOF){

    while(samin.peek()=='@') //ignore SAM file headers
      samin.ignore(10000,'\n');

    samin >> sam.qname
          >> sam.flag
          >> sam.rname
          >> sam.start
          >> sam.mapq
          >> sam.cigar
          >> sam.mrname
          >> sam.mstart
          >> sam.inslen
          >> sam.seq
          >> sam.rqual;

    if(samin.fail() || samin.bad()){
      std::cerr<<"ERROR reading SAM file"<<std::endl;
      exit(EXIT_FAILURE);
    }

    if(samin.peek()!='\n'){
      samin.get(); //get rid of TAB
      getline(samin,sam.opt);
    }

    //Parse sam.qname, erase "/1" or "/2"
    std::size_t found;
    found=sam.qname.rfind("/");
    if(found != std::string::npos &&
       sam.qname.size()>found+1 &&
       (sam.qname.at(found+1) == '1' ||
        sam.qname.at(found+1) == '2' )){
      sam.qname.erase(found);
    }

  }

  return(sam);

}

//----------------------------------- nextBDP --------------------------------//
BDProw SAMfile::nextBDP()
{

  BDProw bdp;
  SAMrow read1, read2;

  if(!samin.eof() &&
     !samin.fail() &&
     !samin.bad() &&
     samin.peek() != EOF ){
 
    read2=nextSAM();

    do{
      read1=read2;
      read2=nextSAM();
      filtSingleReads++;
    }while(read1.qname != read2.qname && samin.good() && samin.peek() != EOF);

    if(read1.qname != read2.qname && (!samin.good() || samin.peek() == EOF)){
      std::cerr<<"ERROR reading read pairs from SAM file"<<std::endl;
      exit(EXIT_FAILURE);
    }

    filtSingleReads--;
    
    if(read1.start <= read2.start){

      bdp.chr1=read1.rname;
      bdp.chr2=read2.rname;
  
      bdp.start1=read1.start;
      bdp.start2=read2.start;
  
      bdp.mapq=std::min(read1.mapq, read2.mapq);

      bdp.end1=read1.start + read1.seq.length();
      bdp.end2=read2.start + read2.seq.length();
  
      bdp.strand1=(read1.flag & 16) ? '-' : '+';
      bdp.strand2=(read2.flag & 16) ? '-' : '+';
  
    } else {
  
      bdp.chr1=read2.rname;
      bdp.chr2=read1.rname;
  
      bdp.start1=read2.start;
      bdp.start2=read1.start;

      bdp.mapq=std::min(read1.mapq, read2.mapq);

      bdp.end1=read2.start + read2.seq.length();
      bdp.end2=read1.start + read1.seq.length();

      bdp.strand1=(read2.flag & 16) ? '-' : '+';
      bdp.strand2=(read1.flag & 16) ? '-' : '+';
    }

    bdp.qname=read1.qname;
  
  }

  bdp.bdpfile=NULL;
  bdp.samfile=this;
  totalReads++;

  return(bdp);

}

//----------------------------------- writeBDP -------------------------------//
std::ostream& SAMfile::writeBDP(std::ostream& out, int minMapQ)
{

  BDProw bdp;

  while(samin.good() && samin.peek()!=EOF){
    bdp=nextBDP();
    if(printFiltered || !bdp.filterOut(minMapQ))
      bdp.print(out);
  }
}

//----------------------------------- writeBDPsort ---------------------------//
std::ostream& SAMfile::writeBDPsort(std::ostream& out, int minMapQ)
{

  BDProw bdp;
  BDPfile bdpfile;

  while(samin.good() && samin.peek() != EOF){
    bdp=nextBDP();
    if(!bdp.filterOut(minMapQ)){
      bdpfile.mpairs.push_back(bdp);
    }
  }
  bdpfile.sort();
  if(! bdpfile.print(out)){
    std::cerr << "ERROR writing bedpe.sort file\n";
    exit(EXIT_FAILURE);
  }
  bdpfile.mpairs.clear(); 
  return(out);
}

//----------------------------------- loadSeqLens ----------------------------//
bool SAMfile::loadSeqLens(std::map<std::string, int >& ctglens)
{

  if(samin.peek() != '@')
    return(false);

  std::string seqName, header;
  int seqLength;

  while(samin.peek() == '@'){
    samin >> header;
    if(header=="@SQ"){
      samin >> header;
      seqName=header.substr(header.find(":")+1,header.size());
      samin >> header;
      header=header.substr(header.find(":")+1,header.size());
      seqLength=atoi(header.c_str());
      ctglens[seqName]=seqLength;
      samin.ignore(10000,'\n');
    } else //other headers
      samin.ignore(10000,'\n');
  }
  return(true);
}

//----------------------------------- printReadReport ------------------------//
std::ostream& SAMfile::printReadReport(std::ostream& out)
{

  out << "\n\n----- Read Report for: "<<sam_file<<std::endl
      << "---------------------------------" << std::endl
      << "Single reads:                      "<<filtSingleReads << std::endl
      << "Total read pairs:                  "<< totalReads << std::endl
      << "Low mapping quality:               "<< filtMapq << std::endl
      << "Map to different contigs:          "<< filtJmpChrom << std::endl
      << "Final read pairs:                  "<< totalReads - filtMapq - filtJmpChrom
      << std::endl << std::endl;

}

//----------------------------------- SAMfile(char * file) -------------------//
SAMfile::SAMfile(char * file)
{

  printFiltered=false;
  filtSingleReads=0;
  filtMapq=0;
  filtJmpChrom=0;
  totalReads=0;

  sam_file=std::string(file);
  samin.open(file);
 
}

//----------------------------------- SAMfile(string file) -------------------//
SAMfile::SAMfile(std::string file)
{

  printFiltered=false;
  filtSingleReads=0;
  filtMapq=0;
  filtJmpChrom=0;
  totalReads=0;

  sam_file=file;
  samin.open(file.c_str());

}

//===================================================== mateAnRow ============//
//----------------------------------- N,zneg,zpos,zneut ----------------------//
int mateAnRow::getN()
{
  return(N);
}

int mateAnRow::getZneg()
{
  return(zneg);
}

int mateAnRow::getZpos()
{
  return(zpos);
}

int mateAnRow::getZneut()
{
  return(getN()-getZneg()-getZpos());
}

//----------------------------------- goodOrient, badOrient ------------------//
int mateAnRow::getGoodOrient()
{
  return(goodOrientCount);
} 

int mateAnRow::getBadOrient()
{
  return(getN()-getGoodOrient());
}

//----------------------------------- sum,start,end,size ---------------------//
int mateAnRow::getSum()
{
  return(sum);
}

int mateAnRow::getStart()
{
  return(start);
}

int mateAnRow::getEnd()
{
  return(end);
}

//----------------------------------- zvalue ---------------------------------//
double mateAnRow::zvalue(double mu0, double sd0)
{

  if(N>0)
    return( (double(sum)/double(N) - mu0) / (sd0/sqrt(N)) );
  else
    return(0);

}

//----------------------------------- newZstat -------------------------------//
double mateAnRow::newZstat(double mu0, double sd0, int delLength, int insLength)
{
  return(((double(sum)/double(N)) - delLength + insLength - mu0) / (sd0/sqrt(N)));
}

//----------------------------------- newDeltaMean ---------------------------//
double mateAnRow::newDeltaMean(double mu0, int delLength, int insLength)
{
  return(-1*(double(sum)/double(N) - mu0 - delLength + insLength));
}

//----------------------------------- print ----------------------------------//
std::ostream& mateAnRow::print(std::ostream& out, std::string chr, double mu0, double sd0)
{
  out << chr << "\t"
      << getStart()-1 << "\t" //bedpe start is 0-based
      << getEnd() << "\t"
      << getN() << "\t"
      << getZneg() << "\t"
      << getZneut() << "\t"
      << getZpos() << "\t"
      << getSum() << "\t"
      << zvalue(mu0,sd0) << "\t"
      << getGoodOrient() << "\t"
      << getBadOrient() << std::endl;
}

//===================================================== mateAn ===============//
//----------------------------------- meanZstat ------------------------------//
double mateAn::meanZstat()
{
  return(double(sumZstat)/double(totalBases));
}

//----------------------------------- sdZstat --------------------------------//
double mateAn::sdZstat()
{
  return(sqrt(double(sumZstatSquare)/double(totalBases) - pow(meanZstat(),2)));
}

//----------------------------------- meanCov --------------------------------//
double mateAn::meanCov()
{ 
  return(double(sumCov)/double(totalBases));
}

//----------------------------------- sdCov ----------------------------------//
double mateAn::sdCov()
{
  return(sqrt(double(sumCovSquare)/double(totalBases) - pow(meanCov(),2)));
}

//----------------------------------- suggestMu0 -----------------------------//
double mateAn::suggestMu0()
{
  return(mu0 + meanZstat() * sd0 / sqrt(meanCov()));
}

//----------------------------------- print ----------------------------------//
std::ostream& mateAn::print(std::ostream& out)
{
  std::map< std::string, std::vector<mateAnRow> >::iterator mit;
  std::vector<mateAnRow>::iterator vit;

  out << "#@mu0=" << mu0 << std::endl
      << "#@sd0=" << sd0 << std::endl;

  for(mit = mateAns.begin(); mit != mateAns.end(); mit++)
    for(vit = mit -> second.begin(); vit < mit -> second.end(); vit++)
      vit -> print(out,mit->first,mu0,sd0);

  return(out);
}

//----------------------------------- search ---------------------------------//
int mateAn::search(std::string chr, int pos)
{
  //Return -1 if the position chr-pos is not found.
  //Otherwise return the index "i" such that mateAns[chr][i] contains
  //the zstat information for position "pos"

  if(mateAns.count(chr) == 0){
    return(-1);
  }

  std::vector<mateAnRow>& vec= mateAns[chr];
  mateAnRow x;

  int start,end,middle;

  start=0;
  end=vec.size();
  middle=(end-start)/2;

  do{
    x=vec[middle];
    if(x.start <= pos && x.end >= pos){
      return(middle);
    }

    if(pos > x.end)
      start=middle+1;
    else if(pos < x.start)
      end=middle-1;

    middle=(end+start)/2;
  }while(start <= end && middle >= 0 && middle <= (int)vec.size()-1);
  return(-1);
}

//----------------------------------- isRangeGoodCov ---------------------------//
bool mateAn::isRangeGoodCov(std::string chr, int start, int end, int minCov)
{

  bool rangeGoodCov;
  int coord1;
  int coord2;

  start=std::min(start,end);
  end=std::max(start,end);

  coord1=search(chr,start);
  coord2=search(chr,end);

  rangeGoodCov = coord1>=0 && mateAns[chr][coord1].N >= minCov &&  
                 coord2>=0 && mateAns[chr][coord2].N >= minCov;

  if(rangeGoodCov){
    int cerange;
    int prevEnd=mateAns[chr][coord1].end;
    for(cerange=coord1+1; cerange <= coord2; cerange++){
      if(mateAns[chr][cerange].start = prevEnd+1){
        rangeGoodCov=rangeGoodCov && mateAns[chr][cerange].N >= minCov;
      }else{
        return(false);
      }
    }
    return(rangeGoodCov);
  }else{
    return(false);
  }
}

//----------------------------------- updateStats ----------------------------//
void mateAn::updateStats(double mu0, double sd0)
{

  totalRanges=0;
  totalBases=0;
  sumZstat=0;
  sumZstatSquare=0;
  minZstat=100000;
  maxZstat=0;
  sumCov=0;
  sumCovSquare=0;
  minCov=100000;
  maxCov=0;

  std::map< std::string, std::vector<mateAnRow> >::iterator mit;
  std::vector<mateAnRow>::iterator vit;
  double zstat;
  int cov;

  for(mit = mateAns.begin(); mit != mateAns.end(); mit++)
    for(vit = mit -> second.begin(); vit < mit -> second.end(); vit++){

      if(vit -> getN() < MIN_N_STAT)
        continue;
      int rangeLength = vit -> end  -  vit -> start + 1;
      totalRanges++;
      totalBases += rangeLength;
      sumZstat += vit -> zvalue(mu0,sd0) * rangeLength;
      sumZstatSquare += pow(vit -> zvalue(mu0,sd0), 2) * rangeLength;
      minZstat = std::min(minZstat, vit -> zvalue(mu0,sd0));
      maxZstat = std::max(maxZstat, vit -> zvalue(mu0,sd0));
      sumCov += vit -> getN() * rangeLength;
      sumCovSquare += pow(vit -> getN(), 2) * rangeLength;
      minCov = std::min(minCov, vit -> getN());
      maxCov = std::max(maxCov, vit -> getN());
    }

}

void mateAn::updateStats()
{
  updateStats(mu0,sd0);
}

//----------------------------------- printStats -----------------------------//
std::ostream& mateAn::printStats(std::ostream& out)
{

  out << "---------- mateAn Stats ----------" << std::endl
      << "----------------------------------" << std::endl
      << "Cov cutoff:             " << MIN_N_STAT << std::endl
      << "Mu0:                    " << mu0 << std::endl
      << "Sd0:                    " << sd0 << std::endl
      << "Total ranges:           " << totalRanges << std::endl
      << "Total bases:            " << totalBases << std::endl
      << "Mean CE-stat:           " << meanZstat() << std::endl
      << "Sd CE-stat:             " << sdZstat() << std::endl
      << "Min CE-stat:            " << minZstat << std::endl
      << "Max CE-stat:            " << maxZstat << std::endl
      << "Mean Cov:               " << meanCov() << std::endl
      << "Sd Cov:                 " << sdCov() << std::endl
      << "Min Cov:                " << minCov << std::endl
      << "Max Cov:                " << maxCov << std::endl << std::endl
      << "Suggested mu0:          " << suggestMu0()
      << std::endl << std::endl;

}

//----------------------------------- addMateAnRow ---------------------------//
void mateAn::addMateAnRow(std::string chr, mateAnRow& mateRow)
{

  int rangeLength= mateRow.end - mateRow.start + 1;
  mateAns[chr].push_back(mateRow);
   if(mateRow.getN() >= MIN_N_STAT){
     totalRanges++;
     totalBases+= rangeLength;
     sumZstat += mateRow.zvalue(mu0,sd0) * rangeLength;
     sumZstatSquare += pow(mateRow.zvalue(mu0,sd0),2) * rangeLength;
     minZstat = std::min(minZstat, mateRow.zvalue(mu0,sd0));
     maxZstat = std::max(maxZstat, mateRow.zvalue(mu0,sd0));
     sumCov += mateRow.getN() * rangeLength;
     sumCovSquare += pow(mateRow.getN(),2) * rangeLength;
     minCov = std::min(minCov, mateRow.getN());
     maxCov = std::max(maxCov, mateRow.getN());
   }

}

//----------------------------------- readMateAn -----------------------------//
void mateAn::readMateAn(std::ifstream& in, std::string prefix)
{

  char c1,c2;
  std::string header,chr, disregardLine;
  double value;
  size_t pos;
  mateAnRow mateRow;
  int disregardInt;
  double disregardDouble;

  totalBases=0;
  totalRanges=0;
  sumZstat=0;
  sumZstatSquare=0;
  minZstat=1000000;
  maxZstat=0;
  sumCov=0;
  sumCovSquare=0;
  minCov=1000000;
  maxCov=0;

  do{

    in.get(c1);
    in.get(c2);

    if(!(c1=='#' && c2=='@')){
      in.putback(c2);
      in.putback(c1);
      break;
    }

    in.putback(c2);
    in.putback(c1);
    std::getline(in,header);
    if(header.find("#@") != std::string::npos){
      pos=header.find("sd0=");
      if(pos != std::string::npos){
        sd0=atof(header.substr((int)pos+4).c_str());
      }
      pos=header.find("mu0=");
      if(pos != std::string::npos){
        mu0=atof(header.substr((int)pos+4).c_str());
      }
    }
  }while(in.peek()=='#');

  while(in.good() && in.peek() != EOF){
    
    in >> chr
       >> mateRow.start
       >> mateRow.end
       >> mateRow.N
       >> mateRow.zneg
       >> disregardInt
       >> mateRow.zpos
       >> mateRow.sum
       >> disregardDouble
       >> mateRow.goodOrientCount
       >> disregardInt;

    if(in.bad()){
      std::cerr<<"Can't read mateAn file. Be sure that it as a TAB\n"
               <<"delmited table with the following columns:\n"
               <<"1)chrom,2)start,3)end,4)Coverage,5)zneg,6)zneut,\n"
               <<"7)zpos,8)ins_sum,9)zvalue,10)goodOrient,11)badOrient\n";
       exit(EXIT_FAILURE);
     }

    std::getline(in,disregardLine);

    chr=std::string(prefix).append(chr);
    mateRow.start+1; //Bedpe start is 0-based  
    addMateAnRow(chr,mateRow);
  }

  if(mateAns.size() <= 0){
    std::cerr<<"No ce-stat data? No can't do"<<std::endl;
    exit(EXIT_FAILURE);
  } 
}

//----------------------------------- clear() --------------------------------//
void mateAn::clear()
{
  totalBases=0;
  totalRanges=0;
  sumZstat=0;
  sumZstatSquare=0;
  minZstat=1000000;
  maxZstat=0;
  sumCov=0;
  sumCovSquare=0;
  minCov=1000000;
  maxCov=0;
  mateAns.clear();

}

//----------------------------------- mateAn() -------------------------------//
mateAn::mateAn(std::ifstream& in, std::string prefix)
{
  mateAns.clear();
  mu0=0;
  sd0=1;
  readMateAn(in,prefix);
}

mateAn::mateAn(std::ifstream& in)
{
  mateAns.clear();
  mu0=0;
  sd0=1;
  readMateAn(in,std::string(""));
}

mateAn::mateAn(std::string file, std::string prefix)
{
  std::ifstream inf (file.c_str(),std::ifstream::in);
  mateAns.clear();
  mu0=0;
  sd0=1;
  readMateAn(inf,prefix);
}

mateAn::mateAn (double mu, double sd) 
{
  totalRanges=0;
  totalBases=0; 
  sumZstat=0;
  sumZstatSquare=0;
  minZstat=1000000;
  maxZstat=0;
  sumCov=0;
  sumCovSquare=0;
  minCov=1000000;
  maxCov=0;
  mu0=mu;
  sd0=sd;
  MIN_N_STAT=30;
  mateAns.clear();
}

//---------------------------------- zsignal ---------------------------------//
char mateAn::zsignal(std::string scf, int pos, char type){
 
  if(mateAns.count(scf) == 0){
   /* cerr<<"--------\n"
        <<"Error: Zsignal: There is no scaffold with name "<<scf<<"\n"
        <<"in file "<<file<<"\n"
        <<"Position: "<<pos<<"\n"
        <<"Type: "<<type<<"\n\n";
    exit(1);*/
    return('X');
  }

  std::vector<mateAnRow >&mateAn = mateAns[scf];

  int found=search(scf,pos);

  if(found < 0){
    return('X');
  }

  if(mateAn[found].getN() < MIN_COV){
   return('C');
  }

  switch(type){
    case 'i':
      if(mateAn[found].getBadOrient() >= MIN_MIS_OR){
        return('i');
      }
      if(fabs(mateAn[found].zvalue(mu0,sd0)) < ZSTAT_THRE ){
        return('n');
      }
      if(mateAn[found].zvalue(mu0,sd0) <= -1*ZSTAT_THRE){
        if(mateAn[found].getZneg() >= MIN_SUP){
          return('-');
        } else {
          return('n');
        }
      }
      if(mateAn[found].zvalue(mu0,sd0) >= ZSTAT_THRE){
        if(mateAn[found].getZpos() >= MIN_SUP){
          return('+');
        } else {
          return('n');
        }
      }
      return('0');
      break;

    case 'z':
      if(mateAn[found].zvalue(mu0,sd0) <= -1*ZSTAT_THRE){
        if(mateAn[found].getZneg() >= MIN_SUP){
          return('-');
        }else {
          return('n');
        }
      }
      if(mateAn[found].zvalue(mu0,sd0) >= ZSTAT_THRE){
	if(mateAn[found].getZpos() >= MIN_SUP){
          return('+');
        }else {
          return('n');
        }
      }
      if(mateAn[found].getBadOrient() >= MIN_MIS_OR){
        return('i');
      }
      if(fabs(mateAn[found].zvalue(mu0,sd0)) < ZSTAT_THRE ){
        return('n');
      }
      return('0');
      break;
  }
  return('U');
}

//---------------------------------- zvalue ----------------------------------//
double mateAn::zvalue(std::string scf, int pos)
{

  int found=search(scf,pos);

  if(found < 0)
    return(0);
  else
    return(mateAns[scf][found].zvalue(mu0,sd0)); 

}
