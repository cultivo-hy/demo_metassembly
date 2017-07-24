///////////////////////////////////////////////////////////////////////
// 
// author: Alejandro Hernandez Wences
// date: 21 Feb 2013
// 
// Class declarations for merging assemblies
//
// see merge.cc
///////////////////////////////////////////////////////////////////////

#include "merge.hh"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <string>
#include <map>
#include <vector>

int SMALL_SEQ=10;

//===================================================== DeltaFile ============//
bool sR_Sort(const AlignStats& pA, const AlignStats& pB){return(pA.sR < pB.sR);}
bool sQ_Sort(const AlignStats& pA, const AlignStats& pB){return(pA.sQ < pB.sQ);} 
//----------------------------------- parseDeltaR ----------------------------//
void DeltaFile::parseDeltaR(std::map< std::string, std::vector<AlignStats> > & Aligns, std::string Id1)
{

  AlignStats aStats;                     //  single alignment region
  bool found = false;

  DeltaReader_t dr;
  dr.open (DeltaName);

  while ( dr.readNext() ){
    if ( dr.getRecord().idR == Id1 ){
      found = true;
      for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ ){
        aStats.sR = dr.getRecord( ).aligns[i].sR;
        aStats.eR = dr.getRecord( ).aligns[i].eR;
        aStats.sQ = dr.getRecord( ).aligns[i].sQ;
        aStats.eQ = dr.getRecord( ).aligns[i].eQ;

        aStats.Delta = dr.getRecord( ).aligns[i].deltas;

        //-- Add the new alignment
        Aligns[dr.getRecord().idQ].push_back(aStats);
      }
    }
  }

  if ( !found ){
    std::cerr<<"ERROR: Could not find any alignments for "<<Id1<<std::endl;
    exit (EXIT_FAILURE);
  }
  dr.close( );
  return;
}

//----------------------------------- parseDeltaQ ----------------------------//
void DeltaFile::parseDeltaQ(std::map< std::string, std::vector<AlignStats> > & Aligns, std::string Id1)
{
  AlignStats aStats;                     //  single alignment region
  bool found = false;

  DeltaReader_t dr;
  dr.open (DeltaName);

  while ( dr.readNext() ){
    if ( dr.getRecord().idQ == Id1 ){
      found = true;
      for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ ){
        aStats.sR = dr.getRecord( ).aligns[i].sR;
        aStats.eR = dr.getRecord( ).aligns[i].eR;
        aStats.sQ = dr.getRecord( ).aligns[i].sQ;
        aStats.eQ = dr.getRecord( ).aligns[i].eQ;

        aStats.Delta = dr.getRecord( ).aligns[i].deltas;

        //-- Add the new alignment
        Aligns[dr.getRecord().idR].push_back(aStats);
      }
    }
  }

  if ( !found ){
    std::cerr<<"ERROR: Could not find any alignments for "<<Id1<<std::endl;
    exit (EXIT_FAILURE);
  }
  dr.close( );
  return;
} 

//----------------------------------- transCoords ----------------------------//
std::vector<int> DeltaFile::transCoords(std::vector<AlignStats> Aligns, std::string Id1, int pos1, std::string Id2, bool isReferenceCoord)
{
  std::vector<AlignStats>::iterator Ap;
  std::vector<int>::iterator Dp;
  std::vector<int>found;

  int strandA;
  int strandB;
  int currPosA;
  int currPosB;
  int endA;
  int endB;
  int pos2;
  int dif;
  int insLen;

  for(Ap=Aligns.begin(); Ap<Aligns.end(); Ap++){

    if(isReferenceCoord){
      currPosA=Ap->sR;
      endA=Ap->eR;
      currPosB=Ap->sQ;
      endB=Ap->eQ;
    } else{
      currPosA=Ap->sQ;
      endA=Ap->eQ;
      currPosB=Ap->sR;
      endB=Ap->eR;
    }

    strandA=(currPosA > endA)?-1:1;
    strandB=(currPosB > endB)?-1:1;

    if(pos1 >= std::min(currPosA,endA) && pos1 <= std::max(currPosA,endA)){
      insLen=0;

      if(!isReferenceCoord)
        for(Dp=Ap->Delta.begin(); Dp < Ap->Delta.end(); Dp++)
          (*Dp)=(*Dp)*-1;

      for(Dp=Ap->Delta.begin(); Dp < Ap->Delta.end(); Dp++){

        int delta=*Dp;

        if(labs(delta)>1)
          insLen=1;
        else if(delta !=0)
          insLen+=1;
        
        if(delta==0){
          dif=labs(pos1 - currPosA);
          pos2=currPosB + dif * strandB;
          found.push_back(pos2);
          break;
        }

        currPosA+=(labs(delta)-1)*strandA;
        currPosB+=(labs(delta)-1)*strandB;

        if(strandA*currPosA > pos1*strandA){
          currPosA-=(labs(delta)-1)*strandA;
          currPosB-=(labs(delta)-1)*strandB;          
          dif=labs(pos1 - currPosA);
          pos2=currPosB + dif * strandB; 
          found.push_back(pos2);
          break;
        }else if(currPosA == pos1){//pos1 falls in ref insertion
          pos2=currPosB-1*strandB;
          found.push_back(pos2);
          break;
        }
        if(delta > 0)
          currPosA+=1*strandA;
        else if(delta <0)
          currPosB+=1*strandB;
      }
    }
  }
  return(found);
}

//----------------------------------- transCoords ----------------------------//
int DeltaFile::transCoords(std::vector<AlignStats> Aligns, 
                               std::string Id1, int pos1, 
                               std::string Id2, int posBs, int posBe, 
                               bool isReferenceCoord)
{

  std::vector<AlignStats>::iterator Ap;
  std::vector<int>::iterator Dp;
  int found=0;

  int strandA;
  int strandB;
  int currPosA;
  int currPosB;
  int endA;
  int endB;
  int pos2;
  int dif;
  int insLen;

  for(Ap=Aligns.begin(); Ap<Aligns.end(); Ap++){

    if(isReferenceCoord){
      currPosA=Ap->sR;
      endA=Ap->eR;
      currPosB=Ap->sQ;
      endB=Ap->eQ;
    } else{
      currPosA=Ap->sQ;
      endA=Ap->eQ;
      currPosB=Ap->sR;
      endB=Ap->eR;
    }

    if(currPosB != posBs || endB != posBe)
      continue;

    strandA=(currPosA > endA)?-1:1;
    strandB=(currPosB > endB)?-1:1;

    if(pos1 >= std::min(currPosA,endA) && pos1 <= std::max(currPosA,endA)){
      insLen=0;

      if(!isReferenceCoord)
        for(Dp=Ap->Delta.begin(); Dp < Ap->Delta.end(); Dp++)
          (*Dp)=(*Dp)*-1;

      for(Dp=Ap->Delta.begin(); Dp < Ap->Delta.end(); Dp++){

        int delta=*Dp;

        if(abs(delta)>1)
          insLen=1;
        else if(delta !=0)
          insLen+=1;
        
        if(delta==0){
          dif=abs(pos1 - currPosA);
          pos2=currPosB + dif * strandB;
          return(pos2);
        }

        currPosA+=(labs(delta)-1)*strandA;
        currPosB+=(labs(delta)-1)*strandB;

        if(strandA*currPosA > pos1*strandA){
          currPosA-=(labs(delta)-1)*strandA;
          currPosB-=(labs(delta)-1)*strandB;          
          dif=labs(pos1 - currPosA);
          pos2=currPosB + dif * strandB; 
          return(pos2);
          break;
        }else if(currPosA == pos1){//pos1 falls in ref insertion
          pos2=currPosB-1*strandB;
          return(pos2);
        }
        if(delta > 0)
          currPosA+=1*strandA;
        else if(delta <0)
          currPosB+=1*strandB;
      }
    }
  }
  return(found);
}

//----------------------------------- query2ref ------------------------------//
int DeltaFile::query2ref(std::string qscf, int qpos, std::string rscf, int rposA, int rposB)
{

  std::map< std::string, std::vector<AlignStats> > Aligns;
  int found=0;
  parseDeltaQ(Aligns, qscf);

  std::map< std::string, std::vector<AlignStats> >::iterator it;
  for(it=Aligns.begin(); it != Aligns.end(); it++)
    std::sort (it->second.begin( ), it->second.end( ), sR_Sort); 

  for(it=Aligns.begin(); it != Aligns.end(); it++)
    if(it->first == rscf)
      return(transCoords (it->second, qscf, qpos, it->first, rposA, rposB, false));

  return(found);
}

//----------------------------------- query2ref ------------------------------//
std::map< std::string, std::vector<int> > DeltaFile::query2ref(std::string qscf, int qpos)
{

  std::map< std::string, std::vector<AlignStats> > Aligns;
  std::map< std::string, std::vector<int> >found;
  parseDeltaQ(Aligns, qscf);

  std::map< std::string, std::vector<AlignStats> >::iterator it;
  for(it=Aligns.begin(); it != Aligns.end(); it++)
    std::sort (it->second.begin( ), it->second.end( ), sR_Sort); 

  for(it=Aligns.begin(); it != Aligns.end(); it++)
    found[it->first]=transCoords (it->second, qscf, qpos, it->first, false);

  return(found);
}

//----------------------------------- query2ref ------------------------------//
std::vector<int> DeltaFile::query2ref(std::string qscf, int qpos, std::string rscf)
{

  std::map< std::string, std::vector<int> >found;
  found=query2ref(qscf,qpos);
  if(found.count(rscf) > 1)
    return(found[rscf]);
  else
    return(std::vector<int>());
}

//----------------------------------- ref2query ------------------------------//
int DeltaFile::ref2query(std::string rscf, int rpos, std::string qscf, int qposA, int qposB)
{

  std::map< std::string, std::vector<AlignStats> > Aligns;
  int found=0;
  parseDeltaR(Aligns, qscf);

  std::map< std::string, std::vector<AlignStats> >::iterator it;
  for(it=Aligns.begin(); it != Aligns.end(); it++)
    std::sort (it->second.begin( ), it->second.end( ), sQ_Sort); 

  for(it=Aligns.begin(); it != Aligns.end(); it++)
    if(it->first == qscf)
      return(transCoords (it->second, rscf, rpos, it->first, qposA, qposB, false));

  return(found);
}

//----------------------------------- ref2query ------------------------------//
std::map< std::string, std::vector<int> > DeltaFile::ref2query(std::string rscf, int rpos)
{

  std::map< std::string, std::vector<AlignStats> > Aligns;
  std::map< std::string, std::vector<int> >found;
  parseDeltaQ(Aligns, rscf);

  std::map< std::string, std::vector<AlignStats> >::iterator it;
  for(it=Aligns.begin(); it != Aligns.end(); it++)
    std::sort (it->second.begin( ), it->second.end( ), sQ_Sort); 

  for(it=Aligns.begin(); it != Aligns.end(); it++)
    found[it->first]=transCoords (it->second, rscf, rpos, it->first, true);

  return(found);
  
}

//----------------------------------- ref2query ------------------------------//
std::vector<int> DeltaFile::ref2query(std::string rscf, int rpos, std::string qscf)
{
  std::map< std::string, std::vector<int> >found;
  found=ref2query(rscf,rpos);
  if(found.count(qscf) > 1)
    return(found[qscf]);
  else
    return(std::vector<int>());
}

//===================================================== CoordsRow ============//
//----------------------------------- print ----------------------------------//
std::ostream& CoordsRow::print(std::ostream& out)
{
  out << start1 << "\t"
      << end1 << "\t"
      << start2 << "\t"
      << end2 << "\t"
      << alnsize1 << "\t"
      << alnsize2 << "\t"
      << ident << "\t"
      << scfLen1 << "\t"
      << scfLen2 << "\t"
      << cover1 << "\t"
      << cover2 << "\t"
      << scf1 << "\t"
      << scf2 << "\t"
      << std::endl;
  return(out);
}

//----------------------------------- CoordsRow() ----------------------------//
CoordsRow::CoordsRow(std::ifstream& MUM)
{

  if(MUM.good() && MUM.peek() != EOF){
    MUM >> start1
        >> end1
        >> start2
        >> end2
        >> alnsize1
        >> alnsize2
        >> ident
        >> scfLen1
        >> scfLen2
        >> cover1
        >> cover2
        >> scf1
        >> scf2;

    scf1=std::string("R.").append(scf1);
    scf2=std::string("Q.").append(scf2);

    if(MUM.bad()){
      std::cerr<< "Be sure that your file is the output of show-coords\n"
          << "and contains a TAB delimited table with the following\n"
          << "clumns:\n"
          << "start1 end1 start2 end2 alnlength1 alnlength2 ident scfLen1 scfLen2 cover1 cover2 scf1 scf2\n"
          << "The file is the output of:\n"
          << "show-coords -H -c -l -r -T\n";
      exit(EXIT_FAILURE);
    }
    MUM.ignore(10000,'\n');
  }
}

//===================================================== MUMfile ==============//
//----------------------------------- compareOrdCoords -----------------------//
bool compareOrdCoords(ordCoords A, ordCoords B){
  if(A.scf == B.scf){
    int Astart=std::min(A.start,A.end);
    int Aend=std::max(A.start,A.end);
    int Bstart=std::min(B.start,B.end);
    int Bend=std::max(B.start,B.end);
    if(Astart == Bstart){
      return(Aend < Bend);
    } else{
      return(Astart < Bstart);
    }
  } else {
    return(A.scf < B.scf);
  }
}

//----------------------------------- loadFile -------------------------------//
std::ifstream& MUMfile::loadFile()
{

  std::vector<ordCoords> refCoords, queCoords;
  ordCoords refOrd, queryOrd;
  int line=0;

  while(coordsin.good() && coordsin.peek() != EOF){
    CoordsRow mumrow(coordsin);
    coords.push_back(mumrow);    
    line++;

    refOrd.scf=mumrow.scf1;
    refOrd.start=mumrow.start1;
    refOrd.end=mumrow.end1;
    refOrd.line=line-1;
    refCoords.push_back(refOrd);

    queryOrd.scf=mumrow.scf2;
    queryOrd.start=mumrow.start2;
    queryOrd.end=mumrow.end2;
    queryOrd.line=line-1;
    queCoords.push_back(queryOrd);
  }
 
  sort(queCoords.begin(), queCoords.end(), compareOrdCoords);
  sort(refCoords.begin(), refCoords.end(), compareOrdCoords);

  std::vector<ordCoords>::iterator it;
  
  for(it=refCoords.begin(); it != refCoords.end(); it++)
    refScfs[it -> scf].push_back(it -> line);

  for(it=queCoords.begin(); it != queCoords.end(); it++)
    queryScfs[it -> scf].push_back(it -> line);
}

//----------------------------------- printMums ------------------------------//
void MUMfile::printMums(std::string out_pre)
{

  std::ofstream OR(std::string(out_pre).append(".mum.orig").c_str());
  std::ofstream REF(std::string(out_pre).append(".mum.refOrd").c_str());
  std::ofstream QUE(std::string(out_pre).append(".mum.queOrd").c_str());

  std::vector<CoordsRow>::iterator row;

  for(row = coords.begin(); row < coords.end(); row++)
    row -> print(OR);

  std::map<std::string, std::vector<int> >::iterator ord;
  std::vector<int>::iterator rowo;
  for(ord = refScfs.begin(); ord != refScfs.end(); ord++)
    for(rowo = ord -> second.begin(); rowo < ord -> second.end(); rowo++)
      coords[*rowo].print(REF);

  for(ord = queryScfs.begin(); ord != queryScfs.end(); ord++)
    for(rowo = ord -> second.begin(); rowo < ord -> second.end(); rowo++)
      coords[*rowo].print(QUE);
}

//----------------------------------- findQueEdgeAln -------------------------//
std::vector<int>::iterator MUMfile::findQueEdgeAln(std::string scf, std::string edge)
{
  std::vector<int>& alns = queryScfs[scf];
  if(edge=="end"){
    std::vector<int>::reverse_iterator rit;
    for(rit=alns.rbegin(); rit<alns.rend();rit++){
      if(coords[*rit].alnsize2 >= minOverAlnLength && 
         coords[*rit].ident*coords[*rit].alnsize2 >= minBasesAligned &&
         coords[*rit].scfLen2 - std::max(coords[*rit].end2,coords[*rit].start2) <= maxEdgeNotAligned
         ){
           return((rit+1).base());
	 }
    }  
  } else if(edge=="start"){
    std::vector<int>::iterator it;
    for(it=alns.begin(); it<alns.end(); it++){
      if(coords[*it].alnsize2 >= minOverAlnLength && 
         coords[*it].ident*coords[*it].alnsize2 >= minBasesAligned &&
         std::min(coords[*it].start2,coords[*it].end2) <= maxEdgeNotAligned
         ){
           return(it);
         } 	
    }	  
  }
  return(alns.end());
}

//----------------------------------- findRefEdgeAln -------------------------//
std::vector<int>::iterator MUMfile::findRefEdgeAln(std::string scf, std::string edge)
{
  std::vector<int>&alns =  refScfs[scf];
  if(edge=="end"){
    std::vector<int>::reverse_iterator rit;
    for(rit=alns.rbegin(); rit<alns.rend();rit++){
      if(coords[*rit].alnsize1 >= minOverAlnLength && 
         coords[*rit].ident*coords[*rit].alnsize1 >= minBasesAligned &&
         coords[*rit].scfLen1 - std::max(coords[*rit].end1,coords[*rit].start1) <= maxEdgeNotAligned
         ){
            return((rit+1).base());
         } 	
    }  
  } else if(edge=="start"){
    std::vector<int>::iterator it;
    for(it=alns.begin(); it<alns.end(); it++){
      if(coords[*it].alnsize1 >= minOverAlnLength && 
         coords[*it].ident*coords[*it].alnsize1 >= minBasesAligned &&
         std::min(coords[*it].start1,coords[*it].end1) -1 <= maxEdgeNotAligned
         ){
           return(it);
         } 	
    }	  
  }
  return(alns.end());
}

//----------------------------------- adjacentAlnQue()
std::vector<int>::iterator MUMfile::adjacentAlnQue(std::vector<int>&alns, 
                                                   std::vector<int>::iterator pos, 
                                                   char direction)
{

  std::vector<int>::iterator quePos;
  for(quePos = alns.begin(); quePos <= alns.end(); quePos++){
    if(quePos == alns.end())
      break;
    if(*quePos == *pos)
      break;
  }

  if(quePos == alns.end()){
    std::cerr<<"Error: adjacentAln:\n"
        <<"Reference MumRow:\n"
        <<*pos << "\t";
    coords[*pos].print(std::cerr);    
    std::cerr<<"\n"
        <<"Query MumRows:\n";
    std::vector<int>::iterator it;
    for(it=alns.begin();it<alns.end();it++){
      std::cerr<<(*it)<<"\t";
      coords[*it].print(std::cerr);
      std::cerr<<"\n";
    }
    std::cerr<<"\n";
    exit(EXIT_FAILURE);
  }

  int startQue1=std::min(coords[*quePos].start2,coords[*quePos].end2);
  int endQue1=std::max(coords[*quePos].start2,coords[*quePos].end2); 

  std::vector<int>::iterator it;
  if(direction=='+'){
    for(it=quePos+1; it<alns.end();it++){
      if(coords[*it].alnsize2 >= minOverAlnLength &&
         coords[*it].ident * coords[*it].alnsize2 >= minBasesAligned &&
         !(coords[*it].start2 >= startQue1 && coords[*it].start2 <= endQue1 &&
           coords[*it].end2 >= startQue1 && coords[*it].end2 <= endQue1)
        ){
        return(it);  
      }
    }
  } else if(direction=='-'){
    for(it=quePos-1;it>=alns.begin();it--){
      if(coords[*it].alnsize2 >= minOverAlnLength &&
         coords[*it].ident * coords[*it].alnsize2 >= minBasesAligned  &&
         !(coords[*it].start2 >= startQue1 && coords[*it].start2 <= endQue1 &&
           coords[*it].end2 >= startQue1 && coords[*it].end2 <= endQue1)
        ){
        return(it);  
      }
    }
  } else {
    std::cerr<< "Error: adjacentAln: durection must be either \"+\" or \"-\"";
    exit(EXIT_FAILURE);
  }
  return(alns.end());
}

//----------------------------------- adjacentAlnRef -------------------------//
std::vector<int>::iterator MUMfile::adjacentAlnRef(std::vector<int>&alns,
                                          std::vector<int>::iterator pos,
                                          char direction)
{
  std::vector<int>::iterator refPos;
  for(refPos = alns.begin(); refPos <= alns.end(); refPos++){
    if(refPos == alns.end())
      break;
    if(*refPos == *pos)
      break;
  }
  if(refPos == alns.end()){
    std::cerr<<"Error: adjacentAln:\n"
        <<"Reference MumRow:\n"
        <<*pos << "\t";
    coords[*pos].print(std::cerr);    
    std::cerr<<"\n"
        <<"Query MumRows:\n";
    std::vector<int>::iterator it;
    for(it=alns.begin();it<alns.end();it++){
      std::cerr<<(*it)<<"\t";
      coords[*it].print(std::cerr);
      std::cerr<<"\n";
    }
    std::cerr<<"\n";
    exit(EXIT_FAILURE);
  }

  std::vector<int>::iterator it; 
  if(direction=='+'){
    for(it=refPos+1; it<alns.end();it++){
      if(!(coords[*it].start1>=coords[*refPos].start1 &&
              coords[*it].end1 <= coords[*refPos].end1)){
        return(it);  
      }
    }
  } else if(direction=='-'){
    for(it=refPos-1;it>=alns.begin();it--){
      if(!(coords[*it].start1>=coords[*refPos].start1 &&
              coords[*it].end1 <= coords[*refPos].end1)){
        return(it);  
      }
    }
  } else {
    std::cerr<< "Error: adjacentAln: direction must be either \"+\" or \"-\"";
    exit(EXIT_FAILURE);
  }
  return(alns.end());
}

//----------------------------------- QscfBefore -----------------------------//
bool MUMfile::QscfBefore(std::string scf, std::vector<int>&alns, std::vector<int>::iterator last)
{
  std::vector<int>::iterator i;
  for(i = alns.begin(); i < last; i++)
    if(coords[*i].scf2 == scf)
      return(true);
  return(false);
}

//----------------------------------- QscfAfter ------------------------------//
std::vector<int>::iterator MUMfile::QscfAfter(std::string scf, 
                                     std::vector<int>& alns, 
                                     std::vector<int>::iterator begin)
{
  std::vector<int>::iterator i;
  for(i = begin+1; i < alns.end(); i++)
    if(! (coords[*i].start1 >= coords[*begin].start1 &&
          coords[*i].end1 <= coords[*begin].end1) &&
       coords[*i].scf2 == scf)
      return(i);

  return(alns.end());
}

//----------------------------------- analyseInterInsert ---------------------//
bool MUMfile::analyseInterInsert(std::vector<int>& alns,
                                 std::vector<int>::iterator A,
                                 std::vector<int>::iterator B,
                                 mateAn& ce)
{
  std::vector<int>::iterator it;
  int alnCount=0;
  bool foundSig=false;
  for(it = A+1; it < B; it++){
    alnCount++;
    if( coords[*it].start1 - coords[*A].end1  > 1000 &&
        coords[*B].start1 - coords[*it].start1 > 1000 &&
        ce.zsignal(coords[*it].scf1, coords[*it].start1, 'z') != 'n'
      )
    {
      foundSig=true;
      if(alnCount>=3){
        return(false);
      }
    } else if (coords[*it].end1 - coords[*A].end1 > 1000 &&
               coords[*B].start1 - coords[*it].end1 > 1000 &&
               ce.zsignal(coords[*it].scf1, coords[*it].end1, 'z') != 'n'
              )
    {
      foundSig=true;
      if(alnCount>=3){
        return(false);
      }
    }
  } 
  if(alnCount>=3 && foundSig){
    return(false);
  }

  return(true);
}

//----------------------------------- MUMfile(std::string in_coords) ---------//
MUMfile::MUMfile(std::string in_coords, std::string inDelta, int iminOver, int iminBas, int imaxEdge)
{
  coordsin.open(in_coords.c_str());
  delta.DeltaName=inDelta;
  loadFile();
  minOverAlnLength=iminOver;
  minBasesAligned=iminBas;
  maxEdgeNotAligned=imaxEdge;
}


//===================================================== Nsfile ===============// 
//----------------------------------- loadFile -------------------------------//
std::ifstream& Nsfile::loadFile(std::ifstream& in, std::string prefix)
{

  std::vector< int>coords(2,0);
  std::string scf;

  while(in.peek() != EOF){
    in >> scf
         >> coords[0]
         >> coords[1];

    if(!in.eof() && in.fail()){
      std::cerr<<"Besure that your Ns file has a TAB delimited\n"
               <<"table with the following columns:\n"
               <<"scf start end"
               <<std::endl;
      exit(EXIT_FAILURE);
    }

    scf=std::string(prefix).append(scf);
    coords[0]+=1; //Ns file is bed, start is 0-based

    if(Ndata.count(scf)>0)
      if(Ndata[scf].back()[1] == coords[0]+1)
        //if range is consecutive merge them
        Ndata[scf].back()[1]=coords[1];
      else 
        Ndata[scf].push_back(coords);  
    else 
      Ndata[scf].push_back(coords);
    
  }
}

//----------------------------------- search ---------------------------------//
int Nsfile::search(std::string scf, int pos1, int pos2)
{
  //Return -1 if the coordinates scf:pos1-pos2 do not overlap any N-insert
  //range in Ndata. Otherwise return the index "i" such that the range Ndata[scf][i]
  //overlaps the range scf:pos1-pos2

  if(Ndata.count(scf) < 1)
    return(-1);

  std::vector< std::vector< int> >& ranges=Ndata[scf];
  int start,end,middle;

  start=0;
  end=ranges.size()-1;

  middle=(end-start)/2;

  int overlap;

  do{
    overlap=std::min(ranges[middle][1],std::max(pos1,pos2))-std::max(ranges[middle][0],std::min(pos1,pos2));
    if(overlap>=0){
      return(middle);
    }
    if(std::min(pos1,pos2) > ranges[middle][1]){
      start=middle+1;
    } else if(std::max(pos1,pos2) < ranges[middle][0]){
      end=middle-1;
    }
    middle=(end+start)/2;
  }while(start<=end && middle>=0 && middle<=(int)ranges.size());
  return(-1);   
}

//----------------------------------- isNinsert ------------------------------//
bool Nsfile::isNinsert(std::string scf,
               int pos1,
               int pos2,
               int minNinsert,
               double minNinsPer
              )
{

  int found=search(scf,pos1,pos2);
  int it;
  std::vector<int>overlapping;
  overlapping.push_back(found);
  if(found>=0){
    std::vector< std::vector< int> >& ranges=Ndata[scf];
    //Find all N ranges that overlap with [std::min(pos1,pos2),std::max(pos1,pos2)]
    for(it=found-1; it>=0 && std::min(ranges[it][1],std::max(pos1,pos2))-std::max(ranges[it][0],std::min(pos1,pos2))>=0; it--)
      overlapping.push_back(it);
    
    for(it=found+1; it<=(int)ranges.size()-1 && std::min(ranges[it][1],std::max(pos1,pos2))-std::max(ranges[it][0],std::min(pos1,pos2))>=0; it++)
      overlapping.push_back(it);
    
    //Determine whether there is at least one row of Ns with length at least minN
    //also compute the total number of Ns present in [std::min(pos1,pos2),std::max(pos1,pos2)]
    bool meetMinN=false;
    std::vector<int>::iterator r;
    int sumNs=0;
    for(r=overlapping.begin(); r<overlapping.end(); r++){
      sumNs+=std::min(ranges[*r][1],std::max(pos1,pos2))-std::max(ranges[*r][0],std::min(pos1,pos2))+1;
      if(std::min(ranges[*r][1],std::max(pos1,pos2))-std::max(ranges[*r][0],std::min(pos1,pos2))+1 >= minNinsert)
        meetMinN=true;
    }
    if(meetMinN && double(sumNs)/double(std::max(pos1,pos2)-std::min(pos1,pos2)+1) >= minNinsPer){
      //cout<<"file:"<<file<<"\tscf:"<<scf<<"\tstart:"<<pos1<<"\tend:"<<pos2<<"\n";
      return(true);
    }
  }
  return(false);
}

//---------------------------------- Nsfile() --------------------------------//
Nsfile::Nsfile(std::string filename, std::string prefix)
{
  ns_file=filename;
  nsin.open(ns_file.c_str());
  loadFile(nsin, prefix);
}

//===================================================== metaLNcoords =========//
//---------------------------------- print -----------------------------------//
void metaLNcoords::print(std::string metaScf, int mstart, int mend,  
                         MetaRow* mrow, Metassembly* meta,std::ofstream& out
                        )
{
  out<<"------------------------------------------"<<std::endl;
  out<<"----Link Coordinates:"<<std::endl;
  mln.print(out);
  int coordindex=meta->mum.refScfs[mln.refScf1][0];
  out<<mln.refScf1<<" length:"<<meta->mum.coords[coordindex].scfLen1<<std::endl;
  coordindex=meta->mum.queryScfs[mln.queScf][0];
  out<<mln.queScf<<" length:"<<meta->mum.coords[coordindex].scfLen2<<std::endl;
  coordindex=meta->mum.refScfs[mln.refScf2][0];
  out<<mln.refScf2<<" length:"<<meta->mum.coords[coordindex].scfLen1<<std::endl;
  
  out<<"----Reference 1 ce-stat data:"<<std::endl;
  int found;
  found=meta->Rce.search(mln.refScf1,mln.refStart1);
  if(found >= 0)
    meta->Rce[mln.refScf1][found].print(out,mln.refScf1,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<mln.refScf1<<"\t"<<mln.refStart1<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(mln.refScf1,mln.refEnd1);
  if(found >= 0)
    meta->Rce[mln.refScf1][found].print(out,mln.refScf1,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<mln.refScf1<<"\t"<<mln.refEnd1<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Query 1 ce-stat data:"<<std::endl;
  found=meta->Qce.search(mln.queScf,mln.queStart1);
  if(found >= 0)
    meta->Qce[mln.queScf][found].print(out,mln.queScf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<mln.queScf<<"\t"<<mln.queStart1<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Qce.search(mln.queScf,mln.queEnd1);
  if(found >= 0)
    meta->Qce[mln.queScf][found].print(out,mln.queScf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<mln.queScf<<"\t"<<mln.queEnd1<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Query 2 ce-stat data:"<<std::endl;
  found=meta->Qce.search(mln.queScf,mln.queStart2);
  if(found >= 0)
    meta->Qce[mln.queScf][found].print(out,mln.queScf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<mln.queScf<<"\t"<<mln.queStart2<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Qce.search(mln.queScf,mln.queEnd2);
  if(found >= 0)
    meta->Qce[mln.queScf][found].print(out,mln.queScf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<mln.queScf<<"\t"<<mln.queEnd2<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Reference 2 ce-stat data:"<<std::endl;
  found=meta->Rce.search(mln.refScf2,mln.refStart2);
  if(found >= 0)
    meta->Rce[mln.refScf2][found].print(out,mln.refScf2,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<mln.refScf2<<"\t"<<mln.refStart2<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(mln.refScf2,mln.refEnd2);
  if(found >= 0)
    meta->Rce[mln.refScf2][found].print(out,mln.refScf2,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<mln.refScf2<<"\t"<<mln.refEnd2<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Description and metassembly info: ";//<<std::endl;
  out<<desc<<"\n";
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";
  if((orientation == mrow->orientation && edge == "start") ||
     (orientation != mrow->orientation && edge == "end")
    ){
    delta=(orientation == mrow->orientation)?delta:-1*delta;
    out<<metaScf<<"\t"<<mstart + delta <<std::endl;  
  }else{
    delta=(orientation == mrow->orientation)?delta:-1*delta;
    out<<metaScf<<"\t"<<mend + delta <<std::endl;
  }
  out<<"------------------------------------------\n"<<std::endl;
}

//===================================================== metaNAcoords =========//
//---------------------------------- print -----------------------------------//
void metaNAcoords::print(std::string metaScf, int mstart, int mend, 
                         MetaRow* mrow, Metassembly* meta,std::ofstream& out)
{
  int found;
  out<<"------------------------------------------"<<std::endl;
  out<<"----Reference N-insert: "<<rscf<<"\t"<<rnstart<<"\t"<<rnend
                                 <<". Length:\t"<<rnend-rnstart+1<<std::endl;
  found=meta->Rce.search(rscf,rnstart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rnstart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,rnend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else 
    out<<rscf<<"\t"<<rnend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Query "<<qoside<<" overhang: "<<qscf<<"\t"<<qostart<<"\t"<<qoend
                   <<". Length:\t"<<std::max(qostart,qoend) - std::min(qostart,qoend) + 1<<std::endl;
  found=meta->Qce.search(qscf,qostart);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qostart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;
  
  found=meta->Qce.search(qscf,qostart);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qoend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Metassembly info: ";//<<std::endl;
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";

  if(qor == mrow->orientation )
    out<<metaScf<<"\t"<<mstart<<std::endl
       <<metaScf<<"\t"<<mend<<std::endl;
  else
    out<<metaScf<<"\t"<<mend<<std::endl
       <<metaScf<<"\t"<<mstart<<std::endl;
  out<<"------------------------------------------\n"<<std::endl;
}

//===================================================== metaNBcoords =========//
//---------------------------------- print -----------------------------------//
void metaNBcoords::print(std::string metaScf, int mstart, int mend,
                         MetaRow* mrow, Metassembly* meta, std::ofstream& out)
{
  int found;  
  out<<"------------------------------------------"<<std::endl;
  if(qoedge == "left"){
    out<<"----Reference left N-insert: "<<rscf<<"\t"<<rnstart<<"\t"<<rnend
       <<". Length: "<<std::max(rnstart,rnend)-std::min(rnstart,rnend) + 1<<std::endl;
                                        
    found=meta->Rce.search(rscf,rnstart);
    if(found >= 0)
      meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
    else
      out<<rscf<<"\t"<<rnstart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

    found=meta->Rce.search(rscf,rnend);
    if(found >= 0)
      meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
    else
      out<<rscf<<"\t"<<rnend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

    out<<"----Reference unmapped right sequence: "<<rscf<<"\t"<<ristart<<"\t"<<riend
       <<". Lenght: "<<std::max(ristart,riend)-std::min(ristart,riend) + 1<<std::endl;
    found=meta->Rce.search(rscf,riend);
    if(found >= 0)
      meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
    else
      out<<rscf<<"\t"<<riend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;      

  }else{

    out<<"----Reference unmapped left sequence: "<<rscf<<"\t"<<ristart<<"\t"<<riend
       <<". Length: "<<std::max(ristart,riend)-std::min(ristart,riend) + 1<<std::endl;
    found=meta->Rce.search(rscf,ristart);
    if(found >= 0)
      meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
    else
      out<<rscf<<"\t"<<ristart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

    out<<"----Reference right N-insert: "<<rscf<<"\t"<<rnstart<<"\t"<<rnend
       <<". Length: "<<std::max(rnstart,rnend)-std::min(rnstart,rnend) + 1<<std::endl;
    found=meta->Rce.search(rscf,rnstart);
    if(found >= 0)
      meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
    else
      out<<rscf<<"\t"<<rnstart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

    found=meta->Rce.search(rscf,rnend);
    if(found >= 0)
      meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
    else
      out<<rscf<<"\t"<<rnend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;
  }

  out<<"----Query "<<qoedge<<" overhang: "<<qscf<<"\t"<<qostart<<"\t"<<qoend
     <<". Length: "<<std::max(qostart,qoend)-std::min(qostart,qoend) + 1<<std::endl;
  found=meta->Qce.search(qscf,qostart);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qostart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Qce.search(qscf,qoend);
  if(found >= 0) 
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qoend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Metassembly info: ";//<<std::endl;
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";

  if((qor == mrow->orientation && qoedge == "left") || 
     (qor != mrow->orientation && qoedge == "right")
    )
    out<<metaScf<<"\t"<<mstart<<std::endl
       <<metaScf<<"\t"<<mend<<std::endl;
  else
    out<<metaScf<<"\t"<<mend<<std::endl
       <<metaScf<<"\t"<<mstart<<std::endl;

  out<<"------------------------------------------\n"<<std::endl;
}

//===================================================== metaNCcoords =========//
//---------------------------------- print -----------------------------------//
void metaNCcoords::print(std::string metaScf, int mstart, int mend,
                         MetaRow* mrow, Metassembly* meta, std::ofstream& out)
{
  int found;
  out<<"------------------------------------------"<<std::endl;
  out<<"----Query left overhang: "<<qscf1<<"\t"<<qlostart<<"\t"<<qloend
     <<". Length: "<<std::max(qlostart,qloend)-std::min(qlostart,qloend) + 1<<std::endl;
  found=meta->Qce.search(qscf1,qlostart);
  if(found >= 0)
    meta->Qce[qscf1][found].print(out,qscf1,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf1<<"\t"<<qlostart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;  

  found=meta->Qce.search(qscf1,qloend);
  if(found >= 0)
    meta->Qce[qscf1][found].print(out,qscf1,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf1<<"\t"<<qloend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Reference left N-insert: "<<rscf<<"\t"<<rlnstart<<"\t"<<rlnend
     <<". Length: "<<std::max(rlnstart,rlnend)-std::min(rlnstart,rlnend) + 1<<std::endl;
  found=meta->Rce.search(rscf,rlnstart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rlnstart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,rlnend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rlnend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Reference unmapped insert: "<<rscf<<"\t"<<ristart<<"\t"<<riend
     <<". Length: "<<std::max(ristart,riend)-std::min(ristart,riend) + 1<<std::endl;
  found=meta->Rce.search(rscf,ristart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<ristart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,riend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<riend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Reference right N-insert: "<<rscf<<"\t"<<rrnstart<<"\t"<<rrnend
     <<". Length: "<<std::max(rrnstart,rrnend)-std::min(rrnstart,rrnend) + 1<<std::endl;
  found=meta->Rce.search(rscf,rrnstart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rrnstart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,rrnend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rrnend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;
 
  out<<"----Query right overhang: "<<qscf2<<"\t"<<qrostart<<"\t"<<qroend
     <<". Length: "<<std::max(qrostart,qroend)-std::min(qrostart,qroend) + 1<<std::endl;
  found=meta->Qce.search(qscf2,qrostart);
  if(found >= 0)
    meta->Qce[qscf2][found].print(out,qscf2,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf2<<"\t"<<qrostart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Qce.search(qscf2,qroend);
  if(found >= 0)
    meta->Qce[qscf2][found].print(out,qscf2,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf2<<"\t"<<qroend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Metassembly info: ";//<<std::endl;
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";

  if(mrow->orientation != '+'
    )
    out << metaScf << "\t" << mstart <<std::endl
        << metaScf << "\t" << mend << std::endl;
  else 
    out << metaScf << "\t" << mend << std::endl
        << metaScf << "\t" << mstart << std::endl;
  out<<"------------------------------------------\n"<<std::endl;
}

//===================================================== metaN1coords =========//
//---------------------------------- print -----------------------------------//
void metaN1coords::print(std::string metaScf, int mstart, int mend,
                         MetaRow* mrow, Metassembly* meta, std::ofstream& out)
{
  out<<"------------------------------------------"<<std::endl;
  out<<"----Reference insert: "<<rscf<<"\t"<<ristart<<"\t"<<riend
     <<". Length: "<<std::max(ristart,riend)-std::min(ristart,riend) + 1<<std::endl;
  int found=meta->Rce.search(rscf,ristart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<ristart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,riend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<riend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Reference N-insert: "<<rscf<<"\t"<<rnstart<<"\t"<<rnend
     <<". Length: "<<std::max(rnstart,rnend)-std::min(rnstart,rnend) + 1<<std::endl;
  found=meta->Rce.search(rscf,rnstart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rnstart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,rnend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rnend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Query insert: "<<qscf<<"\t"<<qistart<<"\t"<<qiend
     <<". Length: "<<std::max(qistart,qiend)-std::min(qistart,qiend) + 1<<std::endl;
  found=meta->Qce.search(qscf,qistart);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qistart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;  

  found=meta->Qce.search(qscf,qiend);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qiend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Metassembly info: ";//<<std::endl;
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";

  if(qor == mrow->orientation)
    out<<metaScf<<"\t"<<mstart<<std::endl
       <<metaScf<<"\t"<<mend<<std::endl;
  else
    out<<metaScf<<"\t"<<mend<<std::endl
       <<metaScf<<"\t"<<mstart<<std::endl;

  out<<"------------------------------------------\n"<<std::endl;
}


//===================================================== metaQIcoords =========//
//---------------------------------- print -----------------------------------//
void metaQIcoords::print(std::string metaScf, int mstart, int mend,
                         MetaRow* mrow, Metassembly* meta, std::ofstream& out)
{
  out<<"------------------------------------------"<<std::endl;
  out<<"----Reference insert: "<<rscf<<"\t"<<ristart<<"\t"<<riend
     <<". Length: "<<std::max(0, 
                              std::max(ristart,riend)-std::min(ristart,riend) + 1)<<std::endl;
  int found=meta->Rce.search(rscf,ristart);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<ristart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Rce.search(rscf,riend);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<riend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Query insert: "<<qscf<<"\t"<<qistart<<"\t"<<qiend;
  if(metaScf == "Corr_ref_insertion"){
    out<<". Length: 0"<<std::endl;
  } else {
    out<<". Length: "<<std::max(1,std::max(qistart,qiend)-std::min(qistart,qiend) + 1)<<std::endl;
  }
  found=meta->Qce.search(qscf,qistart);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qistart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Qce.search(qscf,qiend);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qiend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Metassembly info: ";//<<std::endl;
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";

  if(qor == mrow->orientation)
    out<<metaScf<<"\t"<<mstart<<std::endl
       <<metaScf<<"\t"<<mend<<std::endl;
  else
    out<<metaScf<<"\t"<<mend<<std::endl
       <<metaScf<<"\t"<<mstart<<std::endl;
  out<<"------------------------------------------\n"<<std::endl;
}

//===================================================== metaBKPcoords ========//
//---------------------------------- print -----------------------------------//
void metaBKPcoords::print(std::string metaScf, int mstart, int mend,
                          MetaRow* mrow, Metassembly* meta, std::ofstream& out)
{
  out<<"------------------------------------------"<<std::endl;
  out<<"----Reference breakpoint: "<<rscf<<"\t"<<rbkp<<std::endl;
  int found=meta->Rce.search(rscf,rbkp);
  if(found >= 0)
    meta->Rce[rscf][found].print(out,rscf,meta->Rce.mu0,meta->Rce.sd0);
  else
    out<<rscf<<"\t"<<rbkp<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Query "<<qoside<<" overhang: "<<qscf<<"\t"<<qostart<<"\t"<<qoend
     <<". Length: "<<std::max(qostart,qoend)-std::min(qostart,qoend) + 1<<std::endl;
  found=meta->Qce.search(qscf,qostart);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qostart<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  found=meta->Qce.search(qscf,qoend);
  if(found >= 0)
    meta->Qce[qscf][found].print(out,qscf,meta->Qce.mu0,meta->Qce.sd0);
  else
    out<<qscf<<"\t"<<qoend<<"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"<<std::endl;

  out<<"----Metassembly info: ";//<<std::endl;
  out<<metaScf<<"\t"
     <<mstart<<"\t"
     <<mend<<"\t"
     <<mrow->file<<"\t"
     <<mrow->scf<<"\t"
     <<mrow->start<<"\t"
     <<mrow->end<<"\t"
     <<mrow->orientation<<"\n";

  if((qoside == "left" && qor == mrow->orientation) || 
     (qoside == "right" && qor != mrow->orientation)
    ) 
    out<<metaScf<<"\t"<<mstart<<std::endl
       <<metaScf<<"\t"<<mend<<std::endl;
  else
    out<<metaScf<<"\t"<<mend<<std::endl
       <<metaScf<<"\t"<<mstart<<std::endl;
  out<<"------------------------------------------\n"<<std::endl;
}


//===================================================== MetaRow ==============//
//----------------------------------- changeOrientatio -----------------------//
void MetaRow::changeOrientation()
{
  if(orientation == '+')
    orientation='-';
  else if(orientation == '-')
    orientation='+';
}

//===================================================== MetaLink =============//
//----------------------------------- print ----------------------------------//
std::ostream& MetaLink::print(std::ostream& out)
{
  out << refScf1 << "\t"
      << refStart1 << "\t"
      << refEnd1 << "\t"
      << edgeScf1 << "\t"
      << queScf << "\t"
      << queStart1 << "\t"
      << queEnd1 << "\t"
      << queOrient1 << "\t"
      << queStart2 << "\t"
      << queEnd2 << "\t"
      << queOrient2 << "\t"
      << refScf2 << "\t"
      << refStart2 << "\t"
      << refEnd2 << "\t"
      << edgeScf2
      << std::endl;

  return(out);
}

//===================================================== Metassembly ==========//
//----------------------------------- addMetaRow -----------------------------//
MetaRow* Metassembly::addMetaRow(std::string mscf,
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
                             int aln)
                            
{

  MetaRow* ret;

  MetaRow Coord;
  Coord.scf=seqScf;
  Coord.start=seqStart;
  Coord.end=seqEnd;
  Coord.file=file;
  Coord.orientation=orientation;
  Coord.sigRef1=sigRef1;
  Coord.sigRef2=sigRef2;
  Coord.sigQue1=sigQue1;
  Coord.sigQue2=sigQue2;
  Coord.aln=aln;
  Coord.isInv=false;
  Coord.mlnc=NULL;
  Coord.mNAc=NULL;
  Coord.mNBc=NULL;
  Coord.mNCc=NULL;
  Coord.mN1c=NULL;
  Coord.mQIc=NULL;
  Coord.mBKPc=NULL;

  if(edge == "end"){
    metassem[mscf].push_back(Coord);
    ret=&(*(metassem[mscf].end()-1));
  } else if(edge=="start"){
    std::vector<MetaRow>::iterator it;
    it=metassem[mscf].begin();
    metassem[mscf].insert(it,Coord);
    ret=&(*(metassem[mscf].begin()));
  } else {
    std::cerr<<"Error in addMetaRow. edge is not \"end\" nor \"start\""<<std::endl;
    exit(EXIT_FAILURE);
  }
  metaLengths[mscf]+=abs(Coord.end-Coord.start)+1;
  return(ret);
}

//----------------------------------- findScfLink -----------------------------//
std::string Metassembly::findScfLink(std::vector<int>::iterator refAln, 
                                     std::string edgeScf1 
                                    )
{
  std::string queScf=mum.coords[*refAln].scf2;
  std::vector<int>& queAlns=mum.queryScfs[queScf];
  char scf1OrientQuery=(mum.coords[*refAln].start2 <= mum.coords[*refAln].end2)? '+':'-';
  char scf2OrientQuery;

  //Get the correspondig adj query aln
  std::vector<int>::iterator adjQueAln;
  if(edgeScf1 == "start"){
    if(scf1OrientQuery == '+'){
      adjQueAln=mum.adjacentAlnQue(queAlns,refAln,'-');
    } else { //scf1OrientQuery=='-'
      adjQueAln=mum.adjacentAlnQue(queAlns,refAln,'+');
    }
  } else{ //edgeScf1=="end"
    if(scf1OrientQuery == '+'){
      adjQueAln=mum.adjacentAlnQue(queAlns,refAln,'+');
    } else { //scf1OrientQuery=='-'
      adjQueAln=mum.adjacentAlnQue(queAlns,refAln,'-');
    }
  }

  if(adjQueAln == queAlns.end())
  //queAlns.end() is the flag for the absence of a good adjQueAln
    return("");

  if(mum.coords[*adjQueAln].scf1 != mum.coords[*refAln].scf1  
     //Disregard when overhang aligns to the same ref scf
    ){

    scf2OrientQuery=(mum.coords[*adjQueAln].start2<=mum.coords[*adjQueAln].end2)?'+':'-';
    std::vector<int>::iterator firstScf2aln=mum.findRefEdgeAln(mum.coords[*adjQueAln].scf1,"start"); //first aln of ref scf2
    std::vector<int>::iterator lastScf2aln=mum.findRefEdgeAln(mum.coords[*adjQueAln].scf1, "end"); //last aln of ref scf2

    if( firstScf2aln != mum.refScfs[mum.coords[*adjQueAln].scf1].end() &&
        *adjQueAln == *firstScf2aln &&/* start edge of scf2 */
        ((edgeScf1 == "end" && scf2OrientQuery == scf1OrientQuery) ||
        (edgeScf1 == "start" && scf2OrientQuery != scf1OrientQuery))
      ){
        //Scf1 links to Scf2. Now corroborate that Scf2 links to Scf1.
        std::vector<int>::iterator Scf2adjQueAln;
        if(scf2OrientQuery == '+')
          Scf2adjQueAln=mum.adjacentAlnQue(queAlns,firstScf2aln,'-');
        else  //scf1OrientQuery=='-'
          Scf2adjQueAln=mum.adjacentAlnQue(queAlns,firstScf2aln,'+');

        //Get query coverage
        bool covQ=false;
        std::vector<int> QinserPos(2,0);
        mateAnRow QZstat1;
        mateAnRow QZstat2;
        QinserPos[0]=std::min(std::max(mum.coords[*refAln].start2,mum.coords[*refAln].end2),
                         std::max(mum.coords[*adjQueAln].start2,mum.coords[*adjQueAln].end2));
        QinserPos[1]=std::max(std::min(mum.coords[*refAln].start2,mum.coords[*refAln].end2),
                         std::min(mum.coords[*adjQueAln].start2,mum.coords[*adjQueAln].end2));

        int QZstat1i=Qce.search(mum.coords[*refAln].scf2, QinserPos[0]);
        if(QZstat1i >= 0)
          QZstat1=Qce[mum.coords[*refAln].scf2][QZstat1i];
        bool Qcov1=(QZstat1i >= 0 && QZstat1.N >= minLinkCov)?true:false;

        int QZstat2i=Qce.search(mum.coords[*refAln].scf2, QinserPos[1]);
        if(QZstat2i >= 0)
          QZstat2=Qce[mum.coords[*refAln].scf2][QZstat2i];
        bool Qcov2=(QZstat2i >= 0 && QZstat2.N >=minLinkCov)?true:false;

        if(QinserPos[1] - QinserPos[0] -1 >= 1 && Qcov1 && Qcov2){          
          covQ=true;
          for(int i=QZstat1i; i <= QZstat2i; i++){
            covQ=covQ && Qce[mum.coords[*refAln].scf2][i].N >= minLinkCov;
          }
        }else
          covQ=Qcov1 && Qcov2; //If the queAlns overlap in the query, it is sufficient
                               //that one edge has a good coverage.

        if(*Scf2adjQueAln == *refAln && covQ){
          addLink(*refAln,edgeScf1,*firstScf2aln,"start");
          countLN++;
          return(mum.coords[*adjQueAln].scf1);
        }else
          return("");

    } else if( lastScf2aln != mum.refScfs[mum.coords[*adjQueAln].scf1].end() &&
               *adjQueAln == *lastScf2aln && /* end edge of scf2 */
               ((edgeScf1 == "end" && scf2OrientQuery != scf1OrientQuery ) ||
               (edgeScf1 == "start" && scf2OrientQuery == scf1OrientQuery)) 
             ){

        //Scf1 links to Scf2. Now corroborate that Scf2 links to Scf1.
        std::vector<int>::iterator Scf2adjQueAln;
        if(scf2OrientQuery == '+')
          Scf2adjQueAln=mum.adjacentAlnQue(queAlns,lastScf2aln,'+');
        else  //scf1OrientQuery=='-'
          Scf2adjQueAln=mum.adjacentAlnQue(queAlns,lastScf2aln,'-');

        //Get query coverage
        bool covQ=false;
        std::vector<int> QinserPos(2,0);
        mateAnRow QZstat1;
        mateAnRow QZstat2;
        QinserPos[0]=std::min(std::max(mum.coords[*refAln].start2,mum.coords[*refAln].end2),
                         std::max(mum.coords[*adjQueAln].start2,mum.coords[*adjQueAln].end2));
        QinserPos[1]=std::max(std::min(mum.coords[*refAln].start2,mum.coords[*refAln].end2),
                         std::min(mum.coords[*adjQueAln].start2,mum.coords[*adjQueAln].end2));

        int QZstat1i=Qce.search(mum.coords[*refAln].scf2,QinserPos[0]);
        if(QZstat1i >= 0){
          QZstat1=Qce[mum.coords[*refAln].scf2][QZstat1i];
        }
        bool Qcov1=(QZstat1i >= 0 && QZstat1.N >= minLinkCov)?true:false;

        int QZstat2i=Qce.search(mum.coords[*refAln].scf2,QinserPos[1]);
        if(QZstat2i >= 0){
          QZstat2=Qce[mum.coords[*refAln].scf2][QZstat2i];
        }
        bool Qcov2=(QZstat2i >= 0 && QZstat2.N >=minLinkCov)?true:false;

        if(QinserPos[1] - QinserPos[0] -1 >= 1 && Qcov1 && Qcov2){
          covQ=true;
          for(int i=QZstat1i; i <= QZstat2i; i++){
            covQ=covQ && Qce[mum.coords[*refAln].scf2][i].N >= minLinkCov;
          }
        }else
          covQ=Qcov1 && Qcov2; //If the queAlns overlap in the query, it is sufficient
                               //that one edge has a good coverage.

        if(*Scf2adjQueAln == *refAln && covQ){
          addLink(*refAln,edgeScf1,*lastScf2aln,"end");
          countLN++;
          return(mum.coords[*adjQueAln].scf1);
        } else
          return("");

    } else { // endQueEdge aln is not the first nor the last aln of scf2.
              // No link can be established.
      return("");
    }
  } else { //query scf returns to the same ref scf. Can't link with itself, yet.
    return("");
  }
}

//----------------------------------- joinLinks ------------------------------//
void Metassembly::joinLinks()
{
  std::map<std::string, MetaLink>::iterator lnIt;
  std::string metaScf1,metaScf2;

  for(lnIt=scfsLinks.begin(); lnIt != scfsLinks.end(); lnIt++){
    MetaLink& ln=(*lnIt).second;
    metaScf1=Scf_Meta[ln.refScf1];
    metaScf2=Scf_Meta[ln.refScf2];

//    std::cerr<<"\n\n######################################################################################"<<std::endl;
//    std::cerr<<"Link: "<<std::endl;
//    ln.print(std::cerr);
//    std::cerr<<"######################################################################################"<<std::endl;

    if(metaScf1==metaScf2){
      //each link will be twice in scfsLinks.
      //If a link has already be processed then the corresponding changes
      //in Scf_Meta will have been made already, and metaScf1==metaScf2
      //In this case, the link does not need to be processed. 
      continue;
    }

//    std::cerr<<"Analyse"<<std::endl;
    //Initialize variables
    bool refScf1Over=false; //scf1 overhang?
    bool refScf2Over=false; //scf2 overhang?
    bool queInter=false;  // sequence in the middle of qScf that is not aligned

    int refScf1OverStart=ln.refStart1;
    int refScf1OverEnd=ln.refEnd1;
    int refScf2OverStart=ln.refStart2;
    int refScf2OverEnd=ln.refEnd2;
    int queInterStart;
    int queInterEnd;

    char refScf1OverSupp='l';
    char refScf2OverSupp='l';
    char queInterStartSupp='l';
    char queInterEndSupp='l';

    int startScf2copy;
    int startScf1copy;

//    std::cerr<<"Find edges"<<std::endl;

    std::string edgeMeta1=findEdgeInMeta(metaScf1,ln.refScf1,(ln.edgeScf1=="start")?ln.refStart1:ln.refEnd1,ln,ln.edgeScf1);
    std::string edgeMeta2=findEdgeInMeta(metaScf2,ln.refScf2,(ln.edgeScf2=="start")?ln.refStart2:ln.refEnd2,ln,ln.edgeScf2);

//    std::cerr<<"End find edges"<<std::endl;

    char orientMeta1=(edgeMeta1=="start")?metassem[metaScf1].begin()->orientation:(metassem[metaScf1].end()-1)->orientation;
    char orientMeta2=(edgeMeta2=="start")?metassem[metaScf2].begin()->orientation:(metassem[metaScf2].end()-1)->orientation;

//    std::cerr<<"Infor for scf1"<<std::endl;
    //Information for scf1
    if(ln.edgeScf1=="end"){
      startScf1copy=ln.refEnd1; //copy from startScf1copy to start scf1
      if(ln.refEnd1 < mum.coords[ln.aln1].scfLen1){
        refScf1Over=true;
        refScf1OverStart=ln.refEnd1;
        refScf1OverEnd=mum.coords[ln.aln1].scfLen1;
        refScf1OverSupp=Rce.zsignal(ln.refScf1,ln.refEnd1,'z');
        refScf1OverStart++;
      }
    } else if(ln.edgeScf1 == "start"){
      startScf1copy=ln.refStart1; //copy from startScf1copy to end scf1      
      if(ln.refStart1 > 1){
        refScf1Over=true;
        refScf1OverStart=ln.refStart1;
        refScf1OverEnd=1;
        refScf1OverSupp=Rce.zsignal(ln.refScf1,ln.refStart1,'z');
        refScf1OverStart--;
      }
    }

//    std::cerr<<"Infor for scf2"<<std::endl;
    //information for scf2
    if(ln.edgeScf2=="start"){
      startScf2copy=ln.refStart2;
      if(ln.refStart2 > 1){
        refScf2Over=true;
        refScf2OverStart=ln.refStart2;
        refScf2OverEnd=1;
        refScf2OverSupp=Rce.zsignal(ln.refScf2,ln.refStart2,'z');
        refScf2OverStart--;
      }
    } else if(ln.edgeScf2=="end"){
      startScf2copy=ln.refEnd2;
      if(ln.refEnd2 < mum.coords[ln.aln2].scfLen1){
        refScf2Over=true;
        refScf2OverStart=ln.refEnd2;
        refScf2OverEnd=mum.coords[ln.aln2].scfLen1;
        refScf2OverSupp=Rce.zsignal(ln.refScf2,ln.refEnd2,'z');
        refScf2OverStart++;
      }
    }

//    std::cerr<<"Infor for query scf"<<std::endl;
    //information for query Scf
    if(std::min(std::max(ln.queStart1,ln.queEnd1),std::max(ln.queStart2,ln.queEnd2))-std::max(std::min(ln.queStart1,ln.queEnd1),std::min(ln.queStart2,ln.queEnd2)) < -1){
      queInter=true;
      //queInterStart y queInterEnd estan en direccion 1-N
      queInterStart=std::min(std::max(ln.queStart1,ln.queEnd1),std::max(ln.queStart2,ln.queEnd2));
      queInterEnd=std::max(std::min(ln.queStart1,ln.queEnd1),std::min(ln.queStart2,ln.queEnd2));
      queInterStartSupp=Qce.zsignal(ln.queScf,queInterStart,'z');
      queInterEndSupp=Qce.zsignal(ln.queScf,queInterEnd,'z');
      queInterStart++;
      queInterEnd--;
      if((double)(queInterEnd - queInterStart) > Qce.mu0 - Qce.sd0){
        //If query insertion is too big add it and disregard ref overhangs.
        //Note that for all the bases in [queInterStart, queInterEnd] the coverage is >= minLinkCov
        //as determined by findScfLinks
        queInterStartSupp='n';
        queInterEndSupp='n';
      }
    } else {
      queInterEndSupp=Qce.zsignal(ln.queScf,
                              (std::min(std::max(ln.queStart1,ln.queEnd1),std::max(ln.queStart2,ln.queEnd2))+std::max(std::min(ln.queStart1,ln.queEnd1),std::min(ln.queStart2,ln.queEnd2)))/2,
                              'z');
      queInterStartSupp=queInterEndSupp;//zsignal of the middle of the overlap
    }

    std::vector<MetaRow>::iterator meta1It,meta2It;
    std::string metaScf1=Scf_Meta[ln.refScf1];
    std::string metaScf2=Scf_Meta[ln.refScf2];

//    std::cerr<<"Edge and desc="<<std::endl;

    if(edgeMeta2=="start"){
      meta2It=metassem[metaScf2].begin();
    }else if(edgeMeta2=="end"){
      meta2It=metassem[metaScf2].end()-1;
    }
    addLNcoord(&(*meta2It),ln,ln.edgeScf2,0);
    meta2It->mlnc->desc=std::string("Link. Reference scaffold: ").append(ln.refScf2).append(" breakpoint");
    if(edgeMeta1=="start"){
      meta1It=metassem[metaScf1].begin();
    } else if(edgeMeta1=="end"){
      meta1It=metassem[metaScf1].end()-1;
    }    
    addLNcoord(&(*meta1It),ln,ln.edgeScf1,0);
    meta1It->mlnc->desc=std::string("Link. Reference scaffold: ").append(ln.refScf1).append(" breakpoint");    

//    std::cerr<<"Decide and make links"<<std::endl;
    //Decide and make links (join the corresponding metassem scfs)
    if((refScf1Over && abs(refScf1OverStart-refScf1OverEnd)+1 < SMALL_SEQ &&
       refScf2Over && abs(refScf2OverStart-refScf2OverEnd)+1 < SMALL_SEQ &&
       !queInter) || (!queInter && queInterStartSupp == 'n')
    ){
      //join metassem, disregard small seq overhangs.
      int eovp=std::min(std::max(ln.queStart1,ln.queEnd1),std::max(ln.queStart2,ln.queEnd2));
      int sovp=std::max(std::min(ln.queStart1,ln.queEnd1),std::min(ln.queStart2,ln.queEnd2));
      int overlapAlns=eovp-sovp+1;
      if(overlapAlns>0){
        if(std::min(ln.queStart1,ln.queEnd1) < sovp){
          startScf1copy=mum.delta.query2ref(ln.queScf,sovp-1,ln.refScf1,ln.refStart1,ln.refEnd1);
          if(!startScf1copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<sovp-1<<" in "
                     <<ln.refScf1<<": "<<ln.refStart1<<" - "<<ln.refEnd1<<std::endl;
            exit(EXIT_FAILURE);
          }
          startScf2copy=mum.delta.query2ref(ln.queScf,eovp+1,ln.refScf2,ln.refStart2,ln.refEnd2);
          if(!startScf2copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<eovp+1<<" in "
                     <<ln.refScf2<<": "<<ln.refStart2<<" - "<<ln.refEnd2<<std::endl;
            exit(EXIT_FAILURE);
          }
        } else {
          startScf1copy=mum.delta.query2ref(ln.queScf,eovp+1,ln.refScf1,ln.refStart1,ln.refEnd1);
          if(!startScf1copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<eovp+1<<" in "
                     <<ln.refScf1<<": "<<ln.refStart1<<" - "<<ln.refEnd1<<std::endl;
            exit(EXIT_FAILURE);
          }
          startScf2copy=mum.delta.query2ref(ln.queScf,sovp-1,ln.refScf2,ln.refStart2,ln.refEnd2);
          if(!startScf2copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<sovp-1<<" in "
                     <<ln.refScf2<<": "<<ln.refStart2<<" - "<<ln.refEnd2<<std::endl;
            exit(EXIT_FAILURE);
          } 
        }
      }
      //now join
      meta1It->mlnc->desc=std::string("No query insertion and no overhangs. Reference scaffold: ").append(ln.refScf1).append(" breakpoint");
      meta2It->mlnc->desc=std::string("No query insertion and no overhangs. Reference scaffold: ").append(ln.refScf2).append(" breakpoint");
      joinMetassem(ln,edgeMeta1,startScf1copy,edgeMeta2,startScf2copy,false);
      
    }else if((!queInter && queInterStartSupp=='-')){
      if(refScf1Over && abs(refScf1OverStart-refScf1OverEnd)+1 >= SMALL_SEQ){
        //add Scf1 overhang
        MetaRow* newrow;
        newrow=addMetaRow(metaScf1,
                     edgeMeta1,
                     0,
                     ln.refScf1,
                     std::min(refScf1OverStart,refScf1OverEnd),
                     std::max(refScf1OverStart,refScf1OverEnd),
                     orientMeta1,
                     refScf1OverSupp,
                     'l',
                     queInterStartSupp,
                     queInterEndSupp,
                     meta1It->aln);
        (ln.edgeScf1=="start")?startScf1copy=std::min(refScf1OverStart,refScf1OverEnd):startScf1copy=std::max(refScf1OverStart,refScf1OverEnd);
        addLNcoord(newrow,ln,ln.edgeScf1,0);
        newrow->mlnc->desc=std::string("Reference scaffold: ").append(ln.refScf1).append(" overhang"); 
      }

      if(edgeMeta1=="start"){
        meta1It=metassem[metaScf1].begin();
      } else if(edgeMeta1=="end"){
        meta1It=metassem[metaScf1].end()-1;
      }
      //add Ns
      addMetaRow(metaScf1,
                    edgeMeta1,
                    0,
                    ln.refScf1,
                    1,
                    30,
                    'N',
                    refScf1OverSupp,
                    refScf2OverSupp,
                    queInterStartSupp,
                    queInterEndSupp,
                    meta1It->aln);
      if(edgeMeta2=="start"){
        meta2It=metassem[metaScf2].begin();
      } else if(edgeMeta2=="end"){
        meta2It=metassem[metaScf2].end()-1;
      }
      if(refScf2Over && abs(refScf2OverStart-refScf2OverEnd)+1 >= SMALL_SEQ){
        //add Scf2 overhang
        MetaRow* newrow;
        newrow=addMetaRow(metaScf2,
                     edgeMeta2,
                     0,
                     ln.refScf2,
                     std::min(refScf2OverStart,refScf2OverEnd),
                     std::max(refScf2OverStart,refScf2OverEnd), 
                     orientMeta2,
                     refScf2OverSupp,
                     'l',
                     queInterStartSupp,
                     queInterEndSupp,
                     meta2It->aln);
        (ln.edgeScf2=="start")?startScf2copy=std::min(refScf2OverStart,refScf2OverEnd):startScf2copy=std::max(refScf2OverStart,refScf2OverEnd);
        addLNcoord(newrow,ln,ln.edgeScf2,0);
        newrow->mlnc->desc=std::string("Reference scaffold: ").append(ln.refScf2).append(" overhang");
      }
      //join metassem
      joinMetassem(ln,edgeMeta1,startScf1copy,edgeMeta2,startScf2copy,false);
    } else if(queInter){
      //set the orientation in which queInter must be added.
      char orientQueInter;
      if(orientMeta1=='+'){
        orientQueInter=ln.queOrient1;  
      } else if(orientMeta1=='-'){
        if(ln.queOrient1=='+'){
          orientQueInter='-';
        }else{
          orientQueInter='+';
        }
      }
      meta1It->mlnc->desc.append(std::string(". Query ").append(ln.queScf).append(" insertion start"));
      meta2It->mlnc->desc.append(std::string(". Query ").append(ln.queScf).append(" insertion end"));
      //insert query insertion
      addMetaRow(metaScf1,
                 edgeMeta1,
                 1,
                 ln.queScf,
                 queInterStart,
                 queInterEnd,
                 orientQueInter,
                 refScf1OverSupp,
                 refScf2OverSupp,
                 queInterStartSupp,
                 queInterEndSupp,
                 meta1It->aln);
      //join metassem
      joinMetassem(ln,edgeMeta1,startScf1copy,edgeMeta2,startScf2copy,true);
    } else {
      //join metassem.
      int eovp=std::min(std::max(ln.queStart1,ln.queEnd1),std::max(ln.queStart2,ln.queEnd2));
      int sovp=std::max(std::min(ln.queStart1,ln.queEnd1),std::min(ln.queStart2,ln.queEnd2));
      int overlapAlns=eovp-sovp+1;
      if(overlapAlns>0){
        if(std::min(ln.queStart1,ln.queEnd1) < sovp){
          startScf1copy=mum.delta.query2ref(ln.queScf,sovp-1,ln.refScf1,ln.refStart1,ln.refEnd1);
          if(!startScf1copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<sovp-1<<" in "
                     <<ln.refScf1<<": "<<ln.refStart1<<" - "<<ln.refEnd1<<std::endl;
            exit(EXIT_FAILURE);
          }
          startScf2copy=mum.delta.query2ref(ln.queScf,eovp+1,ln.refScf2,ln.refStart2,ln.refEnd2);
          if(!startScf2copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<eovp+1<<" in "
                     <<ln.refScf2<<": "<<ln.refStart2<<" - "<<ln.refEnd2<<std::endl;
            exit(EXIT_FAILURE);
          }
        } else {
          startScf1copy=mum.delta.query2ref(ln.queScf,eovp+1,ln.refScf1,ln.refStart1,ln.refEnd1);
          if(!startScf1copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<eovp+1<<" in "
                     <<ln.refScf1<<": "<<ln.refStart1<<" - "<<ln.refEnd1<<std::endl;
            exit(EXIT_FAILURE);
          }
          startScf2copy=mum.delta.query2ref(ln.queScf,sovp-1,ln.refScf2,ln.refStart2,ln.refEnd2);
          if(!startScf2copy){
            std::cerr<<"ERROR: Could not find "<<ln.queScf<<": "<<sovp-1<<" in "
                     <<ln.refScf2<<": "<<ln.refStart2<<" - "<<ln.refEnd2<<std::endl;
            exit(EXIT_FAILURE);
          } 
        }
      }
      //now join
      joinMetassem(ln,edgeMeta1,startScf1copy,edgeMeta2,startScf2copy,false);
    }
  }
}

//----------------------------------- addLNcoord -----------------------------//
void Metassembly::addLNcoord(MetaRow* meta, MetaLink imln, std::string iedge, int idelta)
{
  metaLNcoords* mlnc= new metaLNcoords(imln,meta->orientation,iedge,idelta);
  LNs.push_back(mlnc);
  meta->mlnc=mlnc;
}

//---------------------------------- addNAcoord ------------------------------//
void Metassembly::addNAcoord(MetaRow* meta, std::string irscf, int irnstart, int irnend,
                                            std::string iqscf, int iqostart, int iqoend,
                                            std::string iqoside, char iqor)
{
  metaNAcoords* mnac= new metaNAcoords(irscf,irnstart,irnend,iqscf,iqostart,iqoend,iqoside,iqor);
  NAs.push_back(mnac);
  meta->mNAc=mnac;
}

//----------------------------------- addNBcoord -----------------------------//
void Metassembly::addNBcoord(MetaRow* meta, std::string iqscf, int iqostart, int iqoend,
                                            std::string irscf, int irnstart, int irnend,
                                            int iristart, int iriend,
                                            int irbkp,
                                            std::string iqoedge, char iqor)
{
  metaNBcoords* mnbc= new metaNBcoords(iqscf,iqostart,iqoend,irscf,irnstart,irnend,iristart,iriend,irbkp,iqoedge,iqor);
  NBs.push_back(mnbc);
  meta->mNBc=mnbc;
}

//---------------------------------- addNCcoord ------------------------------//
void Metassembly::addNCcoord(MetaRow* meta, std::string iqscf1, int iqlostart, int iqloend, char iqlor,
                                             std::string irscf, int irlnstart, int irlnend,
                                             int iristart, int iriend,
                                             int irrnstart, int irrnend,
                                             std::string iqscf2, int iqrostart, int iqroend, char iqror)
{
  metaNCcoords* mncc= new metaNCcoords (iqscf1,iqlostart,iqloend,iqlor,
                     irscf,irlnstart,irlnend,
                     iristart,iriend,
                     irrnstart,irrnend,
                     iqscf2,iqrostart,iqroend,iqror);
  NCs.push_back(mncc);
  meta->mNCc=mncc;
}

//---------------------------------- addN1coord ------------------------------//
void Metassembly::addN1coord(MetaRow* meta, std::string iqscf, int iqistart, int iqiend, char iqor,
                                            std::string irscf, int iristart, int iriend,
                                            int irnstart, int irnend)
{
  metaN1coords* mn1c= new metaN1coords (iqscf,iqistart,iqiend,iqor,irscf,iristart,iriend,irnstart,irnend);
  N1s.push_back(mn1c);
  meta->mN1c=mn1c;
}

//----------------------------------- addQIcoord -----------------------------//
void Metassembly::addQIcoord(MetaRow* meta, std::string iqscf, int iqistart, int iqiend, char iqor,
                                            std::string irscf, int iristart, int iriend)
{
  metaQIcoords* mqic=new metaQIcoords(iqscf,iqistart,iqiend,iqor,irscf,iristart,iriend);
  QIs.push_back(mqic);
  meta->mQIc=mqic;
}

//----------------------------------- addBKPcoord ----------------------------//
void Metassembly::addBKPcoord(MetaRow* meta, std::string irscf, int irbkp, 
                                             std::string iqscf, int iqostart, int iqoend,
                                             std::string iqoside, char iqor)
{
  metaBKPcoords* bkp= new metaBKPcoords (irscf,irbkp,iqscf,iqostart,iqoend,iqoside,iqor);
  BKPs.push_back(bkp);
  meta->mBKPc=bkp;
}

//----------------------------------- joinMetassem ---------------------------//
void Metassembly::joinMetassem(MetaLink& ln, std::string edgeMeta1, int startScf1copy, std::string edgeMeta2, int startScf2copy, bool addedQueInter)
{

//  std::cerr<<"\n------------------------Join Meta:-------------------"<<std::endl;
//  ln.print(std::cerr);
  std::string metaScf1=Scf_Meta[ln.refScf1];
  std::string metaScf2=Scf_Meta[ln.refScf2];


/*  std::vector<MetaRow>::iterator coords;
  std::cerr << "\n\nOld metassem: "<<metaScf1 <<":"<<std::endl;
  for(coords=metassem[metaScf1].begin(); coords!=metassem[metaScf1].end(); coords++){
    std::cerr << metaScf1 << "\t"
           << (*coords).file << "\t"
           << (*coords).scf << "\t"
           << (*coords).start << "\t"
           << (*coords).end << "\t"
           << (*coords).orientation << "\t"
           << (*coords).sigRef1 << "\t"
           << (*coords).sigRef2 << "\t"
           << (*coords).sigQue1 << "\t"
           << (*coords).sigQue2 << "\t"
           << mum.coords[(*coords).aln].scf1 << "\t"
           << mum.coords[(*coords).aln].start1 << "\t"
           << mum.coords[(*coords).aln].end1 << "\t"
           << mum.coords[(*coords).aln].scf2 << "\t"
           << mum.coords[(*coords).aln].start2 << "\t"
           << mum.coords[(*coords).aln].end2
           << std::endl; 
  }  
  std::cerr << "\n\nOld metassem: "<<metaScf2 <<":"<<std::endl;
  for(coords=metassem[metaScf2].begin(); coords!=metassem[metaScf2].end(); coords++){
      std::cerr << metaScf2 << "\t"
           << (*coords).file << "\t"
           << (*coords).scf << "\t"
           << (*coords).start << "\t"
           << (*coords).end << "\t"
           << (*coords).orientation << "\t"
           << (*coords).sigRef1 << "\t"
           << (*coords).sigRef2 << "\t"
           << (*coords).sigQue1 << "\t"
           << (*coords).sigQue2 << "\t"
           << mum.coords[(*coords).aln].scf1 << "\t"
           << mum.coords[(*coords).aln].start1 << "\t"
           << mum.coords[(*coords).aln].end1 << "\t"
           << mum.coords[(*coords).aln].scf2 << "\t"
           << mum.coords[(*coords).aln].start2 << "\t"
           << mum.coords[(*coords).aln].end2
           << std::endl; 
  }*/


  int addedSeqCount=0; 
 
  if(metassem[metaScf1].size()<=0 ){ //just precaution
    std::cerr<<"Error: "<<metaScf1<<" has size<=0"<<std::endl;
    ln.print(std::cerr);
    exit(EXIT_FAILURE);
  } else if(metassem[metaScf2].size()<=0){
    std::cerr<<"Error: "<<metaScf2<<" has size<=0"<<std::endl;
    ln.print(std::cerr);
    exit(EXIT_FAILURE);
  }
     
  std::vector<MetaRow>::iterator meta1It,meta2It;
  std::vector<MetaRow>::iterator finalMeta2It;
  unsigned int skip=0;

  int eovp=std::min(std::max(ln.queStart1,ln.queEnd1),std::max(ln.queStart2,ln.queEnd2));
  int sovp=std::max(std::min(ln.queStart1,ln.queEnd1),std::min(ln.queStart2,ln.queEnd2));
  int overlapAlns=eovp-sovp+1;
  if(overlapAlns > 0 && addedQueInter){

    if(edgeMeta1 == "start"){
      for(meta1It=metassem[metaScf1].begin();
          meta1It < metassem[metaScf1].end();
          meta1It++){
        skip++;
        if(meta1It -> scf == ln.refScf1 &&
           (meta1It -> start <= startScf1copy &&
            meta1It -> end >= startScf1copy)
          )
        {
          break;
        }
      }
      if(meta1It == metassem[metaScf1].end()){
        std::cerr<<"Error joining links. Could not find mateAnRow for "<<std::endl
                 <<"startScf1copy: "<<startScf1copy<<std::endl;
      }
      if(meta1It != metassem[metaScf1].begin()){
        meta1It -> mlnc = metassem[metaScf1].begin()->mlnc;
      }
      if(skip > 1){
        std::vector<MetaRow>::iterator it;
        for(it=metassem[metaScf1].begin();
            it<meta1It;
            it++){
          addedSeqCount-=it->end - it->start +1;
        }
        metassem[metaScf1].erase(metassem[metaScf1].begin(),meta1It);
        meta1It=metassem[metaScf1].begin();
      }
    }else if(edgeMeta1 == "end"){
      for(meta1It=metassem[metaScf1].end() -1;
          meta1It >= metassem[metaScf1].begin();
          meta1It--){
        skip++;
        if(meta1It -> scf == ln.refScf1 &&
           meta1It -> start <= startScf1copy &&
           meta1It -> end >= startScf1copy){
          break;
        }
      }
      if(meta1It == metassem[metaScf1].begin() &&
         !(meta1It -> scf == ln.refScf1 &&
           meta1It -> start <= startScf1copy &&
           meta1It -> end >= startScf1copy
         )
        ){
        std::cerr<<"Error joining links. Could not find mateAnRow for"<<std::endl
                 <<"startScf1copy: "<<startScf1copy<<std::endl;
      }
      if(meta1It != metassem[metaScf1].end() -1){
         meta1It -> mlnc = (metassem[metaScf1].end()-1)->mlnc;
      }
      if(skip > 1){
        std::vector<MetaRow>::iterator it;
        for(it = metassem[metaScf1].end() -1;
            it > meta1It;
            it --){
          addedSeqCount -= it->end - it->start +1;
        }
        metassem[metaScf1].erase(meta1It+1,metassem[metaScf1].end());
        meta1It=metassem[metaScf1].end()-1;
      }
    }

    if(ln.edgeScf1 == "start"){
      addedSeqCount-=std::max(meta1It->start,startScf1copy)-std::min(meta1It->start,startScf1copy)+1;
      meta1It->start=startScf1copy;
    }else{
      addedSeqCount-=std::max(meta1It->end,startScf1copy)-std::min(meta1It->end,startScf1copy)+1;
      meta1It->end=startScf1copy;
    }

    addMetaRow(metaScf1,edgeMeta1,1,ln.queScf,sovp,eovp,(ln.queStart1<=ln.queEnd1)?'+':'-',
               Qce.zsignal(ln.queScf,sovp,'z'),
               Qce.zsignal(ln.queScf,eovp,'z'),
               'l',
               'l',
               ln.aln1); 
    addedSeqCount+=eovp-sovp+1;  
  }

  skip=0;
  if(edgeMeta2=="start"){
    for(meta2It=metassem[metaScf2].begin();
        meta2It < metassem[metaScf2].end();
        meta2It++){
      skip++;
      if(meta2It -> scf == ln.refScf2 &&
         (meta2It -> start <= startScf2copy &&
          meta2It -> end >= startScf2copy)
        )
      {
        break;
      }
    }
    if(meta2It == metassem[metaScf2].end()){
      std::cerr<<"Error joining links. Could not find mateAnRow for "<<std::endl
               <<"startScf2copy: "<<startScf2copy<<std::endl;
    }
    if(meta2It != metassem[metaScf2].begin()){
      meta2It -> mlnc = metassem[metaScf2].begin()->mlnc;
    }
  }else if(edgeMeta2=="end"){
    for(meta2It=metassem[metaScf2].end() -1;
        meta2It >= metassem[metaScf2].begin();
        meta2It--){
      skip++;
      if(meta2It -> scf == ln.refScf2 &&
         meta2It -> start <= startScf2copy &&
         meta2It -> end >= startScf2copy){
        break;
      }
    }
    if(meta2It == metassem[metaScf2].begin() &&
       !(meta2It -> scf == ln.refScf2 &&
         meta2It -> start <= startScf2copy &&
         meta2It -> end >= startScf2copy
       )
      ){
      std::cerr<<"Error joining links. Could not find mateAnRow for"<<std::endl
               <<"startScf2copy: "<<startScf2copy<<std::endl;
    }
    if(meta2It != metassem[metaScf2].end() -1){
       meta2It -> mlnc = (metassem[metaScf2].end()-1)->mlnc;
    }
  }
  if(edgeMeta1=="start"){
    meta1It=metassem[metaScf1].begin(); 
  } else if(edgeMeta1=="end"){
    meta1It=metassem[metaScf1].end();
  }
        
  //change metassem name for *meta2It
  //if *meta2It.scf is a reference scaffold,
  //that is, if *meta2It = ln.refScf2
  if(Scf_Meta.count(meta2It->scf)>0){ 
    Scf_Meta[meta2It->scf]=metaScf1;
  }
  //add first assemCoords to modify it
  metassem[metaScf1].insert(meta1It,*meta2It); 
      
  //change startScf2copy if necessary
  if(edgeMeta1=="start"){
    meta1It=metassem[metaScf1].begin(); 
  } else if(edgeMeta1=="end"){
   meta1It=metassem[metaScf1].end()-1;
  }
      
  //only change start|end if there is no querInter
  if(!addedQueInter){
    if(ln.edgeScf2 == "start"){
      meta1It->start=startScf2copy;
    } else if(ln.edgeScf2 == "end"){
      meta1It->end=startScf2copy;
    } 
  }
  addedSeqCount+=labs(meta1It->end - meta1It->start)+1;


  //Add metaScf2 to metaScf1
  if(edgeMeta1=="start"){
    meta1It=metassem[metaScf1].begin(); 
  } else if(edgeMeta1=="end"){
    meta1It=metassem[metaScf1].end();
  }

  if(edgeMeta1 != edgeMeta2){ //No need to reverse metaScf2 sequence, just add it to metaScf1
    std::vector<MetaRow>::iterator changeMetaIt;
    if(edgeMeta2 == "start"){
      meta2It=metassem[metaScf2].begin()+skip;
      finalMeta2It=metassem[metaScf2].end();
      if(metassem[metaScf2].size()>1){
        //Change metassem name in Scf_Meta for all ref scaffolds
        //present in metaScf2
        for(changeMetaIt=meta2It;changeMetaIt <= finalMeta2It-1; changeMetaIt++){
          
          addedSeqCount+=labs(changeMetaIt->end - changeMetaIt->start)+1;
          
          if(Scf_Meta.count(changeMetaIt->scf)>0){
            Scf_Meta[changeMetaIt->scf]=metaScf1;
          }
        }
        //Insert metaScf2 in metaScf1
        metassem[metaScf1].insert(meta1It,meta2It,finalMeta2It);
      }
    }else if(edgeMeta2 == "end"){
      meta2It=metassem[metaScf2].begin();
      finalMeta2It=metassem[metaScf2].end()-skip;
      if(metassem[metaScf2].size()>1){
        //Change metassem name in Scf_Meta for all ref scaffolds
        //present in metaScf2
        for(changeMetaIt=meta2It;changeMetaIt <= finalMeta2It-1; changeMetaIt++){
           addedSeqCount+=labs(changeMetaIt->end - changeMetaIt->start)+1;
           if(Scf_Meta.count(changeMetaIt->scf)>0){
            Scf_Meta[changeMetaIt->scf]=metaScf1;
          }
        }
        //Insert metaScf2 in metaScf1
        metassem[metaScf1].insert(meta1It,meta2It,finalMeta2It);
      }
    }  
  } else { //Need to reverse metaScf2 sequence and then add it to metaScf1
    std::vector<MetaRow>::reverse_iterator startMeta2;
    std::vector<MetaRow>::reverse_iterator endMeta2;
    std::vector<MetaRow>::reverse_iterator rit;
    (edgeMeta1=="start")?meta1It->changeOrientation():
                         (meta1It-1)->changeOrientation();
    if(edgeMeta2=="start"){
      startMeta2=metassem[metaScf2].rbegin();
      endMeta2=metassem[metaScf2].rend()-skip;
      //Change orientation in metaScf2 and change metassem name in Scf_Meta
      for(rit=startMeta2; rit<endMeta2; rit++){
        rit->changeOrientation();
        addedSeqCount+=labs(rit->end - rit->start)+1;
        if(Scf_Meta.count(rit->scf)>0){
          Scf_Meta[rit->scf]=metaScf1;
        }
      }
      //Add to metaScf1
      if(metassem[metaScf2].size()>1){
        metassem[metaScf1].insert(meta1It,startMeta2,endMeta2);  
      }
    } else if(edgeMeta2=="end"){
      startMeta2=metassem[metaScf2].rbegin()+skip;
      endMeta2=metassem[metaScf2].rend();
      //Change orietnation in metaScf2 and change metassem name in Scf_Meta
      for(rit=startMeta2;rit<endMeta2;rit++){
        rit->changeOrientation();
        addedSeqCount+=labs(rit->end - rit->start)+1;
        if(Scf_Meta.count(rit->scf)>0){
          Scf_Meta[rit->scf]=metaScf1;
        }
      }
      //Add to metaScf1
      if(metassem[metaScf2].size()>1){
        metassem[metaScf1].insert(meta1It,startMeta2,endMeta2);
      }
    }
  }

//    Scf_Meta[ln.refScf2]=metaScf1; //now refScf2 is in metaScf1
  metaLengths[metaScf1]+=addedSeqCount;
  std::map<std::string, int>::iterator itLen;
  itLen=metaLengths.find(metaScf2);
  metaLengths.erase(itLen);
    
  std::map<std::string, std::vector<MetaRow> >::iterator it;
  it=metassem.find(metaScf2);
  metassem.erase(it); //erase metaScf2 from metassem
  
/*  std::cerr << "\n\nNew metassem: "<<metaScf1 <<":"<<std::endl;
  for(coords=metassem[metaScf1].begin(); coords!=metassem[metaScf1].end(); coords++){
    std::cerr << metaScf1 << "\t"
           << (*coords).file << "\t"
           << (*coords).scf << "\t"
           << (*coords).start << "\t"
           << (*coords).end << "\t"
           << (*coords).orientation << "\t"
           << (*coords).sigRef1 << "\t"
           << (*coords).sigRef2 << "\t"
           << (*coords).sigQue1 << "\t"
           << (*coords).sigQue2 << "\t"
           << mum.coords[(*coords).aln].scf1 << "\t"
           << mum.coords[(*coords).aln].start1 << "\t"
           << mum.coords[(*coords).aln].end1 << "\t"
           << mum.coords[(*coords).aln].scf2 << "\t"
           << mum.coords[(*coords).aln].start2 << "\t"
           << mum.coords[(*coords).aln].end2
           << std::endl; 
    }  */

  return;
}

//----------------------------------- getRefOrientation ----------------------//
char Metassembly::getRefOrientation(std::string scf, int start, int end)
{
  if(inversions[scf][start].count(end) > 0){
    return(inversions[scf][start][end]);
  } else {
    return('+');
  }
}

//----------------------------------- lastPosAdded ---------------------------//
int Metassembly::lastPosAdded2Meta(std::string scf, std::string metaScf)
{
 std::vector<MetaRow>::reverse_iterator rit;
 if(metassem.count(metaScf) ==0){
   std::cerr<< "Error: lastPosAdded:\n"
       << "no metaScf named "<<metaScf<<" found\n";
   exit(EXIT_FAILURE);
 }
 for(rit=metassem[metaScf].rbegin(); rit<metassem[metaScf].rend(); rit++){
   if(rit->scf == scf){
     return(rit->end);
   }    
 }
 return(0); 
}

//----------------------------------- addLink --------------------------------//
void Metassembly::addLink(int alnScf1, std::string edgeScf1, int alnScf2, std::string edgeScf2)
{
  MetaLink ln;
  
  ln.refScf1=mum.coords[alnScf1].scf1;
  ln.refStart1=mum.coords[alnScf1].start1;
  ln.refEnd1=mum.coords[alnScf1].end1;
  ln.edgeScf1=edgeScf1;
  ln.aln1=alnScf1;
  
  ln.queScf=mum.coords[alnScf1].scf2;
  ln.queStart1=mum.coords[alnScf1].start2;
  ln.queEnd1=mum.coords[alnScf1].end2;
  ln.queOrient1=(ln.queStart1 <= ln.queEnd1)? '+' :'-';
  
  ln.queStart2=mum.coords[alnScf2].start2;
  ln.queEnd2=mum.coords[alnScf2].end2;
  ln.queOrient2=(ln.queStart2 <= ln.queEnd2)?'+':'-';
  
  ln.refScf2=mum.coords[alnScf2].scf1;
  ln.refStart2=mum.coords[alnScf2].start1;
  ln.refEnd2=mum.coords[alnScf2].end1;
  ln.edgeScf2=edgeScf2;
  ln.aln2=alnScf2;
  
  std::string nameLink(ln.refScf1);
  (edgeScf1 == "start")?nameLink.append("1"):nameLink.append("2");
  scfsLinks[nameLink]=ln;
}

//----------------------------------- findEdgeInMeta -------------------------//
std::string Metassembly::findEdgeInMeta(std::string metaScf,
                           std::string refScf,
                           int position,
                           MetaLink ln,
                           std::string prefEdge){

  bool start=false;
  bool end=false;
  
  if(metassem[metaScf].begin()->scf==refScf &&
    metassem[metaScf].begin()->start <=position &&
    metassem[metaScf].begin()->end >= position
    ){
    start=true;    
  } 
  if((metassem[metaScf].end()-1)->scf==refScf &&
            (metassem[metaScf].end()-1)->start <=position &&
            (metassem[metaScf].end()-1)->end >= position
    ){
    end=true;          
  }
  
  if(start && end){
	  return(prefEdge);
  } else if(start){
    return("start");
  } else if(end){
    return("end");
  } else {
    std::cerr<<"\nError linking, no edge of metassem of\n"
        <<refScf<<" in "<<metaScf<<":"<<std::endl;
    ln.print(std::cerr);
    mum.coords[ln.aln1].print(std::cerr);
    mum.coords[ln.aln2].print(std::cerr);
    std::cerr<<"\n"
        <<"Start copy is: "<<position<<"\n"<<"\n"
        <<metaScf<<" front:"<<metassem[metaScf].begin()->scf <<"\n"
        <<metaScf<<" start:"<<metassem[metaScf].begin()->start <<"\n"
        <<metaScf<<" end  :"<<metassem[metaScf].begin()->end <<"\n"
        <<"\n"<<"\n"
        <<metaScf<<" back :"<<(metassem[metaScf].end()-1)->scf <<"\n"
        <<metaScf<<" start:"<<(metassem[metaScf].end()-1)->start <<"\n"
        <<metaScf<<" end  :"<<(metassem[metaScf].end()-1)->end <<"\n\n";
    exit(EXIT_FAILURE);
  }   
}

//----------------------------------- analyseAln -----------------------------//
std::vector<int>::iterator Metassembly::analyseAln(std::string refScf,
                                        std::vector<int>& rScfAlns,
                                        std::vector<int>::iterator curScfAln,
                                        std::string metaScf,
                                        bool leftLink,
                                        std::vector<int>::iterator leftAln,
                                        bool rightLink,
                                        std::vector<int>::iterator rightAln
                                       )
{

 if(curScfAln > rScfAlns.begin() && 
    mum.coords[*curScfAln].start1 >= mum.coords[*(curScfAln-1)].start1 && 
    mum.coords[*curScfAln].end1 <= mum.coords[*(curScfAln-1)].end1)
    return(curScfAln+1);

  bool overhangRef=false;
  char overhangRefSupport='o';
  bool overhangQuery=false;
  char overhangQuerySupport='o';
  char orientQuery=(mum.coords[*curScfAln].start2 <= mum.coords[*curScfAln].end2)?'+':'-';
  char orientRef=getRefOrientation(refScf,mum.coords[*curScfAln].start1,mum.coords[*curScfAln].end1);
  std::vector<int> qOverhangPos(2,0);
  
  //Add left overhang sequences to metassembly
  std::vector<int>& queAlns=mum.queryScfs[mum.coords[*curScfAln].scf2];
  
  if(orientQuery == '+'){                                                                      
    if(mum.coords[*curScfAln].start2 > 1+19 && 
       !mum.QscfBefore(mum.coords[*curScfAln].scf2,rScfAlns,curScfAln) && 
       mum.adjacentAlnQue(queAlns,curScfAln,'-') == queAlns.end()
      )
    {
      qOverhangPos[0]=1;
      qOverhangPos[1]=mum.coords[*curScfAln].start2-1;
      int QZstati=Qce.search(mum.coords[*curScfAln].scf2,qOverhangPos[1]+1);
      mateAnRow QZstat;
      if(QZstati >= 0)
        QZstat=Qce[mum.coords[*curScfAln].scf2][QZstati];
      bool Qcov;
      Qcov=(QZstati >= 0 && QZstat.N >= minCov)?true:false;
      overhangQuery= ! QNs.isNinsert(mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],QmaxNinsert,QmaxNinsPer);
      overhangQuery= overhangQuery && Qcov;
      overhangQuerySupport=Qce.zsignal(mum.coords[*curScfAln].scf2,qOverhangPos[1]+1,'z');
    }                                                                                                                         
  } else if(mum.coords[*curScfAln].start2 < mum.coords[*curScfAln].scfLen2-19 && 
            !mum.QscfBefore(mum.coords[*curScfAln].scf2,rScfAlns,curScfAln) && 
            mum.adjacentAlnQue(queAlns,curScfAln,'+') == queAlns.end()
           )
  {
    qOverhangPos[0]=mum.coords[*curScfAln].start2+1;
    qOverhangPos[1]=mum.coords[*curScfAln].scfLen2;
    int QZstati=Qce.search(mum.coords[*curScfAln].scf2,qOverhangPos[0]-1);
    mateAnRow QZstat;
    if(QZstati >= 0)
      QZstat=Qce[mum.coords[*curScfAln].scf2][QZstati];
    bool Qcov;
    Qcov=(QZstati >= 0 && QZstat.N >= minCov)?true:false;
    overhangQuery= !QNs.isNinsert(mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],QmaxNinsert,QmaxNinsPer);
    overhangQuery=overhangQuery && Qcov;
    overhangQuerySupport=Qce.zsignal(mum.coords[*curScfAln].scf2,qOverhangPos[0]-1,'z');    
  }  
  int QovLen=std::max(qOverhangPos[0],qOverhangPos[1])-std::min(qOverhangPos[0],qOverhangPos[1])+1;
  overhangQuery=overhangQuery && QovLen <= maxQinsLength;
 
  if(curScfAln == rScfAlns.begin() && !leftLink){
    if(mum.coords[*curScfAln].start1 > 1+19){
      overhangRef=true;
      overhangRefSupport=Rce.zsignal(refScf,mum.coords[*curScfAln].start1,'z');
    } else if(mum.coords[*curScfAln].start1 > 1){
      //Add little overhang to metassembly
      addMetaRow(metaScf,"end",0,refScf,1,mum.coords[*curScfAln].start1-1,'+','o','o','o','o',*curScfAln);
    }
     
    if(overhangRef){
      if(overhangQuery){
	bool oNinsert=RNs.isNinsert(refScf,std::max( 1,mum.coords[*curScfAln].start1-100),mum.coords[*curScfAln].start1-1,minNinsert,minNinsPer);
	int RNinsLength=0;
	int NinsIdx=-1;
	if(oNinsert){
		NinsIdx=RNs.search(refScf,std::max( 1,mum.coords[*curScfAln].start1-100),mum.coords[*curScfAln].start1-1);
		RNinsLength=RNs[refScf][NinsIdx][1]-RNs[refScf][NinsIdx][0]+1;
	}
        int RZstat1i=Rce.search(refScf,mum.coords[*curScfAln].start1);
        mateAnRow RZstat1;
        if(RZstat1i>=0)
          RZstat1=Rce[refScf][RZstat1i];
        
        int RZstat2i=-1;
	if(NinsIdx >= 0) 
          RZstat2i=Rce.search(refScf,RNs[refScf][NinsIdx][0]);
        mateAnRow RZstat2;
        if(RZstat2i>=0)
          RZstat2=Rce[refScf][RZstat2i];
        

        if(oNinsert && RNinsLength <= maxRinsLength && QovLen <= maxQinsLength &&
           (RZstat1i>=0 && fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RNinsLength,QovLen)) < minZstat) &&
           (RZstat2i>=0 && fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RNinsLength,QovLen)) < minZstat)
            )
         {
            if(RNs[refScf][NinsIdx][0] > 1){
              //insert ref overhang 
              addMetaRow(metaScf,"end",0,refScf,1,RNs[refScf][NinsIdx][0]-1,'+','o','N',overhangQuerySupport,'o',*curScfAln);
              countNrows++;
            }
            addMetaRow(metaScf,"end",0,"InsQueOver",1,30,'N','o','N',overhangQuerySupport,'o',*curScfAln);
            countNrows++;
            //insert query overhang
            MetaRow* newrow;
            newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,'o','N',overhangQuerySupport,'o',*curScfAln);
            addNAcoord(newrow,refScf, RNs[refScf][NinsIdx][0], RNs[refScf][NinsIdx][1],
                       mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                       "right", orientQuery                      
                      );
            countNevents++;
            countNA++;
         }else if( overhangQuerySupport=='n' && RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
	           fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QovLen)) >= minDiffZ &&
                   fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QovLen)) <= minZstat
          ){
          MetaRow* newrow;
          newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,'o',overhangRefSupport,overhangQuerySupport,'o',*curScfAln);
          addBKPcoord(newrow, refScf,mum.coords[*curScfAln].start1, 
                      mum.coords[*curScfAln].scf2,qOverhangPos[0], qOverhangPos[1], "right", orientQuery);
          countBKP++;
        } else {
          addMetaRow(metaScf,"end",0,refScf,1,mum.coords[*curScfAln].start1-1,'+','o',overhangRefSupport,overhangQuerySupport,'o',*curScfAln);
        }
      } else {
        addMetaRow(metaScf,"end",0,refScf,1,mum.coords[*curScfAln].start1-1,'+','o',overhangRefSupport,'o','o',*curScfAln);
      }
    } else if(overhangQuery && overhangQuerySupport=='n'){ 
      MetaRow* newrow;
      newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,'.','.',overhangQuerySupport,'o',*curScfAln);
      addBKPcoord(newrow,refScf,mum.coords[*curScfAln].start1, 
                  mum.coords[*curScfAln].scf2,qOverhangPos[0], qOverhangPos[1],"right",orientQuery);
      countBKP++;
    }
  }

 
  //---Add alignment sequence to metassembly
  if(curScfAln != rScfAlns.begin() && curScfAln > leftAln){
    addMetaRow(metaScf,"end",0,refScf,
                 std::max(mum.coords[*curScfAln].start1,lastPosAdded2Meta(refScf,metaScf)+1), /*avoid inserting the same sequence twice when two contiguous alignmentes overlap*/
                 mum.coords[*curScfAln].end1,
                 orientRef,
                 Rce.zsignal(refScf,mum.coords[*curScfAln].start1,'z'),
                 Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                 Qce.zsignal(mum.coords[*curScfAln].scf2,mum.coords[*curScfAln].start2,'z'),
                 Qce.zsignal(mum.coords[*curScfAln].scf2,mum.coords[*curScfAln].end2,'z'),
                 *curScfAln);
    if(orientRef == '-'){
      (metassem[metaScf].end()-1)->isInv=true;
    }        
  } else {
    addMetaRow(metaScf,"end",0,refScf,
                  mum.coords[*curScfAln].start1,
                  mum.coords[*curScfAln].end1,
                  orientRef,
                  Rce.zsignal(refScf,mum.coords[*curScfAln].start1,'z'),
                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                  Qce.zsignal(mum.coords[*curScfAln].scf2,mum.coords[*curScfAln].start2,'z'),
                  Qce.zsignal(mum.coords[*curScfAln].scf2,mum.coords[*curScfAln].end2,'z'),
                  *curScfAln);
    if(orientRef == '-'){
      (metassem[metaScf].end()-1)->isInv=true;
    }
  }


  std::vector<int>::iterator qScfReturnPos;
  qScfReturnPos=mum.QscfAfter(mum.coords[*curScfAln].scf2,rScfAlns,curScfAln); //whether qScf aligns in a subsequent position within refScf

  if(qScfReturnPos!= rScfAlns.end() && 
     qScfReturnPos <= rightAln  &&
     mum.analyseInterInsert(rScfAlns,curScfAln,qScfReturnPos,Rce)){
    char orientReturnQuery;
    orientReturnQuery=(mum.coords[*qScfReturnPos].start2 <= mum.coords[*qScfReturnPos].end2)?'+':'-';

    if(orientReturnQuery != orientQuery && orientRef=='+'){
	  //Inversions
    
    //  ----|----->-----<-----|------ Query
    //      Q1   Q2     Q4    Q3
    
    //  ----|-----|-----|-----|------ Reference
    //      R1    R2   R3     R4
      
      int Q1=mum.coords[*curScfAln].start2;
      int Q2=mum.coords[*curScfAln].end2;
      int Q3=mum.coords[*qScfReturnPos].start2;
      int Q4=mum.coords[*qScfReturnPos].end2;
      
      int R1=mum.coords[*curScfAln].start1;
      int R2=mum.coords[*curScfAln].end1;
      int R3=mum.coords[*qScfReturnPos].start1;
      int R4=mum.coords[*qScfReturnPos].end1;
      
      //Process overlaps 
      if(orientQuery=='+' && Q4<Q2){
	      Q4=Q2+1;
      } else if (orientQuery=='-' && Q2<Q4){
		    Q2=Q4+1;
      }   
      if(R3 < R2){
        R3=mum.coords[*curScfAln].end1+1;	  
      }
	   
      char Q1sup=Qce.zsignal(mum.coords[*curScfAln].scf2,Q1,'i');
      char Q2sup=Qce.zsignal(mum.coords[*curScfAln].scf2,Q2,'i');
      char Q3sup=Qce.zsignal(mum.coords[*qScfReturnPos].scf2,Q3,'i');
      char Q4sup=Qce.zsignal(mum.coords[*qScfReturnPos].scf2,Q4,'i'); 
      char R1sup=Rce.zsignal(refScf,R1,'i');
      char R2sup=Rce.zsignal(refScf,R2,'i');
      char R3sup=Rce.zsignal(refScf,R3,'i');
      char R4sup=Rce.zsignal(refScf,R4,'i');

      if(!(R1sup=='n' && R2sup=='n' && R3sup=='n' && R4sup=='n')){
        if(R1sup=='i' && R2sup=='i' && (R3sup=='n' || R4sup=='n') && (Q1sup=='n' || Q2sup=='n')){
          //invert R1-r2 in metassem
          (metassem[metaScf].end()-1)->orientation='-';
	  (metassem[metaScf].end()-1)->isInv=true;
	  countInversions++;
        } else if(R3sup=='i' && R4sup=='i' && (R1sup=='n' || R2sup=='n') && (Q3sup=='n' || Q4sup=='n')){
          //R3-R4 inversion,add to inversions
          inversions[mum.coords[*qScfReturnPos].scf1][mum.coords[*qScfReturnPos].start1][mum.coords[*qScfReturnPos].end1]='-';
	  countInversions++;
        }
      }
    }

    //Determine whether insertions exist and get insertion positions
    //Ref info
    bool inserRef=false;
    bool NinserRef=false;
    std::vector< int> RinserPos;
    RinserPos.push_back(mum.coords[*curScfAln].end1+1);
    RinserPos.push_back(mum.coords[*qScfReturnPos].start1-1);
    char RstartSup=Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z');
    char RendSup='v';
    int RZstat1i=Rce.search(refScf,RinserPos[0]-1);
    mateAnRow RZstat1;
    if(RZstat1i>=0){
      RZstat1=Rce[refScf][RZstat1i];
    }
    bool Rcov1;
    Rcov1=(RZstat1i>=0 && RZstat1.N >=minCov)?true:false;
    int RZstat2i=Rce.search(refScf,RinserPos[1]+1);
    mateAnRow RZstat2;
    if(RZstat2i>=0){
      RZstat2=Rce[refScf][RZstat2i];
    }
    bool Rcov2;
    Rcov2=(RZstat2i>=0 && RZstat2.N >=minCov)?true:false;
    int RinsLength=RinserPos[1]-RinserPos[0]+1;
    if(RinsLength<0){
      RinsLength=0;
    }
    
    std::vector< int> Rnrow;
    if(RinserPos[1] - RinserPos[0] +1 >= 1){
      inserRef=true;
      NinserRef=RNs.isNinsert(refScf,RinserPos[0],RinserPos[1],minNinsert,minNinsPer);
      if(NinserRef){
        int i=RNs.search(refScf,RinserPos[0],RinserPos[1]);
        Rnrow.push_back(RNs[refScf][i][0]);
        Rnrow.push_back(RNs[refScf][i][1]);
      }
      RendSup=Rce.zsignal(refScf,RinserPos[1]+1,'z');
    } 
    
    //Que info
    bool inserQuery=false;
    bool NinserQuery=false;
    std::vector<int> QinserPos(2,0);
    char QendSup='v'; 
    char QstartSup='v';
    char GlobalOrientQuery=(mum.coords[*curScfAln].start2 <= mum.coords[*qScfReturnPos].end2)?'+':'-';
    //                   s2                e2
    //              ----------------------------- query
    //              ------|---|------|------|---- ref
    //                    
    
    mateAnRow QZstat1;
    mateAnRow QZstat2;
    bool Qcov1,Qcov2;
    
    if(GlobalOrientQuery=='+'){
      QinserPos[0]=(orientQuery=='+')?mum.coords[*curScfAln].end2+1:mum.coords[*curScfAln].start2+1;
      QinserPos[1]=(orientReturnQuery=='+')?mum.coords[*qScfReturnPos].start2-1:mum.coords[*qScfReturnPos].end2-1;
      QstartSup=Qce.zsignal(mum.coords[*curScfAln].scf2,(orientQuery=='+')?mum.coords[*curScfAln].end2:mum.coords[*curScfAln].start2,'z');
      int QZstat1i=Qce.search(mum.coords[*curScfAln].scf2,QinserPos[0]-1);
      if(QZstat1i>=0){
        QZstat1=Qce[mum.coords[*curScfAln].scf2][QZstat1i];
      }
      Qcov1=(QZstat1i>=0 && QZstat1.N >=minCov)?true:false;
      int QZstat2i=Qce.search(mum.coords[*qScfReturnPos].scf2,QinserPos[1]+1);
      if(QZstat2i>=0){
        QZstat2=Qce[mum.coords[*qScfReturnPos].scf2][QZstat2i];
      }
      Qcov2=(QZstat2i>=0 && QZstat2.N >=minCov)?true:false;

      if(QinserPos[1]-QinserPos[0] +1 >= 1){
        inserQuery=true;
        inserQuery=inserQuery && (Qcov1 && Qcov2);
        NinserQuery=QNs.isNinsert(mum.coords[*curScfAln].scf2,QinserPos[0],QinserPos[1],QmaxNinsert,QmaxNinsPer);
        QendSup=Qce.zsignal(mum.coords[*curScfAln].scf2,QinserPos[1],'z');
      }
    } else if(GlobalOrientQuery == '-'){ 
      QinserPos[0]=(orientQuery=='+')?mum.coords[*curScfAln].start2-1:mum.coords[*curScfAln].end2-1;
      QinserPos[1]=(orientReturnQuery=='+')?mum.coords[*qScfReturnPos].end2+1:mum.coords[*qScfReturnPos].start2+1;
      QstartSup=Qce.zsignal(mum.coords[*curScfAln].scf2,(orientQuery=='+')?mum.coords[*curScfAln].start2:mum.coords[*curScfAln].end2,'z');
      int QZstat1i=Qce.search(mum.coords[*curScfAln].scf2,QinserPos[0]+1);
      if(QZstat1i>=0){
        QZstat1=Qce[mum.coords[*curScfAln].scf2][QZstat1i];
      }
      Qcov1=(QZstat1i>=0 && QZstat1.N >=minCov)?true:false;
      int QZstat2i=Qce.search(mum.coords[*qScfReturnPos].scf2,QinserPos[1]-1);
      if(QZstat2i>=0){
        QZstat2=Qce[mum.coords[*qScfReturnPos].scf2][QZstat2i];
      }
      Qcov2=(QZstat2i>=0 && QZstat2.N >=minCov)?true:false;

      if(QinserPos[0]-QinserPos[1] +1 >= 1){
        inserQuery=true;
        inserQuery=inserQuery && (Qcov1 && Qcov2);     
        NinserQuery=QNs.isNinsert(mum.coords[*curScfAln].scf2,QinserPos[0],QinserPos[1],QmaxNinsert,QmaxNinsPer);
        inserQuery=inserQuery && !NinserQuery;
        QendSup=Qce.zsignal(mum.coords[*curScfAln].scf2,QinserPos[1],'z');
      }
    }
    int QinsLength=abs(QinserPos[1]-QinserPos[0])+1;
    
    if(NinserRef && inserQuery 
       && RinsLength <= maxRinsLength
       && QinsLength <= maxQinsLength &&
       ( ( (RZstat1i<0 && RZstat2i<0) &&
           fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&  
           fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat 
         ) ||
         (
           (RZstat1i>=0 || RZstat2i>=0) &&
           (fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
            fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat
           ) &&
           (RZstat1i>=0 && fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) < minZstat) &&
           (RZstat2i>=0 && fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) < minZstat) 
         )
       )
    ){
      //Close primary assembly gap
      MetaRow* newrow; 
      newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                        std::min(QinserPos[0],QinserPos[1]),
                        std::max(QinserPos[0],QinserPos[1]),
                        GlobalOrientQuery,
                        Rce.zsignal(refScf,RinserPos[0],'z'),
                        Rce.zsignal(refScf,RinserPos[1],'z'),
                        QstartSup,QendSup,*curScfAln);
      addN1coord(newrow, mum.coords[*curScfAln].scf2,QinserPos[0],QinserPos[1], GlobalOrientQuery,
                         refScf,RinserPos[0],RinserPos[1],
                         Rnrow[0],Rnrow[1]);
      countNret++;
      countNevents++;	
      
    } else if(inserRef && Rcov1 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) < minZstat 
	               && Rcov2 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) < minZstat){ 
    
      addMetaRow(metaScf,"end",0,refScf,RinserPos[0],RinserPos[1],'+',RstartSup,RendSup,QstartSup,QendSup,*curScfAln);
      
    } else if(!inserRef && Rcov1 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) < minZstat){
      
      return(qScfReturnPos);
      
    }else if(inserQuery && inserRef 
             && RinsLength <= maxRinsLength
             && QinsLength <= maxQinsLength &&
             (
               ( Rcov1 && Rcov2 &&
		 ( (RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat) &&
                   (RZstat2i >= 0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat)
                  ) &&
                 fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
                 fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
		 (
                   (fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                    fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) >= minDiffZ &&
                    fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength))<= minZstat )
                   &&
                   (fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                    fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) >= minDiffZ &&
                    fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) <= minZstat )
                 )
                ) ||
               (
                 (!Rcov1 || !Rcov2)  &&
                 ( (RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat) && 
                   (RZstat2i >= 0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat)
                  ) && 
                 fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
                 fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
		 ( RZstat1i>=0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) >= minDiffZ &&
                                  fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) <= minZstat && 
		   RZstat2i>=0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) >= minDiffZ &&
                                  fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,QinsLength)) <= minZstat
		  ) 
                 &&
                 ( (RinsLength < QinsLength && QinsLength >= minLenSDs*Qce.sd0) ||
                   (RinsLength > QinsLength && RinsLength >= minLenSDs*Rce.sd0) 
                  )
                ) 
              )){
      MetaRow* newrow;
      newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                        std::min(QinserPos[0],QinserPos[1]),
                        std::max(QinserPos[0],QinserPos[1]),
                        GlobalOrientQuery,RstartSup,RendSup,
                        QstartSup,QendSup,*curScfAln);
      addQIcoord(newrow, mum.coords[*curScfAln].scf2, QinserPos[0], QinserPos[1],GlobalOrientQuery, 
                         refScf, RinserPos[0], RinserPos[1]);
      countQI++;
      
    } else if(inserRef && ( (Rcov1 && RZstat1.zvalue(Rce.mu0,Rce.sd0) < minZstat) || (Rcov2 && RZstat2.zvalue(Rce.mu0,Rce.sd0) < minZstat) ) && 
              !inserQuery && Qcov1 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < -1*minZstat ){
      addMetaRow(metaScf,"end",0,refScf,RinserPos[0],RinserPos[1],'+',RstartSup,RendSup,QstartSup,QendSup,*curScfAln);
      
    } else if(inserQuery && !inserRef &&
              QinsLength <= maxQinsLength &&
              //QinsLength >= minLenSDs*Qce.sd0 && QinsLength <= maxQinsLength &&
              (Rcov1 && Rcov2 &&
	       ( (RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat) &&
                 (RZstat2i >= 0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat)
               ) &&
               fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
               fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
               (
                 (fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                  fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) >= minDiffZ &&
                  fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QinsLength))<= minZstat )
                 &&
                 (fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                  fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) >= minDiffZ &&
                  fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) <= minZstat )
               )
              ) ||
              ((!Rcov1 || !Rcov2) &&
               ( (RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat) && //&&
                 (RZstat2i >= 0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat)
               ) &&
               fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
               fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
               ( RZstat1i>=0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) >= minDiffZ &&
                                fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) <= minZstat &&
                 RZstat2i>=0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) >= minDiffZ &&
                                fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,0,QinsLength)) <= minZstat
                )
               &&
               QinsLength >= minLenSDs*Qce.sd0
	      )
             ){
      MetaRow* newrow;
      newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                        std::min(QinserPos[0],QinserPos[1]),
                        std::max(QinserPos[0],QinserPos[1]),
                        GlobalOrientQuery,RstartSup,RendSup,
                        QstartSup,QendSup,*curScfAln);
      addQIcoord(newrow, mum.coords[*curScfAln].scf2, QinserPos[0], QinserPos[1],GlobalOrientQuery, 
                         refScf,RinserPos[0], RinserPos[1]);
      countQI++;
      
    } else if(inserRef && RinsLength <= maxRinsLength //&& QinsLength <= maxQinsLength 
              && !inserQuery && Qcov1 && Qcov2 && (
              ( Rcov1 && Rcov2 &&
		( (RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat) &&
                   (RZstat2i >= 0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat)
                  ) &&
                 fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
                 fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
                 (
                   (fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                    fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) >= minDiffZ &&
                    fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,0))<= minZstat )
                   &&
                   (fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                    fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) >= minDiffZ &&
                    fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) <= minZstat )
                 )
                ) ||
               (
                 (!Rcov1 || !Rcov2)  &&
                 ( (RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat) && //&&
                   (RZstat2i >= 0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) > minZstat)
                  ) &&
                 fabs(QZstat1.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
                 fabs(QZstat2.zvalue(Qce.mu0,Qce.sd0)) < minZstat &&
                 ( RZstat1i>=0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) >= minDiffZ &&
                                  fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) <= minZstat &&
                   RZstat2i>=0 && fabs(RZstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) >= minDiffZ &&
                                  fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RinsLength,0)) <= minZstat
                  )
                ) 
              )){
      MetaRow* newrow;
      newrow=addMetaRow(metaScf,"end",2,"Corr_ref_insertion",
                        QinserPos[0],
                        QinserPos[1],
                        '+',RstartSup,RendSup,
                        QstartSup,QendSup,*curScfAln);
      addQIcoord(newrow, mum.coords[*curScfAln].scf2, QinserPos[0], QinserPos[1],'+',
                         refScf,RinserPos[0], RinserPos[1]);
      countQI++;
      return(qScfReturnPos);
    
    } else if(inserRef){
      addMetaRow(metaScf,"end",0,refScf,RinserPos[0],RinserPos[1],'+',RstartSup,RendSup,QstartSup,QendSup,*curScfAln);
      
    }
    return(qScfReturnPos);
    
  } else { 
    //if queScf maps only once in refScf, does not returns

    //Add right overhang sequences to metassembly or link scaffolds.
    overhangRef=false;
    overhangRefSupport='o';
    overhangQuery=false;
    overhangQuerySupport='o';
    orientQuery=(mum.coords[*curScfAln].start2 <= mum.coords[*curScfAln].end2)?'+':'-';
    qOverhangPos[0]=0;
    qOverhangPos[1]=0;

    if(mum.coords[*curScfAln].end1 < mum.coords[*curScfAln].scfLen1){
      overhangRef=true;
      overhangRefSupport=Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z');
    }
    mateAnRow QZstat1;
    if(orientQuery=='+'){
      if(mum.coords[*curScfAln].end2 < mum.coords[*curScfAln].scfLen2-19 && mum.adjacentAlnQue(queAlns,curScfAln,'+') == queAlns.end() ){
        qOverhangPos[0]=mum.coords[*curScfAln].end2+1;
        qOverhangPos[1]=mum.coords[*curScfAln].scfLen2;
        int QZstat1i=Qce.search(mum.coords[*curScfAln].scf2,qOverhangPos[0]-1);
        if(QZstat1i>=0){
          QZstat1=Qce[mum.coords[*curScfAln].scf2][QZstat1i];
        }
        bool Qcov1;
        Qcov1=(QZstat1i>=0 && QZstat1.N >=minCov)?true:false;
        
        overhangQuery=!QNs.isNinsert(mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],QmaxNinsert,QmaxNinsPer);
        overhangQuery=overhangQuery && Qcov1;
        overhangQuerySupport=Qce.zsignal(mum.coords[*curScfAln].scf2,qOverhangPos[0]-1,'z');
      } 
    } else if(mum.coords[*curScfAln].end2 > 1+19 && mum.adjacentAlnQue(queAlns,curScfAln,'-') == queAlns.end()){
      qOverhangPos[1]=mum.coords[*curScfAln].end2-1;
      qOverhangPos[0]=1;
      int QZstat1i=Qce.search(mum.coords[*curScfAln].scf2,qOverhangPos[1]+1);
      if(QZstat1i>=0){
        QZstat1=Qce[mum.coords[*curScfAln].scf2][QZstat1i];
      }
      bool Qcov1;
      Qcov1=(QZstat1i>=0 && QZstat1.N >=minCov)?true:false;
        
      overhangQuery=!QNs.isNinsert(mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],QmaxNinsert,QmaxNinsPer);
      overhangQuery=overhangQuery && Qcov1;
      overhangQuerySupport=Qce.zsignal(mum.coords[*curScfAln].scf2,qOverhangPos[1]+1,'z');    
    }
    int QovLen=std::max(qOverhangPos[0],qOverhangPos[1])-std::min(qOverhangPos[0],qOverhangPos[1])+1;
    overhangQuery=overhangQuery && QovLen <= maxQinsLength;
 

    if(mum.adjacentAlnRef(rScfAlns,curScfAln,'+') == rScfAlns.end() && !rightLink){ 
      //-- alignment at the end of ref  (either last good or last spurious curScfAln)
      if(overhangRef){
        if(overhangQuery){
	  bool oNinsert=RNs.isNinsert(refScf,mum.coords[*curScfAln].end1+1,std::min(mum.coords[*curScfAln].end1+100,mum.coords[*curScfAln].scfLen1),minNinsert,minNinsPer);
          int RNinsLength=0;
          int NinsIdx=-1;
          if(oNinsert){
                  NinsIdx=RNs.search(refScf,mum.coords[*curScfAln].end1+1,mum.coords[*curScfAln].end1+100);
                  RNinsLength=RNs[refScf][NinsIdx][1]-RNs[refScf][NinsIdx][0]+1;
          }
          int RZstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);
          mateAnRow RZstat1;
          if(RZstat1i>=0){
            RZstat1=Rce[refScf][RZstat1i];
          }
          int RZstat2i=-1;
          if(NinsIdx >= 0) 
            RZstat2i=Rce.search(refScf,RNs[refScf][NinsIdx][1]);
          mateAnRow RZstat2;
          if(RZstat2i>=0)
            RZstat2=Rce[refScf][RZstat2i];
          

          if(oNinsert && RNinsLength <= maxRinsLength && QovLen <= maxQinsLength &&
             (RZstat1i>=0 && fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,RNinsLength,QovLen)) < minZstat) &&
             (RZstat2i>=0 && fabs(RZstat2.newZstat(Rce.mu0,Rce.sd0,RNinsLength,QovLen)) < minZstat)
            ){
            //insert query overhang
            MetaRow* newrow;
            newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,'N','o',overhangQuerySupport,'o',*curScfAln);
            addNAcoord(newrow,refScf,RNs[refScf][NinsIdx][0],RNs[refScf][NinsIdx][1],
                              mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
            //insert Ns
            addMetaRow(metaScf,"end",0,"InsQueOver",1,30,'N','N','o',overhangQuerySupport,'o',*curScfAln);
            countNrows++;
            if(RNs[refScf][NinsIdx][1]+1 < mum.coords[*curScfAln].scfLen1){
	      //insert ref overhang
              addMetaRow(metaScf,"end",0,refScf,RNs[refScf][NinsIdx][1]+1,mum.coords[*curScfAln].scfLen1,'+','N','o',overhangQuerySupport,'o',*curScfAln);
              countNevents++;
              countNA++;            
	    }
          
          } else if(overhangQuerySupport=='n' && RZstat1i >= 0 && fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) > minZstat &&
                    fabs(RZstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QovLen)) >= minDiffZ &&
                    fabs(RZstat1.newZstat(Rce.mu0,Rce.sd0,0,QovLen)) <= minZstat
            ){
              MetaRow* newrow;
              newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,overhangRefSupport,'o',overhangQuerySupport,'o',*curScfAln);
              addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                                 mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
              countBKP++;
          } else {
            //add overhang ref sequence
            addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*curScfAln].scfLen1,'+',overhangRefSupport,'o',overhangQuerySupport,'o',*curScfAln);
          }     
        } else {
			    //No overhang query, overhang ref well supported
            addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*curScfAln].scfLen1,'+',overhangRefSupport,'o','o','o',*curScfAln);
        }
      } else if(overhangQuery && overhangQuerySupport == 'n'){
        //No overhang ref, overhang query well supported 
        MetaRow* newrow; 
        newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,'o','o',overhangQuerySupport,'o',*curScfAln);
        addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                    mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
        countBKP++;
      }
      return(rScfAlns.end());
    } else if(curScfAln < rightAln){
      //Query maps only once and No end of ref scf: simply add/disregard right query overhang
      //or unaligned ref insertion
                  
      //Unaligned insertion in the reference
      std::vector<int>::iterator next=curScfAln;
      do{
        next++;
      }while(mum.coords[*next].start1 >= mum.coords[*curScfAln].start1 &&
             mum.coords[*next].end1 <= mum.coords[*curScfAln].end1 &&
             next < rScfAlns.end() );
      if(next > rightAln){
        std::cerr<<"Error: next>rightAln\n";
        exit(EXIT_FAILURE);
      }

      //get next overhang information
      bool NoverhangQuery=false;
      char NorientQuery=(mum.coords[*next].start2 <= mum.coords[*next].end2)?'+':'-';
      std::vector<int> NqOverhangPos(2,0);
      char NoverhangQuerySupport='o';
  
     
      std::vector<int>& NqueAlns=mum.queryScfs[mum.coords[*next].scf2];
      mateAnRow NQZstat1;
      if(NorientQuery=='+'){                       
                                               // end of query scf
        if(mum.coords[*next].start2 > 1+19 && !mum.QscfBefore(mum.coords[*next].scf2,rScfAlns,next) && mum.adjacentAlnQue(NqueAlns,next,'-') == NqueAlns.end()){
     
          NqOverhangPos[0]=1;
          NqOverhangPos[1]=mum.coords[*next].start2-1;
          int QZstat1i=Qce.search(mum.coords[*next].scf2,NqOverhangPos[1]+1);
          if(QZstat1i >= 0){
            NQZstat1=Qce[mum.coords[*next].scf2][QZstat1i];
          }
          bool Qcov1;
          Qcov1=(QZstat1i >= 0 && NQZstat1.N >= minCov)?true:false;
        
          NoverhangQuery=!QNs.isNinsert(mum.coords[*next].scf2,NqOverhangPos[0],NqOverhangPos[1],QmaxNinsert,QmaxNinsPer);
          NoverhangQuery= NoverhangQuery && Qcov1;
          NoverhangQuerySupport=Qce.zsignal(mum.coords[*next].scf2,NqOverhangPos[1]+1,'z');
        }                                                                                                                         //end of query scf
      }else if(mum.coords[*next].start2 < mum.coords[*next].scfLen2-19 && !mum.QscfBefore(mum.coords[*next].scf2,rScfAlns,next) && mum.adjacentAlnQue(NqueAlns,next,'+') == NqueAlns.end()){
        NqOverhangPos[0]=mum.coords[*next].start2+1;
        NqOverhangPos[1]=mum.coords[*next].scfLen2;
        int QZstat1i=Qce.search(mum.coords[*next].scf2,NqOverhangPos[0]-1);
        if(QZstat1i>=0){
          NQZstat1=Qce[mum.coords[*next].scf2][QZstat1i];
        }
        bool Qcov1;
        Qcov1=(QZstat1i >= 0 && NQZstat1.N >= minCov)?true:false;
        
        NoverhangQuery=!QNs.isNinsert(mum.coords[*next].scf2,NqOverhangPos[0],NqOverhangPos[1],QmaxNinsert,QmaxNinsPer);
        NoverhangQuery=NoverhangQuery && Qcov1;
        NoverhangQuerySupport=Qce.zsignal(mum.coords[*next].scf2,NqOverhangPos[0]-1,'z');
      }
      int NQovLen=std::max(NqOverhangPos[0],NqOverhangPos[1])-std::min(NqOverhangPos[0],NqOverhangPos[1])+1;
      NoverhangQuery=NoverhangQuery && NQovLen <= maxQinsLength;
 
      if(next <=rightAln){
        if(mum.coords[*next].start1 > mum.coords[*curScfAln].end1+19){
          if(RNs.isNinsert(refScf,mum.coords[*curScfAln].end1+1,std::min(mum.coords[*curScfAln].end1+100,mum.coords[*next].start1-1),minNinsert,minNinsPer) && overhangQuery
            ){
            //left N-insert
            std::vector< std::vector< int> >&ranges=RNs[refScf];
            int found=RNs.search(refScf,mum.coords[*curScfAln].end1+1,std::min(mum.coords[*curScfAln].end1+100,mum.coords[*next].start1-1));
            
            if(ranges[found][1]+1 < mum.coords[*next].start1){
              if((RNs.isNinsert(refScf,std::max(ranges[found][1]+1,mum.coords[*next].start1-100),mum.coords[*next].start1-1,minNinsert,minNinsPer) ||
                 RNs.isNinsert(refScf,std::max(mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-100),mum.coords[*next].start1-1,minNinsert,minNinsPer)) &&
                   NoverhangQuery 
                ){
                //left and right N-inserts
                int rightN=RNs.search(refScf,std::max(ranges[found][1]+1,mum.coords[*next].start1-100),mum.coords[*next].start1-1);
                if(rightN <0){
                  rightN=RNs.search(refScf,std::max(mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-100),mum.coords[*next].start1-1);
                }
                if((ranges[rightN][0]-1 - ranges[found][1]+1 +1 >= 1 && RNs.isNinsert(refScf,ranges[found][1]+1,ranges[rightN][0]-1,1,.60)) ||
                   ranges[rightN][0] <= ranges[found][1]
                  ){
                  //no intermediate ref sequence, ref insert can be seen as being just an N-row
                  int Zstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);
                  mateAnRow Zstat1;
                  if(Zstat1i >= 0){
                    Zstat1=Rce[refScf][Zstat1i];
                  }
                  bool cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;
                  
                  int Zstat2i=Rce.search(refScf,mum.coords[*next].start1);
                  mateAnRow Zstat2;
                  if(Zstat2i >= 0){
                    Zstat2=Rce[refScf][Zstat2i];
                  }
                  bool cov2=(Zstat2i >= 0 && Zstat2.N >= minCov)?true:false;
                
                  int QO1length=qOverhangPos[1]-qOverhangPos[0]+1;
                  int QO2length=NqOverhangPos[1]-NqOverhangPos[1]+1;
                  int Nlength=mum.coords[*next].start1-1 - mum.coords[*curScfAln].end1+1 +1;
                
                  if(((Zstat1i < 0 || Zstat2i < 0) && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                                   && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) ||
                     (Zstat1i < 0 && Zstat2i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                                  && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                                  && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO2length)) < minZstat) ||
                     (Zstat1i >= 0 && Zstat2i < 0 &&  QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                                  && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat  
                                                  && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < minZstat) ||
                     (Zstat1i >= 0 && Zstat2i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                                   && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                                   && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length+QO2length)) < minZstat 
                                                   && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length+QO2length)) < minZstat)
                    ){
                    countNevents++;
                    countNevents++;
                    countNA++;
                    countNA++;
                    MetaRow* newrow;                    

                    // insert left overhang
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                      qOverhangPos[0],qOverhangPos[1],orientQuery,
                                      Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                      Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                      overhangQuerySupport,
                                      'o',
                                      *curScfAln);
                    addNAcoord(newrow, refScf,RNs[refScf][found][0],RNs[refScf][found][1],
                                       mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
                    // insert arbitrary Ns
                    int inferredNlen;
                    if(Zstat1i >= 0 && Zstat2i >= 0){
                      inferredNlen=int((Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length) +
                                    Zstat2.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length))
                                    /double(2));                
                      inferredNlen=std::max(30,inferredNlen);
                      
                    }else if(Zstat1i >= 0){
                      inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else if(Zstat2i >= 0){
                      inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    }
                    addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln);
                    countNrows++;
                    // insert right overhang
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                      NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                      Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                      Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                      'o',
                                      NoverhangQuerySupport,
                                      *curScfAln);
                     addNAcoord(newrow,refScf,RNs[refScf][rightN][0],RNs[refScf][rightN][1],
                                       mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],"right",NorientQuery);
                  
                  } else if( (Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat ) || 
                            (Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                          && ( fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < minZstat ) 
                                          && (Zstat2i < 0 || fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO2length)))
                            )
                            ){
                    countNevents++;
                    countNA++;
                    MetaRow* newrow;
                    // insert left overhang
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                        qOverhangPos[0],qOverhangPos[1],orientQuery,
                                        Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                        Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                        overhangQuerySupport,
                                        'o',
                                        *curScfAln);
                    addNAcoord(newrow,refScf,RNs[refScf][found][0],RNs[refScf][found][1],
                               mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
                    // insert arbitrary Ns
                    int inferredNlen;
                    if(Zstat1i >= 0){
                      inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    }
                    addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               Rce.zsignal(refScf,ranges[rightN][0],'z'),
                               Rce.zsignal(refScf,ranges[rightN][1],'z'),
                               *curScfAln);
                    countNrows++;
                  } else if((Zstat2i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                            (Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                          && (fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO2length)) < minZstat )
                            )
                           ){
                    countNevents++;
                    countNA++;
                    MetaRow* newrow;
                    //insert arbitrary Ns
                    int inferredNlen;
                    if(Zstat2i >= 0){
                      inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    }
                    addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               Rce.zsignal(refScf,ranges[rightN][0],'z'),
                               Rce.zsignal(refScf,ranges[rightN][1],'z'),
                               *curScfAln);
                    countNrows++;
                    //insert right overhang
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                      NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                      Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                      Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                      'o',
                                      NoverhangQuerySupport,
                                      *curScfAln);
                    addNAcoord(newrow,refScf,RNs[refScf][rightN][0],RNs[refScf][rightN][1],
                                      mum.coords[*curScfAln].scf2, NqOverhangPos[0],NqOverhangPos[1],"right",NorientQuery);
      
                  } else {
                    //insert ref Ns insert [ mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1 ]
                    addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,
                               mum.coords[*next].start1-1,'+',
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               Rce.zsignal(refScf,ranges[rightN][0],'z'),
                               Rce.zsignal(refScf,ranges[rightN][1],'z'),
                               *curScfAln);
                  }
                } else if(ranges[rightN][0] -1 - ranges[found][1]+1 +1 >= 1){
                  //intermediate ref sequence,
                  int Zstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);
                  mateAnRow Zstat1;
                  if(Zstat1i >= 0){
                    Zstat1=Rce[refScf][Zstat1i];
                  }
                  bool cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;
                  int Zstat2i=Rce.search(refScf,mum.coords[*next].start1);
                  mateAnRow Zstat2;
                  if(Zstat2i >= 0){
                    Zstat2=Rce[refScf][Zstat2i];
                  }
                  bool cov2=(Zstat2i >= 0 && Zstat2.N >= minCov)?true:false;
                  int Nlength1=ranges[found][1]-ranges[found][0]+1;
                  int Nlength2=ranges[rightN][1]-ranges[rightN][0]+1;
                  int QO1length=qOverhangPos[1]-qOverhangPos[0] +1;
                  int QO2length=NqOverhangPos[1]-NqOverhangPos[1] +1;
                  if(((Zstat1i < 0 || Zstat2i < 0) && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) ||
                     (Zstat1i < 0 && Zstat2i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                                  && ( fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength2, QO2length)) < minZstat)) ||
                     (Zstat2i < 0 && Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                                  && ( fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength1,QO1length)) < minZstat)) ||
                     (Zstat1i >= 0 && Zstat2i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                                   && ( fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength1+Nlength2, QO1length+QO2length)) < minZstat 
                                                   && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength1+Nlength2, QO1length+QO2length)) < minZstat ))
                     ){
                    countNevents++;
                    countNevents++;
                    countNB++;
                    countNB++;
                    countNC++;
                    MetaRow* newrow;
                    int inferredNlen;
                    
                    //Insert left ov 
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                      qOverhangPos[0],qOverhangPos[1],orientQuery,
                                      Rce.zsignal(refScf,ranges[found][0],'z'),
                                      Rce.zsignal(refScf,ranges[found][1],'z'),
                                      overhangQuerySupport,
                                      'o',
                                      *curScfAln);
                    addNBcoord(newrow,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                      refScf,RNs[refScf][found][0],RNs[refScf][found][1],
                                      RNs[refScf][found][1]+1,RNs[refScf][rightN][0]-1,
                                      RNs[refScf][rightN][0],"left",orientQuery);
                    //arbitrary Ns,
                    if(Zstat1i >= 0){
                      inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength1,QO1length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    }
                     
                    addMetaRow(metaScf,"end",0,"RNinsertL",1,inferredNlen,'N',
                                      Rce.zsignal(refScf,ranges[found][0],'z'),
                                      Rce.zsignal(refScf,ranges[found][1],'z'),
                                      Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                      Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                      *curScfAln);
                    countNrows++;
                    //intermediate seq,
                    newrow=addMetaRow(metaScf,"end",0,refScf,ranges[found][1]+1,ranges[rightN][0]-1,'+',
                                      Rce.zsignal(refScf,ranges[found][0],'z'),
                                      Rce.zsignal(refScf,ranges[found][1],'z'),
                                      Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                      Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                      *curScfAln);
                    addNCcoord(newrow,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],orientQuery,
                                      refScf,RNs[refScf][found][0],RNs[refScf][found][1],
                                      RNs[refScf][found][1]+1,RNs[refScf][rightN][0]-1,
                                      RNs[refScf][rightN][0],RNs[refScf][rightN][1],
                                      mum.coords[*next].scf2,NqOverhangPos[0],NqOverhangPos[1],NorientQuery
                                      );
                                 
                    //arbitrary Ns,
                    if(Zstat2i >= 0){
                      inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength2,QO2length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    }
                    addMetaRow(metaScf,"end",0,"RNinsertR",1,inferredNlen,'N',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln); 
                    countNrows++;
                    //and right ov
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 'o',
                                 NoverhangQuerySupport,
                                 *curScfAln);
                    addNBcoord(newrow,mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                      refScf,RNs[refScf][rightN][0],RNs[refScf][rightN][1],
                                      RNs[refScf][found][1]+1,RNs[refScf][rightN][0]-1,
                                      RNs[refScf][found][1],"right",NorientQuery);
                  } else if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat)|| 
                            ( Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                           && (fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength1,QO1length)) < minZstat ) 
                                           && (Zstat2i < 0 || fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength1,QO1length)) < fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength2,QO2length))))
                            ){
                    countNevents++;
                    countNB++;
                    MetaRow* newrow;
                    //Insert left ov, 
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                      qOverhangPos[0],qOverhangPos[1],orientQuery,
                                      Rce.zsignal(refScf,ranges[found][0],'z'),
                                      Rce.zsignal(refScf,ranges[found][1],'z'),
                                      overhangQuerySupport,
                                      'o',
                                      *curScfAln);
                    addNBcoord(newrow,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                      refScf,RNs[refScf][found][0],RNs[refScf][found][1],
                                      RNs[refScf][found][1]+1,mum.coords[*next].start1-1,
                                      mum.coords[*next].start1,"left",orientQuery);
                    //arbitrary Ns,
                    int inferredNlen;
                    if(Zstat1i >= 0){
                      inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength1,QO1length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    } 
                    addMetaRow(metaScf,"end",0,"RNinsertL",1,inferredNlen,'N',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln);
                    countNrows++;
                    //and intermediate seq
                    addMetaRow(metaScf,"end",0,refScf,ranges[found][1]+1,mum.coords[*next].start1-1,'+',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln);
                    
                  } else if((Zstat2i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                            ( Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                           && (fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength2,QO2length)) < minZstat)) 
                            ){
                    countNevents++;
                    countNB++;
                    MetaRow* newrow;
                    //insert intermediate seq, 
                    addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,ranges[rightN][0]-1,'+',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln);
                    //arbitrary Ns,
                    int inferredNlen;
                    if(Zstat2i >= 0){
                      inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength2,QO2length));
                      inferredNlen=std::max(30,inferredNlen);
                    }else{
                      inferredNlen=30;
                    }  
                    addMetaRow(metaScf,"end",0,"RNinsertR",1,inferredNlen,'N',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln);
                    countNrows++;
                    //and right ov
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                               NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                               Rce.zsignal(refScf,ranges[rightN][0],'z'),
                               Rce.zsignal(refScf,ranges[rightN][1],'z'),
                               'o',
                               NoverhangQuerySupport,
                               *curScfAln);
                    addNBcoord(newrow,mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                      refScf,RNs[refScf][rightN][0],RNs[refScf][rightN][1],
                                      mum.coords[*curScfAln].end1+1,RNs[refScf][rightN][0]-1,
                                      mum.coords[*curScfAln].end1,"right",NorientQuery);
                    
                  } else {
                    //insert intermidiate seq only
                    addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                 Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                 *curScfAln);
                  }
                  
                }
                
              } else{
                //only left N-insert, followed by a unmapped ref sequence
                int Zstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);
                mateAnRow Zstat1;
                if(Zstat1i >= 0){
                  Zstat1=Rce[refScf][Zstat1i];
                }
                bool cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;
                int Zstat2i=Rce.search(refScf,ranges[found][1]+1);
                mateAnRow Zstat2;
                if(Zstat2i >= 0){
                  Zstat2=Rce[refScf][Zstat2i];
                }
                int Nlength= std::min(ranges[found][1],mum.coords[*next].start1-1) - std::max(mum.coords[*curScfAln].end1 +1, ranges[found][0]) +1;
                int QO1length=qOverhangPos[1]-qOverhangPos[0] +1;
                if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                   (Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                 && ( fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < minZstat)) 
                  ){
                  countNevents++;
                  countNB++;
                  MetaRow* newrow;
                  //insert left overhang
                  newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                    qOverhangPos[0],qOverhangPos[1],orientQuery,
                                    Rce.zsignal(refScf,ranges[found][0],'z'),
                                    Rce.zsignal(refScf,ranges[found][1],'z'),
                                    overhangQuerySupport,
                                    'o',
                                    *curScfAln);
                  addNBcoord(newrow,mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                    refScf,RNs[refScf][found][0],RNs[refScf][found][1],
                                    RNs[refScf][found][1]+1,mum.coords[*next].start1-1,
                                    mum.coords[*next].start1,"left",orientQuery);
                  if(mum.coords[*curScfAln].end1+1 < ranges[found][0]){
                    //insert little refe sequence before N-row
                    addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,ranges[found][0]-1,'+',
                                 Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                                 Rce.zsignal(refScf,ranges[found][0]-1,'z'),
                                 Rce.zsignal(refScf,ranges[found][0],'z'),
                                 Rce.zsignal(refScf,ranges[found][1],'z'),
                                 *curScfAln);
                  }
                  //insert arbitrary Ns
                  int inferredNlen;
                  if(Zstat1i >= 0){
                      inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length));
                      inferredNlen=std::max(30,inferredNlen);
                  }else{
                      inferredNlen=30;
                  } 
                  addMetaRow(metaScf,"end",0,"RNinsertL",1,inferredNlen,'N',
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               Rce.zsignal(refScf,ranges[found][1]+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               *curScfAln);
                  countNrows++;
                } else {
                  addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,ranges[found][1],'+',
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               Rce.zsignal(refScf,ranges[found][1]+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               *curScfAln);
                }
                if(ranges[found][1]+1 < mum.coords[*next].start1 && !RNs.isNinsert(refScf,ranges[found][1]+1,mum.coords[*next].start1-1,1,0.85)){
                  //insert right unmapped ref seq
                  addMetaRow(metaScf,"end",0,refScf,ranges[found][1]+1,mum.coords[*next].start1-1,'+',
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               Rce.zsignal(refScf,ranges[found][1]+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               *curScfAln);
                }
                //insert right ov?
                if(NoverhangQuery){
                  Zstat1i=Rce.search(refScf,mum.coords[*next].start1);
                  if(Zstat1i >= 0){
                    Zstat1=Rce[refScf][Zstat1i];
                  }
                  cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;  
                  int QO2length=NqOverhangPos[1] - NqOverhangPos[2] +1;
                  if((Zstat1i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                      (Zstat1i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                    && (Zstat1.zvalue(Rce.mu0,Rce.sd0) < -1*minZstat 
                                    && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) < minZstat 
                                    && fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) >= minDiffZ) )
                    ){
                    MetaRow* newrow;
                    //insert arbitrary ns
                    addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                                 Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                 Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                 'o',
                                 NoverhangQuerySupport,
                                 *curScfAln);
                    countNrows++;
                    //insert right ov
                    newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                      NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                      Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                      Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                      'o',
                                       NoverhangQuerySupport,
                                      *curScfAln);
                    addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                                       mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                       "right",NorientQuery);
                    countBKP++;
                  }
                }
              } 
            } else {
            // ref insert is entirely an N-row.
              int Zstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);
              mateAnRow Zstat1;
              if(Zstat1i >= 0){
                Zstat1=Rce[refScf][Zstat1i];
              }
              bool cov1=(Zstat1i>=0 && Zstat1.N>=minCov)?true:false;
              
              int Nlength=mum.coords[*next].start1-1 - mum.coords[*curScfAln].end1+1 +1;
              int QO1length=qOverhangPos[1]-qOverhangPos[0]+1;
              
              if(NoverhangQuery){
                // Left and right overhangs  
                int Zstat2i=Rce.search(refScf,mum.coords[*next].start1);
                mateAnRow Zstat2;
                if(Zstat2i >= 0){
                  Zstat2=Rce[refScf][Zstat2i];
                }
                bool cov2=(Zstat2i>=0 && Zstat2.N>=minCov)?true:false;
                
                int QO2length=NqOverhangPos[1]-NqOverhangPos[1]+1;
                
                if(((Zstat1i < 0 || Zstat2i < 0) && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat  
                                                 && NQZstat1.zvalue(Rce.mu0,Rce.sd0) < minZstat) ||
                   (Zstat1i < 0 && Zstat2i >= 0  && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat   
                                                 && NQZstat1.zvalue(Rce.mu0,Rce.sd0) < minZstat
                                                 && (fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO2length)) < minZstat)) ||
                   (Zstat2i < 0 && Zstat2i >= 0  && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat   
                                                 && NQZstat1.zvalue(Rce.mu0,Rce.sd0) < minZstat
                                                 && (fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < minZstat)) ||
                   (Zstat1i >= 0 && Zstat2i >= 0  && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat   
                                                  && NQZstat1.zvalue(Rce.mu0,Rce.sd0) < minZstat
                                                  && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length+QO2length)) < minZstat 
                                                  && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length+QO2length)) < minZstat 
                                     )
                  ){
                  countNevents++;
                  countNevents++;
                  countNA++;
                  countNA++;  
                  MetaRow* newrow;
                  // insert left overhang
                  newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                    qOverhangPos[0],qOverhangPos[1],orientQuery,
                                    Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                    Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                    overhangQuerySupport,
                                    'o',
                                    *curScfAln);
                  addNAcoord(newrow, refScf,mum.coords[*curScfAln].end1+1, mum.coords[*next].start1-1,
                                     mum.coords[*curScfAln].scf2,qOverhangPos[0], qOverhangPos[1],"left",orientQuery);
                  // insert arbitrary Ns
                  int inferredNlen;
                  if(Zstat1i >= 0 && Zstat2i >= 0){
                    inferredNlen=int((Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length) +
                                      Zstat2.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length)
                                       )/double(2));
                    inferredNlen=std::max(30,inferredNlen);
                  }else if(Zstat1i >= 0){
                    inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length));
                    inferredNlen=std::max(30,inferredNlen);
                  }else if(Zstat2i >= 0){
                    inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength,QO1length+QO2length));
                    inferredNlen=std::max(30,inferredNlen);
                  }else{
                    inferredNlen=30;
                  }
                  addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                               Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               *curScfAln);
                  countNrows++;
                  // insert right overhang
                  newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                    NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                    Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                    Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                    'o',
                                    NoverhangQuerySupport,
                                    *curScfAln);
                  addNAcoord(newrow,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,
                                    mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],"right",NorientQuery);
                  
                } else if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                          ( Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                         && ((Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length) < minZstat) &&  
                                             (Zstat2i < 0 ||
                                              fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO2length)))
                                            )
                          )
                         ){
                  countNevents++;
                  countNA++;
                  MetaRow* newrow;
                  // insert left overhang
                  newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                    qOverhangPos[0],qOverhangPos[1],orientQuery,
                                    Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                    Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                    overhangQuerySupport,
                                    'o',
                                    *curScfAln);
                  addNAcoord(newrow,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,
                                    mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
                  // insert arbitrary Ns
                  int inferredNlen;
                  if(Zstat1i >= 0){
                    inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length));
                    inferredNlen=std::max(30,inferredNlen);
                  }else{
                    inferredNlen=30;
                  }
                  addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                               Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               *curScfAln);
                  countNrows++;
                } else if((Zstat2i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                          (Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                        && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QO2length)) < minZstat )
                         ){
                  countNevents++;
                  countNA++;
                  MetaRow* newrow;
                  // insert arbitrary Ns
                  int inferredNlen;
                  if(Zstat2i >= 0){
                    inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength,QO2length));
                    inferredNlen=std::max(30,inferredNlen);
                  }else{
                    inferredNlen=30;
                  }
                  addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                               Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               *curScfAln);
                  countNrows++;
                  // insert right overhang
                  newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                    NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                    Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                    Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                    'o',
                                    NoverhangQuerySupport,
                                    *curScfAln);
                  addNAcoord(newrow,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,
                                    mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],"right",NorientQuery);
                } else {
                  //insert ref Ns
                  addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                               Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               *curScfAln);
                }
                
              } else {
                // Only left overhang
                if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                    (Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                  && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QO1length)) < minZstat ) 
                  ){
                  countNevents++;
                  countNA++;  
                  MetaRow* newrow;
                  // insert left overhang
                  newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                    qOverhangPos[0],qOverhangPos[1],orientQuery,
                                    Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                    Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                    overhangQuerySupport,
                                    'o',
                                    *curScfAln);
                  addNAcoord(newrow,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,
                                    mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],"left",orientQuery);
                  // insert arbitrary Ns
                  int inferredNlen;
                  if(Zstat1i >= 0){
                    inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QO1length));
                    inferredNlen=std::max(30,inferredNlen);
                  }else{
                    inferredNlen=30;
                  }
                  addMetaRow(metaScf,"end",0,"RNinsert",1,inferredNlen,'N',
                               Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               *curScfAln);
                  countNrows++;
                  
                } else {
                  //insert reference Ns
                  addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                               Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                               Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                               Rce.zsignal(refScf,ranges[found][0],'z'),
                               Rce.zsignal(refScf,ranges[found][1],'z'),
                               *curScfAln);
                }
                
              }  
            }
          } else if(RNs.isNinsert(refScf,std::max(mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-100),mum.coords[*next].start1-1,minNinsert,minNinsPer) &&
                    NoverhangQuery
                    ){
            //only right N-insert, preceded by a unmapped reference sequence
            std::vector< std::vector< int> >&ranges=RNs[refScf];
            int rightN=RNs.search(refScf,std::max(mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-100),mum.coords[*next].start1-1);
            //insert left ov?
            if(overhangQuery){
              int Zstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);     
              mateAnRow Zstat1;
              if(Zstat1i >= 0){
                Zstat1=Rce[refScf][Zstat1i];
              }
              bool cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;
              int QOlength1=qOverhangPos[1]-qOverhangPos[0]+1;
              if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat ) || 
                 (Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                               && Zstat1.zvalue(Rce.mu0,Rce.sd0) < -1*minZstat 
                               && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QOlength1)) < minZstat 
                               && fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QOlength1)) >= minDiffZ)
                ){
                MetaRow* newrow;
                //insert left ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                  qOverhangPos[0],qOverhangPos[1],orientQuery,
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  overhangQuerySupport,
                                  'o',
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                            mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                   "left",orientQuery);
                countBKP++;
                //insert arbitrary Ns
                addMetaRow(metaScf,"end",0,"QOLNinsert",1,30,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].start1,'z'),
                           Rce.zsignal(refScf,mum.coords[*curScfAln].start1,'z'),
                           overhangQuerySupport,
                           'o',
                           *curScfAln);
                countNrows++;
              }
            }
            if(mum.coords[*curScfAln].end1+1 < ranges[rightN][0] && !RNs.isNinsert(refScf,mum.coords[*curScfAln].end1+1,ranges[rightN][0]-1,1,0.85)){
              //insert left unmapped reference sequence
              addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,ranges[rightN][0]-1,'+',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                           Rce.zsignal(refScf,ranges[rightN][0]-1,'z'),
                           Rce.zsignal(refScf,ranges[rightN][0],'z'),
                           Rce.zsignal(refScf,ranges[rightN][1],'z'),
                           *curScfAln);
            }
            int Zstat1i=Rce.search(refScf,ranges[rightN][0]-1);
            mateAnRow Zstat1;
            if(Zstat1i >= 0){
              Zstat1=Rce[refScf][Zstat1i];
            }
            bool cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;
            int Zstat2i=Rce.search(refScf,mum.coords[*next].start1);
            mateAnRow Zstat2;
            if(Zstat2i >= 0){
              Zstat2=Rce[refScf][Zstat2i];
            }
            bool cov2=(Zstat2i >= 0 && Zstat2.N >= minCov)?true:false;
            int QOlength=NqOverhangPos[1]-NqOverhangPos[0]+1;
            int Nlength=std::min(mum.coords[*next].start1-1,ranges[rightN][1])- std::max(mum.coords[*curScfAln].end1+1,ranges[rightN][0]);
            if((Zstat1i < 0 && Zstat2i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0)< minZstat) || 
               (Zstat2i < 0 && Zstat1i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                            && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QOlength)) < minZstat) ||
               (Zstat1i < 0 && Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                            && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,Nlength,QOlength)) < minZstat ) ||
               (Zstat1i >= 0 && Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat  
                                             && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,Nlength,QOlength)) < minZstat ) 
              ){
              countNevents++;
              countNB++;
              MetaRow* newrow;
              //insert arbitrary Ns
              int inferredNlen;
              if(Zstat1i >= 0 && Zstat2i >= 0){
                inferredNlen=int((Zstat1.newDeltaMean(Rce.mu0,Nlength,QOlength) +
                                    Zstat2.newDeltaMean(Rce.mu0,Nlength,QOlength)
                                    )/double(2));
                inferredNlen=std::max(30,inferredNlen);
                    
              }else if(Zstat1i >= 0){
                inferredNlen=int(Zstat1.newDeltaMean(Rce.mu0,Nlength,QOlength));
                inferredNlen=std::max(30,inferredNlen);
              }else if(Zstat2i >= 0){
                inferredNlen=int(Zstat2.newDeltaMean(Rce.mu0,Nlength,QOlength));
                inferredNlen=std::max(30,inferredNlen);
              }else{
                inferredNlen=30;
              }
              addMetaRow(metaScf,"end",0,"RNinsertR",1,inferredNlen,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                           Rce.zsignal(refScf,ranges[rightN][0]-1,'z'),
                           Rce.zsignal(refScf,ranges[rightN][0],'z'),
                           Rce.zsignal(refScf,ranges[rightN][1],'z'),
                           *curScfAln);
              countNrows++;
              if(ranges[rightN][1]+1 < mum.coords[*next].start1){
                //insert little refe sequence after right N-row
                addMetaRow(metaScf,"end",0,refScf,ranges[rightN][1]+1,mum.coords[*next].start1-1,'+',
                             Rce.zsignal(refScf,ranges[rightN][0],'z'),
                             Rce.zsignal(refScf,ranges[rightN][1],'z'),
                             Rce.zsignal(refScf,ranges[rightN][1]+1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1-1,'z'),
                             *curScfAln);
              }
              //insert right ov
              newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                Rce.zsignal(refScf,ranges[rightN][0],'z'),
                                Rce.zsignal(refScf,ranges[rightN][1],'z'),
                                'o',
                                NoverhangQuerySupport,
                                *curScfAln);    
              addNBcoord(newrow,mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                refScf,RNs[refScf][rightN][0],RNs[refScf][rightN][1],
                                mum.coords[*curScfAln].end1+1,RNs[refScf][rightN][0]-1,
                                mum.coords[*curScfAln].end1,"right",NorientQuery);
            } else {
              //insert reference Ns [ std::max(mum.coords[*curScfAln].end1+1 , ranges[rightN][0]), mum.coords[*next].start-1 ]
              addMetaRow(metaScf,"end",0,refScf,ranges[rightN][0],mum.coords[*next].start1-1,'+',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1+1,'z'),
                           Rce.zsignal(refScf,ranges[rightN][0]-1,'z'),
                           Rce.zsignal(refScf,ranges[rightN][0],'z'),
                           Rce.zsignal(refScf,ranges[rightN][1],'z'),
                           *curScfAln); 
            }
          }else{ //No right nor left N-insert
            int Zstat1i=Rce.search(refScf,mum.coords[*curScfAln].end1);
            mateAnRow Zstat1;
            if(Zstat1i >= 0){
              Zstat1=Rce[refScf][Zstat1i];
            }
            bool cov1=(Zstat1i >= 0 && Zstat1.N >= minCov)?true:false;
            int QO1length=std::max(0,qOverhangPos[1]-qOverhangPos[0]+1);
            
            int Zstat2i=Rce.search(refScf,mum.coords[*next].start1);
            mateAnRow Zstat2;
            if(Zstat2i >= 0){
              Zstat2=Rce[refScf][Zstat2i];
            }
            bool cov2=(Zstat2i >= 0 && Zstat2.N >= minCov)?true:false;
            int QO2length=std::max(0,NqOverhangPos[1]-NqOverhangPos[0]+1);
            
            if(overhangQuery && NoverhangQuery){
  
              if(((Zstat1i < 0 || Zstat2i < 0) && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                               && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) ||
                 (Zstat1i < 0 && Zstat2i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                              && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                              && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0, QO2length)) < minZstat
                                              && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0, QO2length)) <= fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) 
                                              && fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0, QO2length)) >= minDiffZ ) ||
                 (Zstat2i < 0 && Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                              && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                              && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) < minZstat 
                                              && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) <= fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) 
                                              && fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) >= minDiffZ) ||
                 (Zstat1i >= 0 && Zstat2i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                               && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                               && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0, QO1length+QO2length)) < minZstat 
                                               && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0, QO1length+QO2length)) <= fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) 
                                               && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0, QO1length+QO2length)) < fabs(minZstat) 
                                               && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0, QO1length+QO2length)) <= fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) 
                                               && fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0, QO1length+QO2length)) >= minDiffZ 
                                               && fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0, QO1length+QO2length)) >= minDiffZ
                  )
                ){
                MetaRow* newrow;
                //add left ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                  qOverhangPos[0],qOverhangPos[1],orientQuery,
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  overhangQuerySupport,
                                  'o',
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                            mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                                              "left",orientQuery);
                countBKP++;
                //add Ns
                addMetaRow(metaScf,"end",0,"QOL",1,30,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           overhangQuerySupport,
                           'o',
                           *curScfAln);
                countNrows++;
                //add RUI
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             *curScfAln);
                //add Ns
                addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           'o',
                           NoverhangQuerySupport,
                           *curScfAln);
                countNrows++;
                //add right ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                  NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                  Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                  'o',
                                  NoverhangQuerySupport,
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                            mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                                            "right",NorientQuery);
                countBKP++;
              } else if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                        (Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                      && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) < minZstat  
                                      && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) <= fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0))  
                                      && fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) >= minDiffZ
                                      && (Zstat2i < 0 ||
                                         fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) < fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length))
                                         )
                        )
                 ){
                MetaRow* newrow;
                //add left ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                  qOverhangPos[0],qOverhangPos[1],orientQuery,
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  overhangQuerySupport,
                                  'o',
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                            mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                                               "left",orientQuery);
                countBKP++;
                //add Ns
                addMetaRow(metaScf,"end",0,"QOL",1,30,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           overhangQuerySupport,
                           'o',
                           *curScfAln);
                countNrows++;
                //add RUI 
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           *curScfAln);
                             
              } else if((Zstat2i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat ) || 
                        (Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                      && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) < minZstat 
                                      && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) <= fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0))
                                      && fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) >= minDiffZ)
                 ){
                MetaRow* newrow;
                //add RUI
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           *curScfAln);
                //add Ns
                addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           'o',
                           NoverhangQuerySupport,
                           *curScfAln);
                countNrows++;
                //add right ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                  NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                  Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                  'o',
                                  NoverhangQuerySupport,
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                            mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                                            "right",NorientQuery);
                countBKP++;
              } else {
                //add RUI
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           *curScfAln);
              }
  
            } else if(overhangQuery){
  
              if((Zstat1i < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                 ( Zstat1i >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) < minZstat 
                                && fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) <= fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0))  
                                && fabs(Zstat1.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat1.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) >= minDiffZ
                 )
                 ){
                MetaRow* newrow;
                //add left ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                  qOverhangPos[0],qOverhangPos[1],orientQuery,
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                  overhangQuerySupport,
                                  'o',
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                            mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                                               "left",orientQuery);
                countBKP++;
               //add Ns
               addMetaRow(metaScf,"end",0,"QOL",1,30,'N',
                          Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                          Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                          overhangQuerySupport,
                          'o',
                          *curScfAln);
               countNrows++;
               //add RUI
               addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                            Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                            Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                            Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                            Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                            *curScfAln);
              } else {
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             *curScfAln);
              }
  
            } else if(NoverhangQuery){
              if((Zstat2i < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat )|| 
                 (Zstat2i >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                               && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) < minZstat 
                               && fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) <= fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) 
                               && fabs(Zstat2.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat2.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) >= minDiffZ
                          )
                 ){
                MetaRow* newrow;
                //add RUI
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             *curScfAln);
                //add Ns
                addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                           Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                           Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                           'o',
                           NoverhangQuerySupport,
                           *curScfAln);
                countNrows++;
                //add right ov
                newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                  NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                  Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                  Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                  'o',
                                  NoverhangQuerySupport,
                                  *curScfAln);
                addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                                   mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                   "right",NorientQuery);
                countBKP++;
              }else {
                //add RUI
                addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                             Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                             *curScfAln);
              }
  
            } else {
              //add RUI
              addMetaRow(metaScf,"end",0,refScf,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                         Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                         Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                         Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                         Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                         *curScfAln);
            }
          }  	
        } else if(mum.coords[*next].start1 - mum.coords[*curScfAln].end1 > 1 && mum.coords[*next].start1 - mum.coords[*curScfAln].end1 <= 19){
			  //add little ref insertion.
          addMetaRow(metaScf,"end",0,mum.coords[*curScfAln].scf1,mum.coords[*curScfAln].end1+1,mum.coords[*next].start1-1,'+',
                     Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                     Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                     Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                     Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                     *curScfAln);
        } else {
        //No ref insertion with negative Zstat
          int Zstati=Rce.search(refScf,mum.coords[*curScfAln].end1);
          mateAnRow Zstat;
          if(Zstati >= 0){
            Zstat=Rce[refScf][Zstati];
          }
          bool cov=(Zstati >= 0 && Zstat.N >= minCov)?true:false;
          int QO1length=0;
          int QO2length=0;
          if(overhangQuery){
            QO1length=qOverhangPos[1]-qOverhangPos[0]+1;
          }
          if(NoverhangQuery){
            QO2length=NqOverhangPos[1]-NqOverhangPos[0]+1;
          }
          if(Zstati >= 0 && Zstat.zvalue(Rce.mu0,Rce.sd0) >= -1*minZstat){
           return(curScfAln+1);
          } else if((Zstati < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat && 
                             NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat ) ||
                     (Zstati >= 0 && overhangQuery && NoverhangQuery
                                  && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                  && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                  && fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length+QO2length)) < minZstat 
                                  && fabs(Zstat.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length+QO2length)) >= minDiffZ)
                    ){
            if(overhangQuery){
              MetaRow* newrow;
              //insert left overhang
              newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                qOverhangPos[0],qOverhangPos[1],orientQuery,
                                Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                overhangQuerySupport,
                                'o',
                                *curScfAln);
              addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                                 mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                 "left",orientQuery);
              countBKP++;
            }
            if(overhangQuery || NoverhangQuery){
              //insert arbitrary Ns
              addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                           overhangQuerySupport,
                           'o',
                           'o',
                           NoverhangQuerySupport,
                           *curScfAln);
              countNrows++;
            }
            if(NoverhangQuery){
              MetaRow* newrow;
              //insert right overhang
              newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                'o',
                                NoverhangQuerySupport,
                                *curScfAln);
              addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                                 mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                 "right",NorientQuery);
              countBKP++;
            }
          } else if(overhangQuery && NoverhangQuery){
            if((Zstati < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
               (Zstati >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                            && fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) < minZstat 
                            && fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) <= fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO2length))
                            && fabs(Zstat.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) >= minDiffZ )
                     ){
              MetaRow* newrow;
              //insert left overhang
              newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                qOverhangPos[0],qOverhangPos[1],orientQuery,
                                Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                overhangQuerySupport,
                                'o',
                                *curScfAln);
              addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                                 mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                 "left",orientQuery);
              countBKP++;
              //insert arbitrary Ns
              addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                         overhangQuerySupport,
                         'o',
                         'o',
                         NoverhangQuerySupport,
                         *curScfAln);
              countNrows++;
            } else if((Zstati < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat ) || 
                      (Zstati >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                   && fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) < minZstat
                                   && fabs(Zstat.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) >= minDiffZ )
                      ){
              MetaRow* newrow;
              //insert arbitrary Ns
              addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                         overhangQuerySupport,
                         'o',
                         'o',
                         NoverhangQuerySupport,
                         *curScfAln);
              countNrows++;
              //insert right overhang
              newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                                NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                                Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                                'o',
                                NoverhangQuerySupport,
                                *curScfAln);
              addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                                 mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                                 "right",NorientQuery);
              countBKP++;
            } //                Aqui
          } else if(overhangQuery && 
                    ( (Zstati < 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                      (Zstati >= 0 && QZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat
                                   && fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) < minZstat 
                                   && fabs(Zstat.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO1length)) >= minDiffZ ) ) 
                    ){
             MetaRow* newrow;  
             //insert left overhang
             newrow=addMetaRow(metaScf,"end",1,mum.coords[*curScfAln].scf2,
                                qOverhangPos[0],qOverhangPos[1],orientQuery,
                                Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                Rce.zsignal(refScf,mum.coords[*curScfAln].end1,'z'),
                                overhangQuerySupport,
                                'o',
                                *curScfAln);
             addBKPcoord(newrow,refScf,mum.coords[*curScfAln].end1,
                                mum.coords[*curScfAln].scf2,qOverhangPos[0],qOverhangPos[1],
                                 "left",orientQuery);
             countBKP++;
             //insert arbitrary Ns
             addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                        overhangQuerySupport,
                        'o',
                        'o',
                        NoverhangQuerySupport,
                        *curScfAln);
             countNrows++;
          } else if(NoverhangQuery && 
                    ((Zstati < 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat) || 
                     (Zstati >= 0 && NQZstat1.zvalue(Qce.mu0,Qce.sd0) < minZstat 
                                  && fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) < minZstat
                                  && fabs(Zstat.zvalue(Rce.mu0,Rce.sd0)) - fabs(Zstat.newZstat(Rce.mu0,Rce.sd0,0,QO2length)) >= minDiffZ ))
                   ){
            MetaRow* newrow;
            //insert arbitrary Ns
            addMetaRow(metaScf,"end",0,"QOR",1,30,'N',
                       overhangQuerySupport,
                       'o',
                       'o',
                       NoverhangQuerySupport,
                       *curScfAln);
            countNrows++;
            //insert right overhang
            newrow=addMetaRow(metaScf,"end",1,mum.coords[*next].scf2,
                              NqOverhangPos[0],NqOverhangPos[1],NorientQuery,
                              Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                              Rce.zsignal(refScf,mum.coords[*next].start1,'z'),
                              'o',
                              NoverhangQuerySupport,
                              *curScfAln);
            addBKPcoord(newrow,refScf,mum.coords[*next].start1,
                               mum.coords[*curScfAln].scf2,NqOverhangPos[0],NqOverhangPos[1],
                               "right",NorientQuery);
            countBKP++;
          } 
        }
      }
    }
    return(curScfAln+1);   
  } 
}

//---------------------------------- printStats ------------------------------//
void Metassembly::printStats()
{
  std::ofstream OUTR (std::string(outp).append(".report").c_str());

  std::map<std::string, int>::iterator it;
  std::vector< int> lengths;

  double lmean=0;
  double lsd=0;
  int lmin=0;
  int lmax=0;
  int lmedian=0;
  long ltotal=0;
  int N50=0;

  //Compute mean, sd, and total length
  for(it=metaLengths.begin(); it != metaLengths.end(); it++){
    lengths.push_back(it->second);
    lmean+=it->second;
    lsd+=pow(it->second,2);
  }
  ltotal=( long)lmean;

  //Compute median
  std::sort(lengths.begin(),lengths.end());
  if((int)lengths.size() % 2 == 0)
    lmedian=lengths[(lengths.size()/2)+1];
  else
    lmedian=lengths[(int)ceil((float)lengths.size()/(float)2)];

  //Compute N50
  std::vector< int>::reverse_iterator rlit;
  long lsum=0;
  for(rlit=lengths.rbegin(); rlit<lengths.rend(); rlit++){
    lsum+=(long)*rlit;
    if(lsum >= ltotal/double(2)){
      N50=*rlit;
      break;
    }
  }

  lmean=lmean/(double)lengths.size();
  lsd=lsd/(double) lengths.size() - pow(lmean,2);
  lsd=sqrt(lsd);
  lmin=lengths[0];
  lmax=lengths.back();
  
  OUTR << "\n\n-----Metassembly "<<outp<<" statistics-----"<<std::endl
          << "---------------------------------"<<std::endl;
  OUTR.precision(3);
  OUTR << "Scaffolds: "<<lengths.size()<<"\n"
      << "length mean:    "<<std::fixed<<lmean<<"\n"
      << "length sd:      "<<std::fixed<<lsd<<"\n"
      << "length min:     "<<lmin<<"\n"
      << "length max:     "<<lmax<<"\n"
      << "length med:     "<<lmedian<<"\n"
      << "N50:            "<<N50<<"\n"
      << "total length:   "<<ltotal<<"\n"
      << std::endl;

  OUTR << "\n\n-----Events---------------"<<std::endl
       << "---------------------------------"<<std::endl
       << "----Ns:"
       << "Total number of events:  "<<countNevents<<"\n"
       << "Case 1:                  " << countNret << "\n"
       << "Case 2.A:                "<< countNA<<"\n"
       << "Case 2.B:                "<< countNB<<"\n"
       << "Case 2.C:                "<< countNC<<"\n"
       << "----\n"
       << "BKP:                     "<< countBKP<<"\n"
       << "INV:                     "<< countInversions<<"\n"
       << "Insertions:              "<< countQI <<"\n"
       << "Scf Joined-Links:        "<< countLN <<" - "<<countLN/2<<"\n";
}

//---------------------------------- printLinks ------------------------------//
void Metassembly::printLinks()
{
  std::ofstream OUTL (std::string(outp).append(".links").c_str());
  std::map<std::string, MetaLink>::iterator it;
  for(it=scfsLinks.begin(); it != scfsLinks.end(); it++){
    it->second.print(OUTL);
  }
}

//---------------------------------- printMetaRow ----------------------------//
std::ofstream& Metassembly::printMetaRow(std::ofstream& OUTF, std::string metScf, int mstart, int mend,
                            MetaRow& mrow
                           )
{

  std::string scf=mrow.scf;
  if(scf != "InsQueOver") 
    scf.erase(0,2);
  
  OUTF << metScf << "\t"
       << mstart << "\t"
       << mend << "\t"
       << mrow.file << "\t"
       << scf << "\t"
       << mrow.start << "\t"
       << mrow.end << "\t"
       << mrow.orientation << "\t"
       << mrow.sigRef1 << "\t"
       << mrow.sigRef2 << "\t"
       << mrow.sigQue1 << "\t"
       << mrow.sigQue2 << "\t"
       << mum.coords[mrow.aln].scf1 << "\t"
       << mum.coords[mrow.aln].start1 << "\t"
       << mum.coords[mrow.aln].end1 << "\t"
       << mum.coords[mrow.aln].scf2 << "\t"
       << mum.coords[mrow.aln].start2 << "\t"
       << mum.coords[mrow.aln].end2
       << std::endl;

  return(OUTF);
}

//---------------------------------- printMetassem --------------------------//
bool compareName_Len(name_len A, name_len B){return(A.length > B.length);}

void Metassembly::printMetassem()
{

  std::ofstream OUTM (std::string(outp).append(".metassem").c_str());
  std::ofstream OUTINV;
  std::ofstream OUTNA;
  std::ofstream OUTNB;
  std::ofstream OUTNC;
  std::ofstream OUTN1;
  std::ofstream OUTQI;
  std::ofstream OUTBKP;
  std::ofstream OUTLN;

  if(outp.find(".noLink") == std::string::npos){
    OUTINV.open(std::string(outp).append(".INVs").c_str());
    OUTNA.open(std::string(outp).append(".NAs").c_str());
    OUTNB.open(std::string(outp).append(".NBs").c_str());
    OUTNC.open(std::string(outp).append(".NCs").c_str());
    OUTN1.open(std::string(outp).append(".N1s").c_str());
    OUTQI.open(std::string(outp).append(".QIs").c_str());
    OUTBKP.open(std::string(outp).append(".BKPs").c_str());
    OUTLN.open(std::string(outp).append(".LNs").c_str());
  }

  std::map<std::string, std::vector<MetaRow> >::iterator itm;
  std::vector<MetaRow>::iterator coords;


  //sort metassembly scaffolds by decreasing length
  std::vector<name_len> nls;
  name_len nl;
  for(itm=metassem.begin(); itm != metassem.end(); itm++){
    nl.name=(*itm).first;
    nl.length=metaLengths[(*itm).first];
    nls.push_back(nl);
  }
  std::sort(nls.begin(),nls.end(),compareName_Len);
  std::vector<name_len>::iterator it;
  int i=0;
  for(it=nls.begin();it<nls.end();it++){
    std::string metScf=std::string(mscfPrefix).append("_");
    std::ostringstream c;
    c << i;
    metScf.append(c.str());
    OUTM << "@SQ\tSN:"<<metScf<<"\tLN:"<<(*it).length<<std::endl;
    i++;
  }
  i=0;
  for(it=nls.begin();it < nls.end(); it++){
    int sum=0;
    int start;
    int end;
    std::string metScf=std::string(mscfPrefix).append("_");
    std::ostringstream c;
    c << i;
    metScf.append(c.str());
    for(coords=metassem[it->name].begin(); coords != metassem[it->name].end(); coords++){
      if(coords->mQIc != NULL && coords -> scf == "Corr_ref_insertion"){
	//Deal with ref deletions
        coords->mQIc->print(metScf,end,end,&(*coords),this,OUTQI);
	continue;
      }
      start=sum+1;
      end=start+ coords->end - coords->start;
      sum+=end-start+1;
      printMetaRow(OUTM, metScf, start, end, *coords);

      if(outp.find(".noLink") != std::string::npos)
        continue;

      if(coords->isInv){
        printMetaRow(OUTINV, metScf, start, end, *coords);
      }

      if(coords->mlnc != NULL){
//        std::cerr<<"mlnc1"<<std::endl;
        coords->mlnc->print(metScf,start,end,&(*coords),this,OUTLN);
//        std::cerr<<"mlnc2"<<std::endl;
      }

      if(coords->mNAc != NULL){
//        std::cerr<<"mNAc1"<<std::endl;
        coords->mNAc->print(metScf,start,end,&(*coords),this,OUTNA);
//        std::cerr<<"mNAc2"<<std::endl;
      }

      if(coords->mNBc != NULL){
//        std::cerr<<"mNBc1"<<std::endl;
        coords->mNBc->print(metScf,start,end,&(*coords),this,OUTNB);
//        std::cerr<<"mNBc2"<<std::endl;
      }

      if(coords->mNCc != NULL){
//        std::cerr<<"mNCc1"<<std::endl;
        coords->mNCc->print(metScf,start,end,&(*coords),this,OUTNC);
//        std::cerr<<"mNCc1"<<std::endl;
      }

      if(coords->mN1c != NULL){
//        std::cerr<<"mN1c1"<<std::endl;
        coords->mN1c->print(metScf,start,end,&(*coords),this,OUTN1);
//        std::cerr<<"mN1c1"<<std::endl;
      }

      if(coords->mQIc != NULL){
//        std::cerr<<"mQIc1"<<std::endl;
        coords->mQIc->print(metScf,start,end,&(*coords),this,OUTQI);   
//        std::cerr<<"mQIc1"<<std::endl;
      }
      if(coords->mBKPc != NULL){
//        std::cerr<<"mBKPc1"<<std::endl;
        coords->mBKPc->print(metScf,start,end,&(*coords),this,OUTBKP);
//        std::cerr<<"mBKPc2"<<std::endl;
      }
    }
    i++;
  }
}  

//---------------------------------- Metassembly ----------------------------//
Metassembly::Metassembly(std::string ioutp,
                         std::string scf_prefix,
                         char * iDelta, 
                         char * mumcoords, //--MUMfile options
                         int iminOver,
                         int iminBas,
                         int imaxEdge, 
                         char * Rcefile, //-- CEstat files 
                         char * Qcefile,
                         int minMisOr,
                         int minSup,
                         char * RNsfile, //-- Ns files
                         char * QNsfile,
                         int iminCov, //-- metassembly
                         double iminZstat,
                         double iminDiffZ,
                         int iminLinkCov,
                         int iminNinsert,
                         double iminNinsPer)
: mum(std::string(mumcoords),std::string(iDelta),iminOver,iminBas,imaxEdge),
  Qce(std::string(Qcefile),"Q."),
  Rce(std::string(Rcefile),"R."),
  QNs(std::string(QNsfile),"Q."), 
  RNs(std::string(RNsfile),"R.")
{
  outp=ioutp;
  mscfPrefix=scf_prefix;
  metassem.clear();
  metaLengths.clear();
  Scf_Meta.clear();
  inversions.clear();
  scfsLinks.clear();

  Qce.MIN_MIS_OR=minMisOr;
  Qce.ZSTAT_THRE=iminZstat;
  Qce.MIN_SUP=minSup;
  Qce.MIN_COV=iminCov;

  Rce.MIN_MIS_OR=minMisOr;
  Rce.ZSTAT_THRE=iminZstat;
  Rce.MIN_SUP=minSup;
  Rce.MIN_COV=iminCov;

  LNs.clear();
  NAs.clear();
  NBs.clear();
  NCs.clear();
  N1s.clear();
  QIs.clear();
  BKPs.clear();

  minLinkCov=iminLinkCov;
  minCov=iminCov;
  minNinsert=iminNinsert;
  minNinsPer=iminNinsPer;
  minDiffZ=iminDiffZ;
  minZstat=iminZstat;
  minLenSDs=2;
  maxQinsLength=(Qce.mu0+5*Qce.sd0)*4;
  maxRinsLength=(Rce.mu0+5*Rce.sd0)*4;
  QmaxNinsert=15;
  QmaxNinsPer=0.02;

  ScfCount=0; 
  countInversions=0;
  countNevents=0;
  countNret=0;
  countNA=0;
  countNB=0;
  countNC=0;
  countLN=0;
  countBKP=0;
  countQI=0;
  countNrows=0;
}

