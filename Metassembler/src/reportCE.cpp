///////////////////////////////////////////////////////////////////////////////
//////// 
//////// author: Alejandro Hernandez Wences
//////// date: 21 Feb 2013
//////// 
//////// Read .mateAn file and report CEstatistic
////////
//////// see CEstat.hh
///////////////////////////////////////////////////////////////////////////////////

#include "CEstat.hh"
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace std;

//-------------------------------------- Globals
vector<string> CEfiles;
vector<string> COfiles;

//-------------------------------------- Functions

void ParseArgs (int argc, char ** argv);
void PrintHelp (const char * s);
void PrintUsage (const char * s);

int main (int argc, char**argv) 
{

  ParseArgs(argc,argv);

  vector<mateAn*>CEs;

  vector<string>::iterator cefile;
  for(cefile=CEfiles.begin(); cefile < CEfiles.end(); cefile++){
    cout<<"---- Loading "<<*cefile<<std::endl
        <<"---------------------------------"<<std::endl;
    mateAn* ce=new mateAn(*cefile,"");
    CEs.push_back(ce);
  }

  vector<string>::iterator cofile;
  int cef=0;
  for(cofile=COfiles.begin(); cofile < COfiles.end(); cofile++){

    ofstream out (string(*cofile).append(".cestat").c_str());
    ifstream in (cofile->c_str());    
    string line;
    vector<string>words;
    while(in.good() && in.peek() != EOF){
      words.clear();
      getline(in,line);
      out << line ;
      istringstream isl (line);  
      string sub;  
      while(isl && isl.peek() != EOF){
        isl >> sub;
        words.push_back(sub);
      }
      if(words.size() == 2 ){
        bool isInteger=true;
        for(int l=0; l < words[1].size(); l++){
          isInteger= isInteger && isdigit(words[1][l]);
        }
        if(isInteger){
          int found=CEs[0]->search(words[0],atoi(words[1].c_str()));
          if(found >= 0){
            mateAnRow m=CEs[0]->mateAns[words[0]][found];
            out <<"\t";
            m.print(out,words[0],CEs[0]->mu0,CEs[0]->sd0);
          } else {
            out << "\tNA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"
                << "NA\t"<<"\n";
          }
        } else {
          out << "\n";
        }
      } else {
        out << "\n";
      }
    }
  }
}

void ParseArgs (int argc, char ** argv)
{

  int ch, errflg =0;
  optarg = NULL;

  if(argc <= 1){
    PrintHelp(argv[0]);
    exit(EXIT_FAILURE);
  }

  while(!errflg && (( ch = getopt (argc, argv, "hC:P:")) != EOF) )
    switch(ch){

      case 'h':
        PrintHelp(argv[0]);
        exit(EXIT_FAILURE);
        break;
    
      case 'C':
        CEfiles.push_back(string(optarg));
        break;

      case 'P':
        COfiles.push_back(string(optarg));        
        break;
 
      default:
        cerr << "ERROR" << std::endl;
        errflg++;
    }

  if(errflg > 0      ||
     CEfiles.empty() ||
     COfiles.empty() 
    ){
    PrintUsage(argv[0]);
    cerr << "Try " << argv[0] << " -h for more information\n";
    exit(EXIT_FAILURE);
  }
}

void PrintUsage(const char * s)
{

  cerr << "\nUSAGE: " << s << " -P INpos1 [-P INpos2 ...] -C in.mateAn1\n\n";
  cerr <<"For now only supports one .mateAn file"<<std::endl<<std::endl;

}

void PrintHelp(const char * s)
{

  PrintUsage(s);

  cerr << "Find and report the compression-expansion (CE) statistic for the positions given\n"
       << "in the -P <in.file> files. Output will be in.file.zstat"<<std::endl<<std::endl;

  cerr << "-P   in.pos file. Multiple -P are allowed."<<std::endl
       << "-C   in.mateAn file."<<std::endl;

}

