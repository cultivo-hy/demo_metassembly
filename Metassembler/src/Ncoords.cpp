/////////////////////////////////////////////////////////////////////////
//// 
//// author: Alejandro Hernandez Wences
//// date: 21 Feb 2013
//// 
//// Get GAP sequence coordinates
////
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

using namespace std;

//-----------------Globals
vector< map<string, vector< vector< int > > > > Ns;
vector<string>fastas;
char * DirOut;
//---------------------------

void ParseArgs (int argc, char ** argv);
void PrintHelp (const char * s);
void PrintUsage (const char * s);

map<string, vector< vector<int> > > fasta2Ns(ifstream& file, ofstream& Nout){

  map<string, vector< vector<int> > >data;
  string scf;
  vector<int> coords(2,0);
  int position=0;
  char c;
  int start,end;
	
  while(file.peek()!=EOF){		
		file.get(c);
    
    if(c=='>'){
      file>>scf;
      file.ignore(1000,'\n');
      position=0;
      continue;
    }
    
		if(c != '\n') {
			position++;
			if(c=='N'){
				start=position;
				end=start;
				while(file.peek()!= EOF){
					file.get(c);
					if(c !='\n'){
						position++;
						if(c=='N'){
							end++;
						} else {
              coords[0]=start;
              coords[1]=end;
							data[scf].push_back(coords);
              Nout<<scf<<"\t"
                  <<start-1<<"\t" //bedpe based
                  <<end<<"\n";
							break;
						}
          }
				}
			}
		}
	}
  return(data);
}

int main(int argc, char ** argv){
  vector<string>::iterator f;
  ifstream fas;
  ofstream Nout;
  string out;
  size_t found;
  
  ParseArgs(argc,argv);
 
  out=DirOut;
  found=out.find_last_of('/');
  if(found == std::string::npos || found != out.size()){
    out.append("/");
  }
 
  for(f=fastas.begin();f<fastas.end();f++){
    fas.open((*f).c_str());
    //out=DirOut;
    found=(*f).find_last_of('/');
    if(found != std::string::npos){
      out.append((*f).substr(found+1));
    } else{
      out.append(*f);
    }
    out.append(".Ns");
    Nout.open(out.c_str());
    Ns.push_back(fasta2Ns(fas,Nout));     
  }
}

void ParseArgs (int argc, char ** argv)
{
  int ch, errflg = 0;
  optarg = NULL;

  while ( !errflg && ((ch = getopt (argc, argv, "hN:D:")) != EOF) )
    switch (ch)
      {
      case 'h':
        PrintHelp (argv[0]);
        exit(1);
        break;

      case 'N':
        fastas.push_back(optarg);
        break;

      case 'D':
        DirOut = optarg;
        break;

      default:
        cout<<"Error"<<std::endl;
        errflg ++;
      }

  if (errflg > 0 || fastas.size()<1 || !DirOut)
  {
    PrintUsage (argv[0]);
    std::cerr << "Try '" << argv[0] << " -h' for more information.\n";
    exit (1);
  }
}

void PrintHelp (const char * s)
{
  PrintUsage (s);
  std::cerr
    << "Generate table with Ns coordinates for each scaffold"
    << "\n\n"
    << "Required:" << std::endl
    << "-N       path of .fasta file. There should be at least one -N file specified,\n"
    << "-D       Directory name where output .Ns file will be written. \n"
    << std::endl;
}


void PrintUsage (const char * s)
{
  std::cerr
    << "\nUSAGE: " << s << "  -N in.fasta -D directory/ \n\n";
  return;
}
