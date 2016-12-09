#include "bedio.h"
#include <sstream>
#include <algorithm>


/*
Scan one line information of bed file, and store the information into rec.
Return 0 if success, -1 if failed
*/
int scanBed(string oneline, bedrec& rec){
  stringstream ss(oneline);
  ss>>rec.chr; //1st
  ss>>rec.start; //2nd
  rec.start++; // convert 0-base to 1-base; plus 1. 
  ss>>rec.end; //3rd
  ss>>rec.name; //4th, name
  ss>>rec.score; //5th, score
  ss>>rec.dir; //6th direction
  long ta,tb;
  ss>>ta;
  ss>>tb; //7th,8th, display range
  ss>>rec.color; //9th, color
  ss>>rec.nsegs; //10th, number of segments
  if(!ss.good())return -1;
  //scanning 11th
  string lenstr,startstr;
  ss>>lenstr;
  ss>>startstr;
  if(lenstr[lenstr.size()-1]==',') lenstr=lenstr.substr(0,lenstr.size()-1);
  if(startstr[startstr.size()-1]==',') startstr=startstr.substr(0,startstr.size()-1);
  if(count(lenstr.begin(),lenstr.end(),',')!=rec.nsegs-1) return -1;
  if(count(startstr.begin(),startstr.end(),',')!=rec.nsegs-1) return -1;
  rec.segstart.clear();rec.segend.clear();
  stringstream ss2(startstr);
  stringstream ss3(lenstr);
  for(int i=0;i<rec.nsegs;i++){
    int tc;
    ss2>>tc;
    char tmp;
    if(i<rec.nsegs-1)ss2>>tmp;
    rec.segstart.push_back(tc+rec.start);
  }
  //12th
  for(int i=0;i<rec.nsegs;i++){
    int tc;
    ss3>>tc;
    char tmp;
    if(i<rec.nsegs-1)ss3>>tmp;
    rec.segend.push_back(rec.segstart[i]+tc-1);
  }
  return 0;
}

/*
Write one bed record to file.
Return 0 if success, -1 if failed.
If rec.nsegs and rec.start, rec.end are inconsistent with rec.segstart/segend, use values of rec.segstart/end 
*/
int writeBed(ostream& out, bedrec rec){
  if(rec.segstart.size() !=rec.segend.size())return -1;
  if(rec.nsegs!=rec.segstart.size())rec.nsegs=rec.segstart.size();
  if(rec.nsegs==0)return -1;
  if(rec.start!=rec.segstart.front())rec.start=rec.segstart.front();
  if(rec.end!=rec.segend.back())rec.end=rec.segend.back();
  out<<rec.chr<<"\t"; //1st
  out<<rec.start-1<<"\t"; //2nd. SInce this is 0 base, minus 1.
  out<<rec.end<<"\t"; //3rd
  out<<rec.name<<"\t"; //4th
  out<<rec.score<<"\t"; //5th, score
  out<<rec.dir<<"\t"; //6th, direction
  out<<rec.start-1<<"\t"; //7th
  out<<rec.end<<"\t"; //8th
  out<<rec.color<<"\t"; //9th
  out<<rec.nsegs<<"\t"; //10th, # of segs
  //11th,len
  for(int i=0;i<rec.segend.size();i++){
    out<<rec.segend[i]-rec.segstart[i]+1;
    if(i<rec.segend.size()-1) out<<",";
  }
  out<<"\t";
  //12th,start
  for(int i=0;i<rec.segstart.size();i++){
    out<<rec.segstart[i]-rec.start;
    if(i<rec.segstart.size()-1) out<<",";
  }
  out<<endl;
  

  
  return 0;
}

