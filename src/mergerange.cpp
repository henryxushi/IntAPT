#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "rangeset.h"
#include "bedio.h"
#include <string>
using namespace std;

bool DEBUG=false;

int MIN_CLUSTER_DIST=50;

int MIN_CLUSTER_READ=5;

/* Minimum overlap rate for two clusters to merge */
double MIN_CLUSTER_OVERLAP=0.3;

/* Check if the bed record are clustered by the whole range, not the individual segments */
bool RANGE_ONLY=false;

/* Check if we further merge ranges from the same direction within a given exon range */
bool MERGE_FROM_EXON=false;

/* Cluster ID */
int NCLUSTER=0;

void processRangeSet(vector<bedrec>& rb, ofstream &outfile){
  if(rb.size()==0)return;
  
  //create rangesets
  vector<RangeSet> rs(rb.size());
  vector<int> rsdir(rb.size());
  for(int i=0;i<rb.size();i++){
    rs[i].add(rb[i].segstart,rb[i].segend);
    int cdir=0;
    if(rb[i].dir=='+')cdir=1;
    if(rb[i].dir=='-')cdir=-1;
    rsdir[i]=cdir;
  }
  //merge these range sets
  vector<RangeSet> merged;
  vector<int> mergeddir;
  mi mergeid;
  for(int ri=0;ri<rs.size();ri++){
    //check if overlapping with existing members in merged
    bool enoughoverlap=false;
    long csize;
    range_t crange=rs[ri].toSingleRange();
    if(RANGE_ONLY){
      csize=crange.second-crange.first+1; 
    }else{
      csize=rs[ri].size();
    }
    int currentdir=rsdir[ri];
    if(DEBUG)cout<<"#BED "<<ri<<", size:"<<csize<<endl;
    //for each bed record, check if it overlaps with current merged ranges?
    for(int rm=0;rm<merged.size();rm++){
      long msize;
      long ovl;
      range_t mrange=merged[rm].toSingleRange();
      if(RANGE_ONLY){
        msize=mrange.second-mrange.first+1;
        ovl=getoverlaplength(crange.first,crange.second,mrange.first,mrange.second);
      }else{
        msize=merged[rm].size();
        ovl=merged[rm].overlapLen(rs[ri]);
      }
      int rmdir=mergeddir[rm];
      if(DEBUG)cout<<"Checking rangeset "<<rm<<", overlap:"<<ovl<<", range size:"<<msize<<",currentdir:"<<currentdir<<endl;
      //have to deal with direction
      if( (ovl*1.0/msize>MIN_CLUSTER_OVERLAP || ovl*1.0/csize>MIN_CLUSTER_OVERLAP) 
         && (rmdir*currentdir>0  //this will cause all undirectional ranges unmerged
            || (rmdir==0 && currentdir==0 ))
      ){
        if(DEBUG)cout<<"Merge to rangeset "<<rm<<",dir="<<rmdir<<endl;
        enoughoverlap=true;
        merged[rm].add(rs[ri]);
        mergeid[rm].push_back(ri);
        //mergeddir[rm]=((rmdir+currentdir)>0?1:-1);
        break;
      }
    }
    if(enoughoverlap==false){
      if(DEBUG)cout<<"Add to rangeset as "<<merged.size()<<endl;
      merged.push_back(rs[ri]);
      mergeid.push_back(vi(1,ri));
      mergeddir.push_back(rsdir[ri]);
    }
  }
  
  //for the merged ranges, 
  //check all undirected ranges. If it has enough overlap with others,
  //delete them.
  //check directed ranges to see if they can merge
  vi validmerged(merged.size(),1);
  vector<range_t> exonrange; //record the range of the exons
  for(int i=0;i<merged.size();i++){
    if(mergeddir[i]==0){
      long msize;
      range_t mrange=merged[i].toSingleRange();
      //if(RANGE_ONLY){ msize=mrange.second-mrange.first+1; }else{
        msize=merged[i].size();
      //}     
      //check directed ranges one by one
      long overlapsize=0;
      for(int j=0;j<merged.size();j++){
        if(mergeddir[j]!=0){
          long ovl;
          range_t jrange=merged[j].toSingleRange();
          //if(RANGE_ONLY){   ovl=getoverlaplength(mrange.first,mrange.second,jrange.first,jrange.second);  }else{
            ovl=merged[i].overlapLen(merged[j]);
          //}
          overlapsize+=ovl;
        }
      }
      if(DEBUG)cout<<"Merged "<<i<<" is undirected. Overlap size:"<<overlapsize<<", msize:"<<msize<<".\n";
      if(overlapsize*1.0/msize>MIN_CLUSTER_OVERLAP) {
        if(DEBUG)cout<<"Delete this range.\n";
        validmerged[i]=0;
        vector<range_t> drg=merged[i].toRange();
        for(int k=0;k<drg.size();k++)exonrange.push_back(drg[k]);
      }
    }else{
      //for the directed ranges, check if they can be merged.
      if(validmerged[i]==0)continue;
      for(int j=i+1;j<merged.size();j++){
        if(validmerged[j]==0 || mergeddir[j]!=mergeddir[i])continue;
        long m2ovl=merged[i].overlapLen(merged[j]);
        if(m2ovl*1.0/merged[i].size()>MIN_CLUSTER_OVERLAP 
           || m2ovl*1.0/merged[j].size()>MIN_CLUSTER_OVERLAP){
          if(DEBUG)cout<<"Merge merged range "<<i<<" and "<<j<<", overlap:"<<m2ovl<<", i,j size: "<<merged[i].size()<<","<<merged[j].size()<<endl;
          merged[i].add(merged[j]);
          validmerged[j]=0;
          for(int k=0;k<mergeid[j].size();k++)mergeid[i].push_back(mergeid[j][k]);
        }
      }
    }
  }

  
  //Make use of the exon range (of the undirected transcripts or reference transcripts), to merge two merged ones into one
  //only merge directed ones
  if(MERGE_FROM_EXON){
    for(int i=0;i<merged.size();i++){
      if(validmerged[i]==0 || mergeddir[i]==0)continue;
      range_t mri=merged[i].toSingleRange();
      for(int j=i+1;j<merged.size();j++){  
        if(validmerged[j]==0 || mergeddir[j]!=mergeddir[i])continue;
        range_t mrj=merged[j].toSingleRange();
        //search the exon range list. First, identify the end of the 1st range and the beginning of the 2nd range
        long t1end=mri.second; long t2start=mrj.first;
        if(t1end>t2start){t1end=mrj.second;t2start=mri.first; }
        for(int k=0;k<exonrange.size();k++){
          if(exonrange[k].first<t1end && t1end<t2start && t2start<exonrange[k].second){
            //merge both range and the exon
            merged[i].add(merged[j]);
            merged[i].add(exonrange[k],true);
            //disable the second range
            validmerged[j]=0;     
          }//end if
        }//end for k
      }//end for j
    }//end for i
  }//end if

  //write to bed
  bedrec rec;
  rec.chr=rb[0].chr;
  rec.color="255,0,0";
  //cout<<"--RESULTS--\n";
  NCLUSTER++;
  for(int i=0;i<merged.size();i++){
    if(validmerged[i]==0)continue;
    //set up bed record
    if(RANGE_ONLY){
      range_t crg=merged[i].toSingleRange();
      rec.segstart.clear();rec.segend.clear();
      rec.segstart.push_back(crg.first);
      rec.segend.push_back(crg.second);
    }else{
      rec.segstart=merged[i].getStart();
      rec.segend=merged[i].getEnd();
    }
    rec.start=rec.segstart.front();
    rec.end=rec.segend.back();
    rec.nsegs=rec.segstart.size();
    switch(mergeddir[i]){
      case 1: rec.dir='+';break;
      case -1:rec.dir='-';break;
      default: rec.dir='.';
    }
    //set up score
    rec.score=0;
    stringstream nss;
    nss<<"Inst_"<<NCLUSTER<<"_"<<i<<"_m"<<mergeid[i].size();
    nss>>rec.name;
    for(int j=0;j<mergeid[i].size();j++){
      rec.score+=rb[mergeid[i][j]].score;
    }
    writeBed(outfile,rec);
  }
  //cout<<"--END RESULTS--\n";
}

void mergebed(string filename,string outputfile){
  ifstream ifs(filename.c_str());
  if(!ifs.is_open()){
    cerr<<"Error opening BED file "<<filename<<endl;
    exit(-1);
  }
  ofstream outfile(outputfile.c_str());
  int linecount=0;
  string prevchr="";
  range_t currentrange;
  vector<bedrec> rs;
  while(true){
    string oneline;
    getline(ifs,oneline);
    linecount++;
    bedrec rec;
    if(ifs.eof())break;
    if(scanBed(oneline,rec)==-1)continue;
    if(rec.score<MIN_CLUSTER_READ)continue;
    if(prevchr!=rec.chr || rec.start> currentrange.second+MIN_CLUSTER_DIST){
      //save the current rangeset
      processRangeSet(rs,outfile);
      //update the record
      currentrange.first=rec.start;
      currentrange.second=rec.end;
      rs.clear();
    }
    if(rec.end>currentrange.second){
      currentrange.second=rec.end;
    }
    rs.push_back(rec);
    prevchr=rec.chr;
  }
  if(rs.size()>0) {
    processRangeSet(rs,outfile);
    rs.clear();
  }
  ifs.close();
  outfile.close();
}