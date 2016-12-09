#include "rangeset.h"
//#include <cmath>
#include <climits>
#include <iostream>

inline long abs2(long a){
  return a>0?(a):(-a);
}

ostream& operator<<(ostream& out,const RangeSet& rs){
  //for(map<long,int>::const_iterator mit=rs.cvg.begin();mit!=rs.cvg.end(); mit++){
  //  out<<"("<<mit->first<<","<<mit->second<<") ";
  //}
  vector<range_t> vr=rs.toRange();
  for(int i=0;i<vr.size();i++){
    out<<"["<<vr[i].first<<","<<vr[i].second<<"] ";
  }
  return out;
}


/*
Add a range to range set.
It will try to merge existing ranges
*/
int RangeSet::add(range_t r, bool ismerge){
  if(cvg[r.first]<r.second) cvg[r.first]=r.second;
  ////OLD CODE
  bool removelast=false;
  if(isOverlap(r.second+1))removelast=true;
  cvg[r.first]++;
  cvg[r.second]++;
  cvg[r.second+1]=cvg[r.second+1];
  map<long,int>::iterator i1=cvg.find(r.first), i2=cvg.find(r.second);
  if(i1!=i2){
    i1++;
    while(i1!=i2){
      cvg.erase(i1++);
    };
  }
  if(removelast) cvg.erase(r.second+1);
  ismerged=false;
  //merge, if necessary
  if(ismerge)
      merge();
  return 0;
}
int RangeSet::add(const vector<range_t>& v_r){
  for(int i=0;i<v_r.size();i++){
    add(v_r[i],false);
  }
  merge();
  return 0;
}
int RangeSet::add(pos_t & start, pos_t & end){
  for(int i=0;i<start.size();i++){
    int a=start[i],b=end[i];
    add(range_t(a,b),false);
  }
  merge();
  return 0;
}
int RangeSet::add(const RangeSet& r){
  return add(r.toRange());
}
/*
Condense the range set
*/
void RangeSet::merge(){
  bool prevv=false;
  for(map<long,int>::iterator mitr=cvg.begin();mitr!=cvg.end();){
    if(prevv!=(mitr->second>0)){
      prevv=(mitr->second>0);
      mitr->second=(mitr->second>0?1:0);
      mitr++;
    }
    else{
      cvg.erase(mitr++);
    }
  }
  ismerged=true;
}
/*
Check if an existing range overlaps with current ranges
Only works if cvg are merged
*/
bool RangeSet::isOverlap(range_t r){
  map<long,int>::iterator ihow=cvg.upper_bound(r.first);
    if(ihow==cvg.end()) return false;
    if(ihow->second==0) return true;
    else{
     if(ihow->first>r.second)
       return false;
     else 
       return true;
    }
}
long RangeSet::minDist(range_t r){
  map<long,int>::iterator ihow=cvg.upper_bound(r.first);
    if(ihow==cvg.end()) return abs2(r.first+1-(cvg.rbegin()->first));
    if(ihow->second==0) return 0;
    else{
     if(ihow->first>r.second){
       long k=abs2(ihow->first-r.second);
       if(ihow!=cvg.begin()){
        ihow--;//should point to 0
        long r2=abs2(r.first+1-ihow->first);
        if(r2<k)k=r2;
      }   
      return k;        
     }
     else 
       return 0;
    }
}
bool RangeSet::isOverlap(long x){
  map<long,int>::iterator ihow=cvg.upper_bound(x);
  if (ihow==cvg.end())return false;
  if (ihow->second==0) return true;
  else return false;
}
long RangeSet::minDist(long x){
  map<long,int>::iterator ihow=cvg.upper_bound(x);
  if (ihow==cvg.end())return abs2(x+1-(cvg.rbegin()->first));
  if (ihow->second==0) return 0;
  long r=abs2(x-ihow->first);
  if(ihow!=cvg.begin()){
    ihow--;//should point to 0
    long r2=abs2(x+1-ihow->first);
    if(r2<r)r=r2;
  }   
  return r;
}

/* This range is defined by a vector of start, end points */
bool RangeSet::isOverlap(pos_t& start, pos_t& end){
  for(int i=0;i<start.size();i++){
    if(isOverlap(range_t(start[i],end[i]))) return true;
  }
  return false;
}
long RangeSet::minDist(pos_t& start, pos_t& end){
  long m=LONG_MAX;
  for(int i=0;i<start.size();i++){
    long x=minDist(range_t(start[i],end[i]));
    if(x==0) return 0;
    if(x<m)m=x;
  }
  return m;
}

bool RangeSet::isOverlap(RangeSet& r){
  map<long,int>::iterator mit;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) 
      {if (r.isOverlap(mit->first) )return true;}
    else
      {if (r.isOverlap(mit->first-1)) return true;}
  }
  for(mit=r.cvg.begin();mit!=r.cvg.end();mit++){
    if(mit->second>0) 
      {if (isOverlap(mit->first) )return true;}
    else
      {if (isOverlap(mit->first-1)) return true;}
  }
  return false;
}
long RangeSet::minDist(RangeSet& r){
  map<long,int>::iterator mit;
  long m=LONG_MAX;
  long x;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
     if(mit->second>0) {
       x=r.minDist(mit->first);
      if (x==0)return 0;
     }
     else{
      x=r.minDist(mit->first-1);
      if (x==0) return 0;
     }
     if(x<m)m=x;
  }
  for(mit=r.cvg.begin();mit!=r.cvg.end();mit++){
     if(mit->second>0){
       x=minDist(mit->first);
       if (x==0 )return 0;
     }
     else{
      x=minDist(mit->first-1);
      if (x==0) return 0;
     }
     if(x<m)m=x;
  }
  return m;
}
/* Convert to a set of ranges
Must be merged before calling it
*/
vector<range_t> RangeSet::toRange() const{
  vector<range_t> vr;
  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      vr.push_back(range_t(prev, mit->first-1));
    }
  }
  return vr;
}

/* Convert to start and end positions
  Must be merged before calling it
*/
pos_t RangeSet::getStart() const{
  pos_t pt;
  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      pt.push_back(prev);
    }
  }
  return pt;
}
pos_t RangeSet::getEnd() const{
  pos_t pt;
  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      pt.push_back( mit->first-1);
    }
  }
  return pt;
}
/* Convert to one range */
range_t RangeSet::toSingleRange() const{
  range_t t(-1,-1);
  if(cvg.size()==0)return t;
  t.first=cvg.begin()->first;
  t.second=cvg.rbegin()->first;
  return t;
}
/* 
Get the length of the range
*/
long RangeSet::size() const{
  long t=0;
  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      t+=mit->first-prev;
    }
  }
  return t;
} 

/*
Get the length of overlapping ranges
*/
long RangeSet::overlapLen(range_t r){
  long t=0;
  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      if(prev>r.second)break;
      t+=getoverlaplength(prev,mit->first-1,r.first,r.second);
    }
  }
  return t;
}
/*
Get the length of overlapping ranges
*/
long RangeSet::overlapLen(pos_t& start, pos_t& end){
  long t=0;
  int bi=0;

  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      while(bi<start.size()&& end[bi]<prev)bi++;
      if(bi>=start.size())break;
      while(bi<start.size() && start[bi] <=mit->first){
        t+=getoverlaplength(prev,mit->first-1,start[bi],end[bi]);
        bi++;
      }
    }
  }
  return t;
}
/*
Get the length of overlapping ranges
*/
long RangeSet::overlapLen(const vector<range_t>& vr){
  long t=0;
  int bi=0;

  map<long,int>::const_iterator mit;
  long prev=-1;
  for(mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->second>0) prev=mit->first;
    else{
      while(bi<vr.size()&& vr[bi].second<prev)bi++;
      if(bi>=vr.size())break;
      while(bi<vr.size() && vr[bi].first <= mit->first){
        long k=getoverlaplength(prev,mit->first-1,vr[bi].first,vr[bi].second);
        //cout<<"Checking ["<<prev<<","<<mit->first-1<<"] with ["<<vr[bi].first<<","<<vr[bi].second<<"], overlap:"<<k<<endl;
        t+=k;
        bi++;
      }
    }
  }
  return t;
}

 
/*
Get the length of overlapping ranges
*/
long RangeSet::overlapLen(RangeSet& rs){
  long t=0;
  return overlapLen(rs.toRange());
}


/* Get the sub range set within a specified range */
RangeSet RangeSet::subRangeSet(range_t r){
  RangeSet rs;
  for(map<long,int>::iterator mit=cvg.begin();mit!=cvg.end();mit++){
    if(mit->first>=r.first && mit->first<=r.second) rs.cvg[mit->first]=mit->second;
  }
  //check the first and last element
  if(rs.cvg.size()>0){
    if(rs.cvg.begin()->second==0) rs.cvg[r.first]=1;
    if(rs.cvg.rbegin()->second==1) rs.cvg[r.second]=0; 
  }else{
    rs.cvg[r.first]=1;
    rs.cvg[r.second]=0;  
  }
  return rs;
}
