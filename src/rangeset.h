/*
Range set structure.
Different to a normal range (a,b),
a range set is a set of disjoint ranges.
All ranges are 1-base inclusive.
*/

#ifndef RANGESET_H
#define RANGESET_H

#include "align.h"

/*
 * Check if two ranges [x,y] and [u,v] overlap.
 * x must be <= y and u must be <=v
 */
inline bool checkoverlap(long x, long y, long u, long v){
  if(u>y || x>v)return false;
  return true;
}

inline bool checkoverlap(range_t r1, range_t r2){
  return checkoverlap(r1.first,r1.second,r2.first,r2.second);
}
inline long max2(long x,long y){
  return x>y?x:y;
}
inline long min2(long x,long y){
  return x<y?x:y;
}
inline long getoverlaplength(long x, long y,long u, long v){
  long dist=min2(y,v)-max2(x,u)+1;
  return dist>0?dist:0;
}


class RangeSet{
  /*
  DATA, a vector of ranges
  The ranges MUST be sorted from the beginning to the end
  */
  map<long,int> cvg;
  bool ismerged;
  /*
  Function
  */
public:
  
  /*
  Add a range to range set.
  It will try to merge existing ranges
  */
  int add(range_t r, bool ismerge=true);
  int add(const vector<range_t>& v_r);
  int add(pos_t & start, pos_t & end);
  int add(const RangeSet& r);
  /*
  Condense the range set
  */
  void merge();
  /*
  Check if an existing range overlaps with current ranges
  Only works if cvg are merged
  */
  bool isOverlap(range_t r);
  bool isOverlap(long x);
  /* This range is defined by a vector of start, end points */
  bool isOverlap(pos_t& start, pos_t& end);
  bool isOverlap(RangeSet& r);
  /*
  Get the smallest distance for a range set.
  0 represents a given point(range) is within this range.
  */
  long minDist(range_t r);
  long minDist(long x);
  long minDist(pos_t& start, pos_t& end);
  long minDist(RangeSet& r);

  /* 
  Get the length of the range
  */
  long size() const;
  
  /*
  Get the length of overlapping ranges
  */
  long overlapLen(range_t r);
  long overlapLen(pos_t& start, pos_t& end);
  long overlapLen(const vector<range_t>& vr);
  long overlapLen(RangeSet& rs);
  

  /* Convert to a set of ranges
  Must be merged before calling it
  */
  vector<range_t> toRange() const;
  /* Convert to start and end positions
  Must be merged before calling it
  */
  pos_t getStart() const;
  pos_t getEnd() const;
  
  /* Convert to one range */
  range_t toSingleRange() const;


  friend ostream& operator<<(ostream& out,const RangeSet& rs);

  void clear(){cvg.clear();ismerged=false;}

  /* Get the sub range set within a specified range.
     Must use merge() in advance 
   */
  RangeSet subRangeSet(range_t r);

};

#endif
