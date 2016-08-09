#ifndef ALIGN_H
#define ALIGN_H


#include "commontype.h"
#include <iostream>

using namespace std;
/**
isoform (read) definitions.
*/
/*
An alignment structure in the SAM/BAM
*/
class Align{
private:
  pos_t start;
  pos_t end;
  void parsecigar();
public:
  int splicedir;//the direction of the splicing, if this is a splicing element
  //SAM fields in one record
  string qname;
  int flag;
  string rname;
  long pos;
  int mapq;
  string cigar;
  string rnext;
  long pnext;
  long plen;

  int nmismatch;
  int nhits;

  /*
  parse fields from a line of string
  Return 0 if success, -1 if any error occurs.
  */
  int parse(string oneline);
  /*
  Return the current range of the reads. FOr paired-end reads, return the range of the whole fragment.
  */
  range_t getRange(bool forcesingle=false);
 
  int getReadLen() const;
  
  /*
  Return if the alignment is paired-end alignment
  */ 
  bool isPairedEnd(){return (flag & 0x0001) & (!(flag & 0x0008));}
  
  bool isValid(){return ! (flag & 0x4);}
  
  pos_t & s(){return start;}
  pos_t & e(){return end;}
  friend ostream& operator<<(ostream& out, Align& a);
};


#endif
