/*
Parsing BED file
*/
#ifndef BEDIO_H
#define BEDIO_H

#include "commontype.h"
#include <string>
#include <ostream>
using namespace std;

/*
A BED file record
The coordinates are converted to 1-base inclusive: [start,end]
*/
struct bedrec{
  string chr;
  long start,end;
  string name;
  int score;
  char dir;
  string color;
  int nsegs;
  vl segstart;
  vl segend;
};

/*
Scan one line information of bed file, and store the information into rec.
Return 0 if success, -1 if failed
*/
int scanBed(string oneline, bedrec& rec);
/*
Write one bed record to file.
Return 0 if success, -1 if failed
*/
int writeBed(ostream& out, bedrec rec);

#endif 
