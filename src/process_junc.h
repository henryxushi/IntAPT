#ifndef PROCESS_JUNC_H
#define PROCESS_JUNC_H

using namespace std;

void process_junc(vector<string> &instancelist, string &outfilename, string &regoutfilename, string &intronfileprefix);
void filter_bed(string &inputbed, string &outputbed);
void filter_2nd(string &inputinst, string &outputinst);
#endif
