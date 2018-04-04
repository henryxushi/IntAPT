/* 
Defining common types 
*/

#ifndef COMMONTYPE_H
#define COMMONTYPE_H

#include <vector>
#include <string>
#include <map>
using namespace std;

/**
Structures
**/
typedef pair<long,long> range_t;

inline range_t make_range_t(long a,long b){return make_pair<long,long>(a,b);}


/*
A collection of positions
*/
typedef vector<long> pos_t;

/*
A collection of pos_t
*/
typedef vector<pos_t> vpos_t;


/*
vector and matrix
*/
typedef vector<int> vi;
typedef vector<long> vl;
typedef vector<double> vd;

typedef vector<vi> mi;
typedef vector<vl> ml;
typedef vector<vd> md;


/* templates for creating vector sums and matrix sums */
template<class T>
T getVectorSum(vector<T> & t){
  T res=0;
  for(int i=0;i<t.size();i++) res+=t[i];
  return res;
}

template<class T>
vector<T> getRowSum(vector<vector<T> >& t){
  vector<T> res(t.size(),0);
  for(int i=0;i<t.size();i++)
    for(int j=0;j<t[i].size();j++)
      res[i]+=t[i][j];
  return res;
}




#endif
