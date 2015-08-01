#ifndef BCTREE
#define BCTREE

#include<cstdio>
#include<cstring>
#include<cassert>

#include "MultireadCSEMPosPair.h"

static long tinsert = 0;
static long tcount = 0;

class bcTree {
 public:
  bcTree() { n = 0; axis = NULL; sum = NULL; }
  bcTree(int, PosPair*);
  ~bcTree();
  inline void clear();
  inline double count(int);
  inline void insert(int, double);
  inline void del(int, double);
  inline void showTree();
 private:
  int n; // n elements;
  int *axis; // the axis of elements
  double *sum; // the sum of corresponding substree;
};

bcTree::bcTree(int len, PosPair* arr) {
  n = len;
  axis = new int[n];
  sum = new double[n];
  for (int i = 0; i < n; i++) axis[i] = arr[i].pos;
  memset(sum, 0, sizeof(double) * n);
}

bcTree::~bcTree() {
  if (n > 0) {
    delete[] axis;
    delete[] sum;
  }
  n = 0;
}

inline void bcTree::clear() {
  if (n > 0) {
    memset(sum, 0, sizeof(double) * n);
  }
}

inline double bcTree::count(int pos) {
  int l, r, mid;
  double res = 0.0;

  l = 0; r = n - 1;
  while (l <= r) {

    tcount++;

    mid = (l + r) >> 1;
    if (axis[mid] == pos) {
      res += sum[mid];
      break;
    }
    else if (axis[mid] < pos) {
      res += sum[mid];
      l = mid + 1;
    }
    else {
      r = mid - 1;
    }
  }

  return res;
}

inline void bcTree::insert(int pos, double val) {
  int l, r, mid;
  
  l = 0; r = n - 1;
  while (l <= r) {

    tinsert++;

    mid = (l + r) >> 1;
    
    if (axis[mid] == pos) {
      sum[mid] += val;
      return;
    }
    else if (axis[mid] < pos) {
      l = mid + 1;
    }
    else {
      sum[mid] += val;
      r = mid - 1;
    }
  }
  assert(false);
}

inline void bcTree::del(int pos, double val) {
  int l, r, mid;
  
  l = 0; r = n - 1;
  while (l <= r) {
    mid = (l + r) >> 1;
    if (axis[mid] == pos) {
      sum[mid] -= val;
      return;
    }
    else if (axis[mid] < pos) {
      l = mid + 1;
    }
    else {
      sum[mid] -= val;
      r = mid - 1;
    }
  }
  assert(false);
}

inline void bcTree::showTree() {
  if (n == 0) return;
  for (int i = 0; i < n; i++) 
    printf("%.6lf\t", sum[i]);
  printf("\n");
}

#endif
