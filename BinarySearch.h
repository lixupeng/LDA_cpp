//
// Created by Xupeng Li on 04/04/2017.
//

#ifndef LDA_BINARYSEARCH_H
#define LDA_BINARYSEARCH_H

class BinarySearch {
public:
  /*! search the first x which x > u from [start,end] of p */
  static int binarySearch(float* p, float u, int start, int end) {
    int pstart = start, pend = end;

    while (pstart < pend) {
      if (pstart + 1 == pend) {
        if (p[pstart] > u)
          return pstart;
        else if (p[end] > u)
          return pend;
        else
          return -1;
      }

      int mid = (pstart + pend) / 2;
      double value = p[mid];
      if (value == u) {
        return mid + 1;
      }
      if (value < u) {
        pstart = mid + 1;
      } else {
        pend = mid;
      }
    }
    return pstart;
  }

  static int binarySearch(double* p, double u, int start, int end) {
    int pstart = start, pend = end;

    while (pstart < pend) {
      if (pstart + 1 == pend) {
        if (p[pstart] > u)
          return pstart;
        else if (p[end] > u)
          return pend;
        else
          return -1;
      }

      int mid = (pstart + pend) / 2;
      double value = p[mid];
      if (value == u) {
        return mid + 1;
      }
      if (value < u) {
        pstart = mid + 1;
      } else {
        pend = mid;
      }
    }
    return pstart;
  }
};

#endif //LDA_BINARYSEARCH_H
