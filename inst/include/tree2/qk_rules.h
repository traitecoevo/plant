// -*-c++-*-
#ifndef TREE_QK_RULES_H_
#define TREE_QK_RULES_H_

#include <stddef.h>

namespace quadrature {

struct QK15 {
  static const size_t n  = 8;
  static const size_t ng = 4;
  static const double xgk[8]; // abscissae of the 15-point kronrod rule
  static const double wg[4];  // weights of the 7-point gauss rule
  static const double wgk[8]; // weights of the 15-point kronrod rule
};

struct QK21 {
  static const size_t n  = 11;
  static const size_t ng = 5;
  static const double xgk[11]; // abscissae of the 21-point kronrod rule
  static const double wg[5];   // weights of the 10-point gauss rule
  static const double wgk[11]; // weights of the 21-point kronrod rule
};

struct QK31 {
  static const size_t n = 16;
  static const size_t ng = 8;
  static const double xgk[16]; // abscissae of the 31-point kronrod rule
  static const double wg[8];   // weights of the 15-point gauss rule
  static const double wgk[16]; // weights of the 31-point kronrod rule
};

struct QK41 {
  static const size_t n = 21;
  static const size_t ng = 11;
  static const double xgk[21]; // abscissae of the 41-point kronrod rule
  static const double wg[11];  // weights of the 20-point gauss rule
  static const double wgk[21]; // weights of the 41-point kronrod rule
};

struct QK51 {
  static const size_t n = 26;
  static const size_t ng = 13;
  static const double xgk[26]; // abscissae of the 51-point kronrod rule
  static const double wg[13];  // weights of the 25-point gauss rule
  static const double wgk[26]; // weights of the 51-point kronrod rule
};

struct QK61 {
  static const size_t n = 31;
  static const size_t ng = 15;
  static const double xgk[31]; // abscissae of the 61-point kronrod rule
  static const double wg[15];  // weights of the 30-point gauss rule
  static const double wgk[31]; // weights of the 61-point kronrod rule
};

}

#endif
