#ifndef PSEUDO_RANDOM_H
#define PSEUDO_RANDOM_H

#include <stdint.h>

/*
⭐2026-11-12
* Random is a helper class generating pseudo random numbers starting from a seed
*
* Here are some usefull abbreviations:
* LCG - Linear Congruential Generator
*
*/
class PseudoRandom {
 public:
  inline static void setSeed(unsigned long s) { seed_ = s; }
  static unsigned long lcg_rand();
  static double lcg_uniform_rand();
  const static uint64_t kRandMax = 4294967291u;
 protected:
  static uint64_t seed_;
};


#endif // RANDOM_H
