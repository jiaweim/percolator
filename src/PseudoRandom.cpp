#include "PseudoRandom.h"

uint64_t PseudoRandom::seed_ = 1u;

// Park–Miller random number generator
// from wikipedia
unsigned long PseudoRandom::lcg_rand() {
  seed_ = (seed_ * 279470273u) % 4294967291u;
  return seed_;
}

// Generates a random double between 0 and 1
double PseudoRandom::lcg_uniform_rand() {
  return (double)PseudoRandom::lcg_rand() / ((double)PseudoRandom::kRandMax + (double)1);
}