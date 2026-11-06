#include <iostream>
#ifdef WIN32
#include <float.h>
#define isfinite _finite
#endif
#include <cmath>
#include <vector>
#include <string>

using namespace std;

#include "Normalizer.h"
#include "NoNormalizer.h"

// ⭐
NoNormalizer::NoNormalizer() {
}

NoNormalizer::~NoNormalizer() {
}

void NoNormalizer::unnormalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  unsigned int i = 0;
  for (; i < numFeatures+1; i++) {
    out[i] = in[i];
  }
}

void NoNormalizer::normalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  size_t i = 0;
  for (; i < numFeatures+1; i++) {
    out[i] = in[i];
  }
}

void NoNormalizer::setSet(std::vector<double*>& /* featuresV */,
                            std::vector<double*>& /* rtFeaturesV */, size_t nf,
                            size_t nrf) {
  numFeatures = nf;
  numRetentionFeatures = nrf;
  sub.resize(nf + nrf, 0.0);
  div.resize(nf + nrf, 0.0);
  size_t ix;
  for (ix = 0; ix < numFeatures; ++ix) {
    sub[ix] = 0.;
    div[ix] = 1.0;
  }
}


void NoNormalizer::updateSet(vector<double*> & /* featuresV */, size_t /* offset */,
                               size_t numFeatures) {
  size_t ix;
  for (ix = 0; ix < numFeatures; ++ix) {
    sub[ix] = 0.;
    div[ix] = 1.0;
  }
}
