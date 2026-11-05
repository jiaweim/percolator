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
#include "StdvNormalizer.h"
#include "Globals.h"


StdvNormalizer::StdvNormalizer() {
}

StdvNormalizer::~StdvNormalizer() {
}

void StdvNormalizer::unnormalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  double sum = 0;
  unsigned int i = 0;
  for (; i < numFeatures; i++) {
    out[i] = in[i] / div[i];
    sum += sub[i] * in[i] / div[i];
  }
  out[i] = in[i] - sum;
}

void StdvNormalizer::normalizeweight(const std::vector<double>& in, 
    std::vector<double>& out) const {
  double sum = 0;
  size_t i = 0;
  for (; i < numFeatures; i++) {
    out[i] = in[i] * div[i];
    sum += sub[i] * in[i];
  }
  out[i] = in[i] + sum;
}

/**
 * 计算每个特征的 sub 和 div⭐
 */
void StdvNormalizer::setSet(std::vector<double*>& featuresV,
                            std::vector<double*>& rtFeaturesV, size_t nf,
                            size_t nrf) {
  numFeatures = nf;
  numRetentionFeatures = nrf;
  sub.resize(nf + nrf, 0.0);
  div.resize(nf + nrf, 0.0);


  double n = 0.0;
  double* features;
  size_t ix;
  vector<double*>::iterator it = featuresV.begin();
  for (; it != featuresV.end(); ++it) {
    features = *it;
    n++;
    for (ix = 0; ix < numFeatures; ++ix) {
      sub[ix] += features[ix];
    }
  }
  // 这里假设保留时间特征数与 PSM 数一样
  for (it = rtFeaturesV.begin(); it != rtFeaturesV.end(); ++it) {
    features = *it;
    for (ix = numFeatures; ix < numFeatures + numRetentionFeatures; ++ix) {
      sub[ix] += features[ix - numFeatures];
    }
  }
  if (VERB > 2) {
    cerr.precision(2);
    cerr << "Normalization factors" << endl << "Avg ";
  }
  for (ix = 0; ix < numFeatures + numRetentionFeatures; ++ix) {
    if (n > 0.0) {
      sub[ix] /= n;
    }
    if (VERB > 2) {
      cerr << "\t" << sub[ix];
    }
  }
  for (it = featuresV.begin(); it != featuresV.end(); ++it) {
    features = *it;
    for (ix = 0; ix < numFeatures; ++ix) {
      if (!isfinite(features[ix])) {
        cerr << "Reached strange feature with val=" << features[ix]
            << " at col=" << ix << endl;
      }
      double d = features[ix] - sub[ix];
      div[ix] += d * d;
    }
  }
  for (it = rtFeaturesV.begin(); it != rtFeaturesV.end(); ++it) {
    features = *it;
    for (ix = numFeatures; ix < numFeatures + numRetentionFeatures; ++ix) {
      if (!isfinite(features[ix-numFeatures])) {
        cerr << "Reached strange feature with val=" << features[ix
            - numFeatures] << " at col=" << ix << endl;
      }
      double d = features[ix - numFeatures] - sub[ix];
      div[ix] += d * d;
    }
  }
  if (VERB > 2) {
    cerr << endl << "Stdv";
  }
  for (ix = 0; ix < numFeatures + numRetentionFeatures; ++ix) {
    if (div[ix] <= 0 || n == 0) {
      div[ix] = 1.0;
    } else {
      div[ix] = sqrt(div[ix] / n);
    }
    if (VERB > 2) {
      cerr << "\t" << div[ix];
    }
  }
  if (VERB > 2) {
    cerr << endl;
  }
}


void StdvNormalizer::updateSet(vector<double*> & featuresV, size_t offset,
                               size_t numFeatures) {
  double n = 0.0;
  double* features;
  size_t ix;
  vector<double*>::iterator it = featuresV.begin();
  for (; it != featuresV.end(); ++it) {
    features = *it;
    n++;
    for (ix = 0; ix < numFeatures; ++ix) {
      sub[offset + ix] += features[ix];
    }
  }
  if (VERB > 2) {
    cerr.precision(2);
    cerr << "Normalization factors" << endl << "Avg ";
  }
  for (ix = 0; ix < numFeatures; ++ix) {
    if (n > 0.0) {
      sub[offset + ix] /= n;
    }
    if (VERB > 2) {
      cerr << "\t" << sub[offset + ix];
    }
  }
  for (it = featuresV.begin(); it != featuresV.end(); ++it) {
    features = *it;
    for (ix = 0; ix < numFeatures; ++ix) {
      if (!isfinite(features[ix])) {
        cerr << "Reached strange feature with val=" << features[ix]
            << " at col=" << ix << endl;
      }
      double d = features[ix] - sub[offset + ix];
      div[offset + ix] += d * d;
    }
  }
  
  if (VERB > 2) {
    cerr << endl << "Stdv";
  }
  for (ix = 0; ix < numFeatures; ++ix) {
    if (div[offset + ix] <= 0 || n == 0) {
      div[offset + ix] = 1.0;
    } else {
      div[offset + ix] = sqrt(div[offset + ix] / n);
    }
    if (VERB > 2) {
      cerr << "\t" << div[offset + ix];
    }
  }
  if (VERB > 2) {
    cerr << endl;
  }
}
