#ifndef STDVNORMALIZER_H_
#define STDVNORMALIZER_H_

class StdvNormalizer : public Normalizer { // virtual Normalizer
 public:
  StdvNormalizer();
  virtual ~StdvNormalizer();
  virtual void setSet(vector<double*> & featuresV,
                      vector<double*> & rtFeaturesV, size_t numFeatures,
                      size_t numRetentionFeatures);
  virtual void updateSet(vector<double*> & featuresV, size_t offset,
                         size_t numFeatures);
  void unnormalizeweight(const vector<double>& in, vector<double>& out) const;
  void normalizeweight(const vector<double>& in, vector<double>& out) const;
};

#endif /*STDVNORMALIZER_H_*/
