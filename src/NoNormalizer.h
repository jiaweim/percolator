#ifndef NONORMALIZER_H_
#define NONORMALIZER_H_
// ⭐
class NoNormalizer : public Normalizer { // virtual Normalizer
 public:
  NoNormalizer();
  virtual ~NoNormalizer();
  virtual void setSet(vector<double*> & featuresV,
                      vector<double*> & rtFeaturesV, size_t numFeatures,
                      size_t numRetentionFeatures);
  virtual void updateSet(vector<double*> & featuresV, size_t offset,
                         size_t numFeatures);
  void unnormalizeweight(const vector<double>& in, vector<double>& out) const;
  void normalizeweight(const vector<double>& in, vector<double>& out) const;
};

#endif /*NONORMALIZER_H_*/
