#ifndef SETHANDLER_H_
#define SETHANDLER_H_

#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <map>
#include <vector>
#include <string>
#include <cctype>
#include <locale>
#include <queue>
#include <climits>

#include "ResultHolder.h"
#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "SanityCheck.h"
#include "PseudoRandom.h"
#include "FeatureMemoryPool.h"

struct PSMDescriptionPriority {
  PSMDescription* psm;
  size_t priority;
  LabelType label;
  
  bool operator<(const PSMDescriptionPriority& psmp) const {
    return (priority < psmp.priority);
  }
};

/*
* ScanId is a tuple of (ScanNr, ExpMass).
*/
typedef std::pair<int, double> ScanId;

/**
 * 该类负责处理训练集、测试集和验证集，以及读写 PSM 文件
 */
class SetHandler {
 public:
  SetHandler(unsigned int maxPSMs);
  virtual ~SetHandler();

  void push_back_dataset(DataSet* ds);

  void setDecoyPrefix(std::string& decoyPrefix) { decoyPrefix_ = decoyPrefix; }

     
  //const double* getFeatures(const int setPos, const int ixPos) const; 
  size_t getMaxPSMs() { return maxPSMs_; }
  
  // Reads in tab delimited stream and returns a SanityCheck object based on
  // the presence of default weights. Returns 0 on error, 1 on success.
  int readTab(std::istream& dataStream, SanityCheck*& pCheck);
  int readAndScoreTab(std::istream& dataStream, 
    std::vector<double>& rawWeights, Scores& allScores, SanityCheck*& pCheck);
  void addQueueToSets(std::priority_queue<PSMDescriptionPriority>& subsetPSMs,
    DataSet* targetSet, DataSet* decoySet);
  
  void writeTab(const std::string& dataFN, SanityCheck* pCheck);
  void populateScoresWithPSMs(std::vector<ScoreHolder> &scores, LabelType label);
  void normalizeFeatures(Normalizer*& pNorm);
  
  LabelType const getLabel(int setPos);
  /**
   * ⭐2026-11-09
   */
  inline int getSizeFromLabel(LabelType label) {
    return static_cast<int>(subsets_[getSubsetIndexFromLabel(label)]->getSize());
  }
  
  inline DataSet* getSubset(unsigned int ix) { return (subsets_[ix]); }
  inline DataSet* getSubsetFromLabel(LabelType label) {
    return (subsets_[getSubsetIndexFromLabel(label)]);
  }
  
  static void deletePSMPointer(PSMDescription* psm);
  
  FeatureMemoryPool& getFeaturePool() { return featurePool_; }
  
  void reset();

 protected:
 /**
  * 目标采样大小，0 表示不采用
  */
  size_t maxPSMs_;
  /**
   * 通常为 [targetSet, decoySet]
   */
  std::vector<DataSet*> subsets_;
  /**
   * 所有 Dataset 共享的特征内存池，用于管理特征向量的内存
   */
  FeatureMemoryPool featurePool_;
  /**
   * 用于校验 decoy PSM 的蛋白质前缀
   */
  std::string decoyPrefix_; // Used to determine if a psm is a decoy
  
  unsigned int getSubsetIndexFromLabel(LabelType label);
  static inline std::string &rtrim(std::string &s);
  
  int getOptionalFields(const std::string& headerLine, 
    std::vector<OptionalField>& optionalFields);
  bool isDefaultDirectionLine(const std::string& defaultDirectionLine);
  int getNumFeatures(const std::string& line, int optionalFieldCount);
  void getFeatureNames(const std::string& headerLine, int numFeatures, 
    int optionalFieldCount, FeatureNames& featureNames);
  bool getInitValues(const std::string& defaultDirectionLine, 
    int optionalFieldCount, std::vector<double>& init_values);
  ScanId getScanId(const std::string& psmLine, int& label,
    std::vector<OptionalField>& optionalFields, unsigned int lineNr);
    
  void readPSMs(std::istream& dataStream, std::string& psmLine,
    bool hasInitialValueRow, bool& separateSearches,
    std::vector<OptionalField>& optionalFields);
  void readAndScorePSMs(std::istream& dataStream, std::string& psmLine,
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields, 
    std::vector<double>& rawWeights, Scores& allScores);
};

#endif /*SETHANDLER_H_*/
