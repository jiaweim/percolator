#ifndef DATASET_H_
#define DATASET_H_

#include <string>
#include <cassert>
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <cerrno>
#include <random>

#include "Scores.h"
#include "ResultHolder.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "FeatureNames.h"
#include "FeatureMemoryPool.h"
#include "ProteinProbEstimator.h"
#include "TabFileValidator.h"


// ⭐2026-11-12

/**
 * PIN 文件中的可选字段⭐
 */
enum OptionalField {
  SCANNR, EXPMASS, CALCMASS, RETTIME, DELTAMASS, FILENAME
};

/**
 * Dataset 表示一组同标签的 PSM。⭐
 * 所有数据集共享特征名称和特征数量。
 */
class DataSet {
 public:
 /**
  * 构造函数，默认空的
  */
  DataSet();
  virtual ~DataSet();
  
  void initFeatureTables(const unsigned int numFeatures);
  
  void inline setLabel(LabelType l) { label_ = l; }
  LabelType inline getLabel() const { return label_; }
  
  unsigned int inline getSize() const { return static_cast<unsigned int>(psms_.size()); }
    
  static FeatureNames& getFeatureNames() { return featureNames_; }
  /**
   * 清空特征名称和特征数量，通常用于读取新数据集或调用 SetHandler::reset() 时
   */
  static void resetFeatureNames() { 
    featureNames_ = FeatureNames();
    FeatureNames::resetNumFeatures();
  }
  /**
   * 返回特征数量，用于分配特征数组和遍历特征列。
   * @return 特征数量
   */
  static unsigned getNumFeatures() { return static_cast<unsigned>(featureNames_.getNumFeatures()); }
  
  bool writeTabData(std::ofstream& out);
  
  void print_10features();
  void print_features();

  void fillScores(std::vector<ScoreHolder>& scores);
  void fillFeatures(std::vector<double*>& features);
  
  void readPsm(const std::string& line, const unsigned int lineNr,
               const std::vector<OptionalField>& optionalFields, 
               FeatureMemoryPool& featurePool, std::string decoyPrefix);
  static LabelType readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool, std::string decoyPrefix);
  
  void registerPsm(PSMDescription* myPsm);
  
 protected:   
  /**
   * 所有 PSM 的指针，由 DataSet 析构函数清理
   */
  std::vector<PSMDescription*> psms_;
  /**
   * 数据集标签，Percolator 将 target 和 decoy 各自放在一个 DataSet 中
   */
  LabelType label_;
  /**
   * static，PSM 特征名称表，所有 DataSet 共享同一个特征名称表
   */
  static FeatureNames featureNames_;
  /**
   * decoy 前缀不匹配的警告是否已经触发过，只输出一次 decoy 前缀警告。
   */
  static bool decoyWarningTripped_;
};

#endif /*DATASET_H_*/
