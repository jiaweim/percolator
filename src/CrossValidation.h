#ifndef CROSSVALIDATION_H_
#define CROSSVALIDATION_H_

#include <iomanip>
#include <iostream>
#include <vector>

#include "DataSet.h"
#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "Globals.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "ssl.h"

struct CandidateCposCfrac {
  double cpos;
  double cfrac;
  unsigned int set;
  int nestedSet;
  vector<double> ww;
  int tp;
};

/**
 * 负责：
 * 1. 将数据划分成多个 fold。
 * 2. 在训练集中生成正负类。
 * 3. 使用 nested cross-validation 选择正负类别权重。
 * 4. 迭代重新选择训练正样本。
 * 5. 合并每个 fold 的测试分数。
 * 6. 输出原始空间或归一化空间权重。
 */
class CrossValidation {
 public:
  CrossValidation(bool quickValidation,
                  bool reportPerformanceEachIteration,
                  double testFdr,
                  double selectionFdr,
                  double initialSelectionFdr,
                  double selectedCpos,
                  double selectedCneg,
                  unsigned int niter,
                  bool usePi0,
                  unsigned int nestedXvalBins,
                  bool trainBestPositive,
                  unsigned int numThreads,
                  bool skipNormalizeScores,
                  double decoyFractionTraining,
                  unsigned int numFolds);

  ~CrossValidation();

  /**
   * 准备每一轮训练的输入
   * @param fullset
   */
  int preIterationSetup(Scores& fullset,
                        SanityCheck* pCheck,
                        const Normalizer* pNorm,
                        FeatureMemoryPool& featurePool);

  /**
   * 主训练循环
   */
  void train(const Normalizer* pNorm);
  
  /**
   * 训练后处理
   */
  void postIterationProcessing(Scores& fullset, SanityCheck* pCheck);

  /**
   * 输出最终权重
   */
  void printAllWeights(ostream& weightStream, const Normalizer* pNorm);

  void getAvgWeights(std::vector<double>& weights, const Normalizer* pNorm);

  void inline setSelectedCpos(double cpos) { selectedCpos_ = cpos; }
  double inline getSelectedCpos() { return selectedCpos_; }
  void inline setSelectedCneg(double cneg) { selectedCneg_ = cneg; }
  double inline getSelectedCneg() { return selectedCneg_; }

  void inline setSelectionFdr(double fdr) { selectionFdr_ = fdr; }
  double inline getSelectionFdr() { return selectionFdr_; }
  void inline setTestFdr(double fdr) { testFdr_ = fdr; }

  void inline setNiter(unsigned int n) { niter_ = n; }
  unsigned int inline getNiter() { return niter_; }
  void inline setQuickValidation(bool on) { quickValidation_ = on; }
  void inline setReportPerformanceEachIteration(bool on) {
    reportPerformanceEachIteration_ = on;
  }

 protected:
  std::vector<AlgIn*> svmInputs_;
  /**
   * 每个 fold 学到的 SVM 权重
   */
  std::vector<std::vector<double> > weights_;  // svm weights for each fold
  /**
   * 保存每个 fold / nested fold 下的 (Cpos, Cfrac) 参数组合及对应权重
   */
  std::vector<CandidateCposCfrac>
      classWeightsPerFold_;  // cpos, cneg pairs to train for each nested CV
                             // fold
  /**
   * 控制是否启用快速收敛判断。
   * 例如连续两轮改进很小，如 <= 1%，则提前停止训练，避免无意义迭代。
   */
  bool quickValidation_;
  bool usePi0_;
  /**
   * 是否在每轮结束后报告当前验证集性能。
   */
  bool reportPerformanceEachIteration_;

  unsigned int numThreads_;

  /**
   * 测试阶段的 FDR
   */
  double testFdr_;       // fdr used for cross validation performance measuring
  /**
   * 选择 positive 训练样本时用的 FDR
   */
  double selectionFdr_;  // fdr used for determining positive training set
  /**
   * 第一轮训练时用的 FDR
   */
  double initialSelectionFdr_;  // fdr used for determining positive training
                                // set in first iteration
  /**
   * 用户手动指定的正负类惩罚参数
   */
  double selectedCpos_;  // soft margin parameter for positive training set
  double selectedCneg_;  // soft margin parameter for negative training set

  /**
   * 最大迭代轮数
   */
  unsigned int niter_;
  /**
   * 嵌套交叉验证的 fold 数
   */
  unsigned int nestedXvalBins_;

  bool trainBestPositive_;
  bool skipNormalizeScores_;

  double decoyFractionTraining_;
  /**
   * 快速停止条件的阈值，值是 0.01
   */
  const static double requiredIncreaseOver2Iterations_;

  /**
   * 外层交叉验证的 fold 数
   */
  unsigned int numFolds_;  // number of folds for cross validation
  /**
   * 每个 fold 的训练集和测试集分数
   */
  std::vector<Scores> trainScores_, testScores_;
  std::vector<double> candidatesCpos_, candidatesCfrac_;

  void initializeGridSearch(double targetDecoySizeRatio);
  void trainCpCnPair(CandidateCposCfrac& cpCnFold,
                     options& pOptions,
                     AlgIn* svmInput);

  int mergeCpCnPairs(double selectionFdr,
                     options& pOptions,
                     std::vector<std::vector<Scores> >& nestedTestScoresVec,
                     const vector<double>& cpos_vec,
                     const vector<double>& cfrac_vec);
  int doStep(const Normalizer* pNorm, double selectionFdr);

  void printSetWeights(ostream& weightStream, unsigned int set);
  void printRawSetWeights(ostream& weightStream,
                          unsigned int set,
                          const Normalizer* pNorm);

  void printAllWeightsColumns(ostream& weightStream);
  void printAllRawWeightsColumns(ostream& weightStream,
                                 const Normalizer* pNorm);
  void printAllWeightsColumns(
      std::vector<std::vector<double> > weightMatrix,
      ostream& weightStream);
};

#endif /*CROSSVALIDATION_H_*/
