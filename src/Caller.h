#ifndef CALLER_H_
#define CALLER_H_
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Timer.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Enzyme.h"
#include "Globals.h"
#include "Normalizer.h"
#include "ProteinProbEstimator.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "SetHandler.h"

#define NO_BOOST_DATE_TIME_INLINE
#include <algorithm>
#include <boost/asio.hpp>
#include <boost/functional/hash_fwd.hpp>

/*
* Percolator 的主控制器，负责解析命令行到结果输出的全流程协调。
* 
 * Main class that starts and controls the calculations.
 * 
 * In addition to the calculation Caller also handles command line
 * initialization, parsing input parameters, handling output of the
 * calculation results.
 *
 */
class Caller {
 public:
  Caller();
  virtual ~Caller();

  /// <summary>
  /// 版本信息
  /// </summary>
  /// <returns>版本</returns>
  static string greeter();
  /// <summary>
  /// 版本 + 命令 + 参数
  /// </summary>
  /// <returns></returns>
  string extendedGreeter();
  /// <summary>
  /// 解析命令行参数
  /// </summary>
  /// <param name="argc"></param>
  /// <param name="argv"></param>
  /// <returns></returns>
  bool parseOptions(int argc, char** argv);
 
  int run();

 protected:
  /// <summary>
  /// 特征归一化
  /// </summary>
  Normalizer* pNorm_;
  /// <summary>
  /// 训练前检查：默认方向、初始权重、过拟合诊断
  /// </summary>
  SanityCheck* pCheck_;
  /// <summary>
  /// 蛋白质推断
  /// </summary>
  ProteinProbEstimator* protEstimator_;
  /// <summary>
  /// 酶切模型（picked-protein 使用）
  /// </summary>
  Enzyme* enzyme_;

  // 文件输入参数
  bool tabInput_;
  bool readStdIn_;
  std::string inputFN_;
  std::vector<std::string> inputFNs_;
  /**
   * decoy 蛋白质的前缀
   */
  std::string protEstimatorDecoyPrefix_;

  // 文件输出参数
  std::string tabOutputFN_, pepXMLOutputFN_;
  std::string weightOutputFN_;
  std::string psmResultFN_, peptideResultFN_, proteinResultFN_;
  std::string decoyPsmResultFN_, decoyPeptideResultFN_, decoyProteinResultFN_;
  bool outputRT_ = false;

  // report level parameters
  bool reportUniquePeptides_;
  bool reportPepXML_;
  bool targetDecoyCompetition_;
  bool useMixMax_;
  std::string inputSearchType_;

  // SVM/CV 参数
  double selectionFdr_, initialSelectionFdr_, testFdr_;
  unsigned int numIterations_, maxPSMs_, nestedXvalBins_, numThreads_;
  double selectedCpos_, selectedCneg_;
  bool reportEachIteration_, quickValidation_, trainBestPositive_,
      skipNormalizeScores_, analytics_, useResetAlgorithm_,
      useCompositionMatch_, useIrlsPep_, useInterpolatingPep_, usePavaPep_;

  // reporting parameters
  std::string call_;

  /// <summary>
  /// 计时器
  /// </summary>
  Timer timer;

  std::istream& getDataInStream(std::ifstream& fileStream);
  bool loadAndNormalizeData(std::istream& dataStream,
                            SetHandler& setHandler,
                            Scores& allScores);
  void calcAndOutputResult(Scores& allScores);

  void calculatePSMProb(Scores& allScores, bool uniquePeptideRun);
  void writeResults(Scores& allScores, bool unique, bool writeOutput);
  void calculateProteinProbabilities(Scores& allScores);
  void checkIsWritable(const std::string& filePath);
};

#endif /*CALLER_H_*/
