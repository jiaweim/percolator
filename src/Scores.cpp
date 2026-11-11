#include <boost/assign.hpp>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "DataSet.h"
#include "Globals.h"
#include "MassHandler.h"
#include "Normalizer.h"
#include "PosteriorEstimator.h"
#include "Scores.h"
#include "SetHandler.h"
#include "ssl.h"
#include "IsotonicPEP.h"

/**
 * 把多个交叉验证 fold 的测试结果汇总到当前 Scores 对象中，并在合并前对各 fold
 * 的分数和统计量做必要的处理。
 * @param sv 每个 fold 的 Scores 结果集合
 * @param fdr 计算 q-value 时使用的 FDR 参数。
 * @param skipNormalizeScores 是否跳过各 fold 之间的分数归一化。
 * @param all_w 每个 fold 对应的 SVM 权重；若进行归一化，也会用于同步调整权重。
 */
void Scores::merge(std::vector<Scores>& sv,
                   double fdr,
                   bool skipNormalizeScores,
                   std::vector<std::vector<double> >& all_w) {
  scores_.clear();
  std::vector<Scores>::iterator cvBinScores = sv.begin();
  std::vector<std::vector<double> >::iterator weights = all_w.begin();
  for (; cvBinScores != sv.end(); cvBinScores++, weights++) {
    sort(cvBinScores->begin(), cvBinScores->end(), greater<ScoreHolder>());
    cvBinScores->checkSeparationAndSetPi0();
    cvBinScores->calcQvals(fdr);
    if (!skipNormalizeScores) {
      cvBinScores->normalizeScores(fdr, *weights);
    }
    copy(cvBinScores->begin(), cvBinScores->end(), back_inserter(scores_));
  }
  postMergeStep();
}

void Scores::postMergeStep() {
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder>());
  totalNumberOfDecoys_ = static_cast<unsigned int>(
      count_if(scores_.begin(), scores_.end(), mem_fn(&ScoreHolder::isDecoy)));
  totalNumberOfTargets_ = static_cast<unsigned int>(
      count_if(scores_.begin(), scores_.end(), mem_fn(&ScoreHolder::isTarget)));
  targetDecoySizeRatio_ =
      totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);
  checkSeparationAndSetPi0();
}

/**
 * 计算 SVM score ⭐2026-11-11
 */
double Scores::calcScore(const double* feat,
                         const std::vector<double>& w) const {
  std::size_t ix = FeatureNames::getNumFeatures();
  double score = w[ix];
  for (; ix--;) {
    score += feat[ix] * w[ix];
  }
  return score;
}

void Scores::scoreAndAddPSM(ScoreHolder& sh,
                            const std::vector<double>& rawWeights,
                            FeatureMemoryPool& featurePool) {
  const unsigned int numFeatures =
      static_cast<unsigned int>(FeatureNames::getNumFeatures());

  for (unsigned int j = 0; j < numFeatures; j++) {
    sh.score += sh.pPSM->features[j] * rawWeights[j];
  }
  sh.score += rawWeights[numFeatures];

  featurePool.deallocate(sh.pPSM->features);
  sh.pPSM->deleteRetentionFeatures();

  if (sh.isTarget()) {
    ++totalNumberOfTargets_;
  } else if (sh.isDecoy()) {
    ++totalNumberOfDecoys_;
  }

  if (!sh.isTarget() && !sh.isDecoy()) {
    std::cerr << "Warning: the PSM " << sh.pPSM->getId()
              << " has a label not in {1,-1} and will be ignored." << std::endl;
    PSMDescription::deletePtr(sh.pPSM);
  } else {
    scores_.push_back(sh);
  }
}

void Scores::print(LabelType label, std::ostream& os) {
 std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  os << "PSMId\t";
  if (PSMDescription::hasSpectrumFileName()) {
    os << "filename\t";
  }
  os << "score\tq-value\tposterior_error_prob\tpeptide\tproteinIds";
  if (is_output_rt_) {
    os << "\trt\n";
  } else {
    os << "\n";
  }

  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->label == label) {
      std::ostringstream out;
      scoreIt->pPSM->printProteins(out);
      ResultHolder rh(scoreIt->score, scoreIt->q, scoreIt->pep,
                      scoreIt->pPSM->getId(), scoreIt->pPSM->getFullPeptide(),
                      out.str(), scoreIt->pPSM->getSpectrumFileName());
      if (is_output_rt_) {
        rh.retentionTime = scoreIt->pPSM->getRetentionTime();
        rh.outputRT =
            true;  // Ideally this should be set once outside this loop, but
                   // ResultHolder needs to be initialized before.
      }

      os << rh << std::endl;
    }
  }
}

/**
 * 用 SetHandler 的数据填充 Scores⭐2026-11-09
 */
void Scores::populateWithPSMs(SetHandler& setHandler) {
  // 清空当前数据
  scores_.clear();
  // 读入所有 target 和 decoy 
  setHandler.populateScoresWithPSMs(scores_, LabelType::TARGET);
  setHandler.populateScoresWithPSMs(scores_, LabelType::DECOY);
  totalNumberOfTargets_ =
      static_cast<unsigned int>(setHandler.getSizeFromLabel(LabelType::TARGET));
  totalNumberOfDecoys_ =
      static_cast<unsigned int>(setHandler.getSizeFromLabel(LabelType::DECOY));
  targetDecoySizeRatio_ = (double)totalNumberOfTargets_ / totalNumberOfDecoys_;

  if (VERB > 1) {
    cerr << "Train/test set contains " << totalNumberOfTargets_
         << " positives and " << totalNumberOfDecoys_
         << " negatives, size ratio=" << targetDecoySizeRatio_
         << " and pi0=" << pi0_ << endl;
  }

  if (totalNumberOfTargets_ == 0) {
    ostringstream oss;
    oss << "Error: no target PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }

  if (totalNumberOfDecoys_ == 0) {
    ostringstream oss;
    oss << "Error: no decoy PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }

  // check for the minimum recommended number of positive and negative hits
  if (totalNumberOfTargets_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of positive samples read is too small "
                 "to perform a correct classification.\n"
              << std::endl;
  }
  if (totalNumberOfDecoys_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of negative samples read is too small "
                 "to perform a correct classification.\n"
              << std::endl;
  }
}

/**
 * ⭐2026-11-11
 * 交叉验证划分。
 * 
 * 将排序 PSM 结合划分为 xval_fold 个训练集/测试集对，并保证同一张谱图的所有 PSM 落在同一个 fold>
 * 同时处理 RESET 算法所需的 pseudo-target 转换为 nullTargetWinProb 计算。
 * 
 * Divides the PSMs from pin file into xval_fold cross-validation sets based on
 * their spectrum scan number
 * @param train vector containing the training sets of PSMs 输出参数，保存每个 fold 的训练集
 * @param test vector containing the test sets of PSMs，输出参数，保存每个 fold 的测试集
 * @param xval_fold number of folds in train and test， fold 数量
 * @param featurePool pool for storing features to optimize memory access in SVM
 * @param decoyFractionTraining fraction of decoys to use for training, will， 
 * be 1 unless RESET algorithm is used
 */
void Scores::createXvalSetsBySpectrum(std::vector<Scores>& train, // 输出 xal_fold 个训练集
                                      std::vector<Scores>& test, // 输出 xval_fold 个测试集
                                      unsigned int xval_fold, // fold 数
                                      FeatureMemoryPool& featurePool, // 特征内存池
                                      double decoyFractionTraining, // RESET: decoy 进训练集的比例
                                      unsigned int decoysPerTarget) { // RESET：每个 target 对应的 decoy 数
  // set the number of cross validation folds for train and test to xval_fold
  train.resize(xval_fold, Scores(usePi0_)); // 初始化，填充为默认 Scores 对象
  test.resize(xval_fold, Scores(usePi0_));

  // 预留大致容量
  size_t approxTestSize = scores_.size() / xval_fold + 1; // 每个测试集大于 n/k+1 个 PSM
  size_t approxTrainSize = (approxTestSize < scores_.size()) ? scores_.size() - approxTestSize : scores_.size(); // 每个训练集大约 n-approxTestSize 个PSM
  // 给每个 fold 的训练集和测试集预留容量
  for (unsigned int i = 0; i < xval_fold; ++i) {
    train[i].reserve(approxTrainSize);
    test[i].reserve(approxTestSize);
  }

  // 均分余数算法：从组后一个 fold 往前，每次分剩余PSM/剩余 PSM，保证每个 fold 大小差不多
  // remain keeps track of residual space available in each fold
  std::vector<int> remain(xval_fold); // 例如，10 分为 3 3 4,11 分为 3 4 
  // set values for remain: initially each fold is assigned (tot number of
  // scores_ / tot number of folds)
  int fold = static_cast<int>(xval_fold), ix = static_cast<int>(scores_.size());
  while (fold--) {
    remain[static_cast<std::size_t>(fold)] = ix / (fold + 1);
    ix -= remain[static_cast<std::size_t>(fold)];
  }

  // 按谱图 hash 排序，排序后，相同 specFileNr + scan 的所有 PSM 相邻
  // 如果直接按 scan 排序，scan 相邻的谱图会落入相邻的 fold，scan 号通常与色谱梯度、采集顺序相关，可能引入系统偏差
  // 训练集和测试集的谱图在时间上分布不均
  // hash 打散 scan，让 fold 分配看起来随机，但保证同谱图相邻，这是随机化与分组的这种
  std::sort(scores_.begin(), scores_.end(), OrderScanHash());

  // 没有 PSM，无法划分
  if (scores_.size() == 0) {
    ostringstream oss;
    oss << "Error: no scored PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  // 分配循环
  // put scores into the folds; choose a fold (at random) and change it only
  // when scores from a new spectra are encountered
  unsigned int previousSpectrum = scores_.begin()->pPSM->scan; // 上一个 PSM 的 scan
  size_t randIndex = PseudoRandom::lcg_rand() % xval_fold;
  for (std::vector<ScoreHolder>::iterator it = scores_.begin();
       it != scores_.end(); ++it) {
    const unsigned int curScan = (*it).pPSM->scan;
    ScoreHolder sh = (*it);

    // if current score is from a different spectra than the one encountered in
    // the previous iteration, choose new fold
    if (previousSpectrum != curScan) { // 遇到新 scan 时（谱图变化），重新随机选一个 fold
      randIndex = PseudoRandom::lcg_rand() % xval_fold;
      // allow only indexes of folds that are non-full
      while (remain[randIndex] <= 0) { // remain的作用：如果随机到的 fold 已满，就换一个
        randIndex = PseudoRandom::lcg_rand() % xval_fold;
      }
    }

    // if we use multiple folds with RESET, assign 1-decoyFractionTraining as
    // pseudo targets.
    // 满足以下条件，执行 RESET 的 pseudo-target 转换
    // RESET 算法的核心思想是让一部分 decoy "伪装"成 target 进入训练集，平衡训练集的类分布
    if (xval_fold > 1u // 多折模式
      && decoyFractionTraining < 1.0 // 不是标准模式
      && sh.label == LabelType::DECOY // 当前是 decoy
      && PseudoRandom::lcg_uniform_rand() > decoyFractionTraining) { // 以 1-decoyFractionTraining 的概率转换，大约 50% 的 decoy 转换为 pesudo-target
      // From Algorithm S3 of the percolator-RESET supplementary material
      // decoyFractionTraining - the probability of assigning a decoy to the
      // training set
      sh.label = LabelType::PSEUDO_TARGET;
    }

    // 分配 train/test
    // insert
    for (unsigned int i = 0; i < xval_fold; ++i) {
      if (i == randIndex) { // 一个进入 test
        test[i].addScoreHolder(sh);
      }
      // xval_fold-1 个进入 train
      if (i != randIndex || xval_fold == 1) {
        train[i].addScoreHolder(sh);
      }
    }
    // update number of free position for used fold
    --remain[randIndex];
    // set previous spectrum to current one for next iteration
    previousSpectrum = curScan;
  }

  // without RESET: decoysPerTarget = 1 and decoyFractionTraining = 1, then
  // factor = 1 and trainNullTargetWinProb = 1/2 (and decoyFactor = 1)
  //
  // with RESET 1 decoy: decoysPerTarget = 1 and decoyFractionTraining = 0.5,
  // then factor = 1.5 and trainNullTargetWinProb = 3/4 (and decoyFactor = 3)
  //
  // with RESET 2 decoys: decoysPerTarget = 2 and decoyFractionTraining = 0.5,
  // then factor = 2 and trainNullTargetWinProb = 2/3 (and decoyFactor = 2)
  double factor =
      (decoysPerTarget + 1) - decoyFractionTraining * decoysPerTarget;
  double trainNullTargetWinProb = factor / (1.0 + decoysPerTarget);

  // calculate ratios of target over decoy for train and test set
  for (unsigned int i = 0; i < xval_fold; ++i) {
    train[i].recalculateSizes();
    test[i].recalculateSizes();

    train[i].setNullTargetWinProb(trainNullTargetWinProb);
    test[i].setNullTargetWinProb(trainNullTargetWinProb);
  }

  if (featurePool.isInitialized()) {
    boost::unordered_map<double*, double*> movedAddresses;
    size_t idx = 0;
    for (unsigned int i = 0; i < xval_fold; ++i) {
      bool isTarget = true;
      test[i].reorderFeatureRows(featurePool, isTarget, movedAddresses, idx);
      isTarget = false;
      test[i].reorderFeatureRows(featurePool, isTarget, movedAddresses, idx);
    }
  }
}

/**
 * ⭐2026-11-11
 * 重新计算该数据集中包含的 target 和 decoy PSM 数量
 */
void Scores::recalculateSizes() {
  totalNumberOfTargets_ = 0;
  totalNumberOfDecoys_ = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      ++totalNumberOfTargets_;
    } else {
      ++totalNumberOfDecoys_;
    }
  }
  targetDecoySizeRatio_ = totalNumberOfTargets_ / (double)totalNumberOfDecoys_;
}

/**
 * ⭐2026-11-11
 * 把指定类别（target 或 decoy）的 PSM 特征数据，从池中"分散"的原始位置，原地交换到"连续"的新位置，从而在后续 SVM 训练时提升 CPU 缓存命中率。
 * 这是 C++ 版 Percolator 为 FeatureMemoryPool 设计的内存布局优化。在 Java 中特征以 double[] 直接保存在 PSMDescription 里，每个数组本身就是连续的，这个方法在 Java 版中不存在，也不需要。
 */
void Scores::reorderFeatureRows(
    FeatureMemoryPool& featurePool,
    bool isTarget,
    boost::unordered_map<double*, double*>& movedAddresses,
    size_t& idx) {
  size_t numFeatures = FeatureNames::getNumFeatures();
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() == isTarget) {
      double* newAddress =
          featurePool.addressFromIdx(static_cast<unsigned int>(idx++));
      double* oldAddress = scoreIt->pPSM->features;
      while (movedAddresses.find(oldAddress) != movedAddresses.end()) {
        oldAddress = movedAddresses[oldAddress];
      }
      std::swap_ranges(oldAddress, oldAddress + numFeatures, newAddress);
      scoreIt->pPSM->features = newAddress;
      if (oldAddress != newAddress) {
        movedAddresses[newAddress] = oldAddress;
      }
    }
  }
}

/**
 * Linear rescaling of SVM scores such that q=fdr = 0 and the median decoy = -1.
 * @param fdr fdr-cutoff to use, typically --trainFDR
 * @param weights SVM weights used to compute SVM scores
 */
void Scores::normalizeScores(double fdr, std::vector<double>& weights) {
  unsigned int medianIndex = std::max(0u, totalNumberOfDecoys_ / 2u),
               decoys = 0u;

  if (scores_.size() == 0) {
    ostringstream oss;
    oss << "Error: no scored PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  std::vector<ScoreHolder>::iterator it = scores_.begin();
  double fdrScore = it->score;
  double medianDecoyScore = fdrScore + 1.0;

  for (; it != scores_.end(); ++it) {
    if (it->q < fdr)
      fdrScore = it->score;
    if (it->isDecoy()) {
      if (++decoys == medianIndex) {
        medianDecoyScore = it->score;
        break;
      }
    }
  }

  double diff = fdrScore - medianDecoyScore;
  if (diff <= 0) {
    ostringstream oss;
    oss << "Error: median decoy score <= score at " << fdr * 100
        << "% FDR. Cannot rescale scores to merge cross validation bins, try "
        << "lowering --trainFDR.\n";
    if (NO_TERMINATE) {
      cerr << oss.str()
           << "No-terminate flag set: apply offset such that median "
           << "decoy has score -1, but skipping rescaling." << std::endl;
      diff = 1.0;
      fdrScore += 1.0;
    } else {
      throw MyException(oss.str());
    }
  }

  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score -= fdrScore;
    scoreIt->score /= diff;
  }
  Normalizer::endScoreNormalizeWeights(weights, weights, fdrScore, diff);
}

/**
 * ⭐2026-11-11
 * Calculates the SVM cost/score of each PSM and sorts them
 * and calculate q values for each PSM
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcScoresAndQvals(std::vector<double>& w,
                               double fdr,
                               bool skipDecoysPlusOne) {
  calcScores(w);
  return calcQvals(fdr, skipDecoysPlusOne);
}

/**
 * ⭐2026-11-11
 * 计算 SVM 决策值 (wx+b)，按分数降序排序，可选地打印 top10, bottom 10 的调试信息。
 * 
 * 这是 SVM 训练完成后的打分阶段：w 是训练好的权重向量，方法用它给所有 PSM 打分，为后续 q-value 
 * 计算准备有序输入。
 * 
 * Calculates the SVM cost/score of each PSM and sorts them
 * 
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return 0 on successfull execution
 */
int Scores::calcScores(std::vector<double>& w) {

  // 1. 计算每个 PSM 的分数
  std::size_t ix;
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score = calcScore(scoreIt->pPSM->features, w);
  }
  // 2. 降序排序
  sort(scores_.begin(), scores_.end(), greater<ScoreHolder>());
  if (VERB > 3) {
    if (scores_.size() >= 10) {
      cerr << "10 best scores and labels" << endl;
      for (ix = 0; ix < 10; ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
      cerr << "10 worst scores and labels" << endl;
      for (ix = scores_.size() - 10; ix < scores_.size(); ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
    } else {
      cerr << "Too few scores to display top and bottom PSMs ("
           << scores_.size() << " scores found)." << endl;
    }
  }
  return 0;
}

/**
 * ⭐2026-11-09
 */
void Scores::getScoreLabelPairs(std::vector<pair<double, bool> >& combined) {
  combined.resize(scores_.size());
  transform(scores_.begin(), scores_.end(), combined.begin(),
            mem_fn(&ScoreHolder::toPair));
}

/**
 * ⭐2026-11-10
 * Calculates the q-value for each psm in scores_: the q-value is the minimal
 * FDR of any set that includes the particular psm
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcQvals(double fdr, bool skipDecoysPlusOne) {
  std::vector<pair<double, bool> > combined;
  return calcQvals(fdr, skipDecoysPlusOne, combined);
}

/**
 * 计算 q-value ⭐2026-11-10
 */
int Scores::calcQvals(double fdr, bool skipDecoysPlusOne,
                      std::vector<pair<double, bool> >& combined) {
  assert(totalNumberOfDecoys_ + totalNumberOfTargets_ == size());

  // 从 scores_ 构建 combined 列表
  getScoreLabelPairs(combined);

  // 计算 q-values
  std::vector<double> qvals;
  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
  PosteriorEstimator::getQValues(pi0_, combined, qvals, skipDecoysPlusOne,
                                 nullTargetWinProb_);


  // set q-values and count number of positives
  std::vector<double>::const_iterator qIt = qvals.begin();
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();

  // 计算 q<fdr 的 target 数
  int numPos = 0;
  for (; qIt != qvals.end(); ++qIt, ++scoreIt) {
    scoreIt->q = *qIt;
    if (scoreIt->q < fdr && scoreIt->isTarget())
      ++numPos;
  }

  return numPos;
}

void Scores::generateNegativeTrainingSet(AlgIn& data, const double cneg) {
  std::size_t ix2 = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isDecoy()) {
      data.vals[ix2] = scoreIt->pPSM->features;
      data.Y[ix2] = -1;
      // data.C[ix2] = cneg;
      ix2++;
    }
  }
  data.negatives = static_cast<int>(ix2);
}

void Scores::generatePositiveTrainingSet(AlgIn& data,
                                         const double fdr,
                                         const double cpos,
                                         const bool trainBestPositive) {
  std::size_t ix2 = static_cast<std::size_t>(data.negatives);
  int p = 0;

  std::vector<ScoreHolder>::iterator lastUniqueIt = scores_.end();
  if (trainBestPositive) {
    std::sort(scores_.begin(), scores_.end(), OrderScanLabel());
    lastUniqueIt =
        std::unique(scores_.begin(), scores_.end(), UniqueScanLabel());
    std::sort(scores_.begin(), lastUniqueIt, greater<ScoreHolder>());
  }

  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != lastUniqueIt; ++scoreIt) {
    if (scoreIt->isTarget()) {
      if (scoreIt->q <= fdr) {
        data.vals[ix2] = scoreIt->pPSM->features;
        data.Y[ix2] = 1;
        // data.C[ix2] = cpos;
        ix2++;
        ++p;
      }
    }
  }
  data.positives = p;
  data.m = static_cast<int>(ix2);
}

void Scores::weedOutRedundant() {
  std::map<std::string, unsigned int> peptideSpecCounts;
  double specCountQvalThreshold = -1.0;
  weedOutRedundant(peptideSpecCounts, specCountQvalThreshold);
}

/**
 * Routine that sees to that only unique peptides are kept (used for analysis
 * on peptide-fdr rather than psm-fdr)
 */
void Scores::weedOutRedundant(
    std::map<std::string, unsigned int>& peptideSpecCounts,
    double specCountQvalThreshold) {
  // lexicographically order the scores_ (based on peptides names,labels and
  // scores)
  std::sort(scores_.begin(), scores_.end(), lexicOrderProb());

  /*
   * much simpler version but it does not fill up the peptide-PSM map:
   * scores_.erase(std::unique(scores_.begin(), scores_.end(), mycmp),
   * scores_.end());
   */

  std::string previousPeptide = "";
  LabelType previousLabel = LabelType::UNDEFINED;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx) {
    std::string currentPeptide = scores_[idx].pPSM->getPeptideSequence();
    LabelType currentLabel = scores_[idx].label;
    if (currentPeptide != previousPeptide || currentLabel != previousLabel) {
      // insert as a new score
      scores_[lastWrittenIdx++] = scores_[idx];
      previousPeptide = currentPeptide;
      previousLabel = currentLabel;
    }
    // append the psm
    peptidePsmMap_[scores_[lastWrittenIdx - 1].pPSM].push_back(
        scores_[idx].pPSM);
    if (specCountQvalThreshold > 0.0 &&
        scores_[idx].q < specCountQvalThreshold) {
      ++peptideSpecCounts[currentPeptide];
    }
  }
  scores_.resize(lastWrittenIdx);
  postMergeStep();
}

/**
 * Routine that sees to that only unique spectra are kept for TDC
 */
void Scores::weedOutRedundantTDC() {
  // order the scores (based on spectra id and score)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassCharge());
  scores_.erase(
      std::unique(scores_.begin(), scores_.end(), UniqueScanMassCharge()),
      scores_.end());

  /* does not actually release memory because of memory fragmentation
  double previousExpMass = 0.0;
  unsigned int previousScan = 0u;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx){
    double currentExpMass = scores_.at(idx).pPSM->expMass;
    int currentScan = scores_.at(idx).pPSM->scan;
    if (currentExpMass != previousExpMass || currentScan != previousScan) {
      // insert as a new score
      scores_.at(lastWrittenIdx++).swap(scores_.at(idx));
      previousScan = currentScan;
      previousExpMass = currentExpMass;
    } else {
      PSMDescription::deletePtr(scores_.at(idx).pPSM);
    }
  }
  scores_.resize(lastWrittenIdx);
  */
  postMergeStep();
}

/**
 * Routine that sees to that only 1 target and 1 decoy spectra are kept for
 * mix-max when using multiple hits per spectrum and separate searches
 */
void Scores::weedOutRedundantMixMax() {
  // order the scores (based on spectra id and score)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassLabelCharge());
  scores_.erase(
      std::unique(scores_.begin(), scores_.end(), UniqueScanMassLabelCharge()),
      scores_.end());

  postMergeStep();
}

/**
 * ⭐2026-11-10
 * 逐个检查每个特征，把该特征分别按“越大越好”和“越小越好”当做初始分类器，计算每种方向能筛选出多少正类，最后
 * 选择表现最好的特征和方向。
 * @param initialSelectionFdr 初筛正样本使用的 FDR 阈值
 * @param direction 输出参数，函数会在其中写入选中的特征方向
 * @return 按选中的方向和 FDR 阈值统计出的 target 数量。
 */
int Scores::getInitDirection(const double initialSelectionFdr,
                             std::vector<double>& direction) {

  // 初始化当前最佳选择
  int bestPositives = -1; // 目前找到的最大 target 数
  int bestFeature = -1; // 表现最好的特征编号
  bool lowBest = false; // 最佳方向是否越低越好

  // 计算 q-value 时跳过传统的 +1 修正
  // for determining the initial direction, the decoys+1 in the FDR estimates
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true; // 关闭 +1 平滑

  // pre-allocate once and reuse across all 2*numFeatures calcQvals calls
  std::vector<pair<double, bool> > combined;
  combined.reserve(scores_.size());

  // 遍历每个特征
  for (unsigned int featNo = 0; featNo < FeatureNames::getNumFeatures();
       featNo++) {
    for (std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
         scoreIt != scores_.end(); ++scoreIt) {
      scoreIt->score = scoreIt->pPSM->features[featNo];
    }
    // 用当前特征覆盖 scores
    // 升序
    sort(scores_.begin(), scores_.end());
    // check once in forward direction (i = 0, higher scores are better) and
    // once in backward direction (i = 1, lower scores are better)
    for (int i = 0; i < 2; i++) {
      if (i == 1) {
        reverse(scores_.begin(), scores_.end());
      }
      int positives = calcQvals(initialSelectionFdr, skipDecoysPlusOne, combined);
      if (positives > bestPositives) {
        bestPositives = positives;
        bestFeature = static_cast<int>(featNo);
        lowBest = (i == 0);
      }
    }
  }
  for (std::size_t ix = FeatureNames::getNumFeatures(); ix--;) {
    direction[ix] = 0;
  }

  if (bestPositives <= 0) {
    ostringstream oss;
    oss << "Error in the input data: cannot find an initial direction with "
        << "positive training examples. "
        << "Consider setting/raising the initial training FDR threshold "
           "(--train-fdr-initial)."
        << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str();
      std::cerr << "No-terminate flag set: setting initial direction to the "
                << "first feature and ignoring the error." << std::endl;
      bestFeature = 0;
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  }

  if (bestFeature >= 0) {
    direction[static_cast<std::size_t>(bestFeature)] = (lowBest ? -1 : 1);
  }

  if (VERB > 1) {
    cerr << "Selected feature " << bestFeature + 1
         << " as initial direction. Could separate " << bestPositives
         << " training set positives with q<" << initialSelectionFdr
         << " in that direction." << endl;
  }
  return bestPositives;
}

void Scores::checkSeparationAndSetPi0() {
  std::vector<pair<double, bool> > combined;
  getScoreLabelPairs(combined);

  std::vector<double> pvals;
  PosteriorEstimator::getPValues(combined, pvals);

  pi0_ = 1.0;
  bool tooGoodSeparation = PosteriorEstimator::checkSeparation(pvals);
  if (tooGoodSeparation) {
    ostringstream oss;
    oss << "Error in the input data: too good separation between target "
        << "and decoy PSMs.\n";
    if (NO_TERMINATE) {
      cerr << oss.str();
      if (usePi0_) {
        std::cerr
            << "No-terminate flag set: setting pi0 = 1 and ignoring error."
            << std::endl;
      } else {
        std::cerr << "No-terminate flag set: ignoring error." << std::endl;
      }
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  } else if (usePi0_) {
    pi0_ = PosteriorEstimator::estimatePi0(pvals);
  }
}

void Scores::calcPep(const bool spline, const bool interp, const bool pava) {
    if (!spline) {
        // Derive PEPs from q-values so that average PEP ≈ q at every
        // threshold.  Use PAVA when --pava-pep is set, otherwise I-spline
        // for smooth output.
        const bool use_ispline = !pava;
        std::vector<double> target_q, sc;
        for (auto& sh : scores_) {
            if (sh.isTarget()) {
                target_q.push_back(sh.q);
                sc.push_back(sh.score);
            }
        }
        InferPEP reg(use_ispline);
        auto target_pep = interp
                            ? reg.qns_to_pep(target_q, sc)
                            : reg.q_to_pep(target_q);
        // Move PEPs to scoreholders. The PEPs are only defined for target,
        // We use interpolation for decoys.
        // Add elements avoiding overflow problems if last sh is a decoy
        target_pep.push_back(1.0); target_q.push_back(1.0);
        auto it_pep = target_pep.begin();
        auto it_q = target_q.begin();
        double l_q(0.0), l_pep(0.0);
        for (auto& sh : scores_) {
            if (sh.isTarget()) {
                sh.pep = *it_pep;
                // remember last (l_) pep and q for interpolation
                l_pep = *it_pep;
                l_q = *it_q;
                it_pep++; it_q++;
            } else {
                double pep = reg.interpolate(sh.q,l_q,*it_q,l_pep,*it_pep);
                sh.pep = std::max(0.0, std::min(1.0, pep));
            }
        }
    } else {
        std::vector<pair<double, bool> > combined;
        getScoreLabelPairs(combined);

        std::vector<double> peps;
        // Logistic regression on the data
        PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, true);
        for (size_t ix = 0; ix < scores_.size(); ix++) {
            scores_[ix].pep = peps[ix];
        }
    }
}

unsigned Scores::getQvaluesBelowLevel(double level) {
  unsigned hits = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() && scoreIt->q < level) {
      hits++;
    }
  }
  return hits;
}

/**
 * 是否使用 pi0⭐2026-11-11
 */
void Scores::setUsePi0(bool usePi0) {
  usePi0_ = usePi0;
}
