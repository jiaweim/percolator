#ifndef SCORE_HOLDER_H_
#define SCORE_HOLDER_H_

#include "LabelType.h"
#include "PSMDescription.h"

class Scores;  // forward declaration

/** 
 * 将 PSM 与其打分、统计量绑定在一起的包装类。
 * 
 * ScoreHolder is a class that provides a way to assign score value to a
 * PSMDescription and have a way to compare PSMs based on the assigned
 * score value and output them into the stream.
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 */
class ScoreHolder {
 public:
 /**
  * SVM 输出的决策值，未归一化分数
  */
  double score;
  /**
   * q-value，即该 PSM 被绑定为错误时的最小 FDR 阈值
   */
  double q;
  /**
   * PEP, 后验错误概率
   */
  double pep;
  /**
   * p-value
   */
  double p;
  /**
   * 指向底层 PSMDescription 对象的指针
   */
  PSMDescription* pPSM;
  /**
   * 标签，表示该 PSM 是目标还是诱饵
   */
  LabelType label;

  ScoreHolder()
      : score(0.0),
        q(0.0),
        pep(0.0),
        p(0.0),
        label(LabelType::UNDEFINED),
        pPSM(NULL) {}
  ScoreHolder(const double s, const LabelType l, PSMDescription* psm = NULL)
      : score(s), q(0.0), pep(0.0), p(0.0), label(l), pPSM(psm) {}
  virtual ~ScoreHolder() {}

  /**
   * (打分,isTarget)，在 java 中用 Tuple 更合适。⭐2026-11-09
   */
  std::pair<double, bool> toPair() const {
    return pair<double, bool>(score, isTarget());
  }

  /**
   * 检查该 PSM 是否为 target PSM
   */
  inline bool isTarget() const {
    return label == LabelType::TARGET || label == LabelType::PSEUDO_TARGET;
  }
  /**
   * 检查该 PSM 是否为 decoy PSM
   */
  inline bool isDecoy() const { return label == LabelType::DECOY; }
  inline PSMDescription* getPSM() const { return pPSM; }
  void printPSM(ostream& os, bool printDecoys, bool printExpMass);
  void printPepXML(ostream& os, map<char, float>& aaWeight, int index);
  void printPeptide(ostream& os,
                    bool printDecoys,
                    bool printExpMass,
                    Scores& fullset);
  std::string getCharge(std::string id);
};
/**
 * 按 PSM-ID 排序
 * 用途：输出前按 ID 排序，让结果文件稳定可读
 */
struct lessThanBaseName {
  inline bool operator()(const ScoreHolder& struct1,
                         const ScoreHolder& struct2) {
    std::string id1 = struct1.pPSM->getId();
    std::string id2 = struct2.pPSM->getId();
    return (id1 < id2);
  }
};
bool operator>(const ScoreHolder& one, const ScoreHolder& other);

bool operator<(const ScoreHolder& one, const ScoreHolder& other);

/**
 * 用于计算 PEP
 */
struct lexicOrderProb {
  /**
   * 肽段序列字典升序排序
   */
  static int compStrIt(const std::string& sv1, const std::string& sv2) {
    auto it1 = sv1.begin(), end1 = sv1.end();
    auto it2 = sv2.begin(), end2 = sv2.end();

    for (; it1 != end1 && it2 != end2; ++it1, ++it2) {
      if (*it1 < *it2)
        return 1;
      if (*it2 < *it1)
        return -1;
    }
    if (it2 != end2)
      return 1;
    if (it1 != end1)
      return -1;
    return 0;
  }

  /**
   * key: pepSeq+label(降序)+score(降序)⭐
   * 用途：计算 PEP 的预处理
   */
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    const std::string& fullSeqX = x.pPSM->getFullPeptideSequence();
    const std::string& fullSeqY = y.pPSM->getFullPeptideSequence();

    std::string peptSeqX = fullSeqX.substr(2, fullSeqX.size() - 4);
    std::string peptSeqY = fullSeqY.substr(2, fullSeqY.size() - 4);

    int peptCmp = compStrIt(peptSeqX, peptSeqY);
    return ((peptCmp == 1) || ((peptCmp == 0) && (x.label > y.label)) ||
            ((peptCmp == 0) && (x.label == y.label) && (x.score > y.score)));
  }
};

/**
 * 经典的整数哈希函数。
 * 因为某些标准库中整数的哈希就是自身，这样 hash(scan)=scan，排序后就是按 scan 排序，没有随机化效果。
 * 
 * Computes a fast hash for unsigned int using the function suggested at:
 * https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
 * Note that std::hash is sometimes implemented as the identity function,
 * defeating the purpose of random shuffling.
 */
inline unsigned int fast_uint_hash(unsigned int x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x;
}

/**
 * 用于 CV 的数据划分。⭐2026-11-09
 * Percolator 需要把 PSM 随机分配到 K 折，但如果直接按顺序分配，同一张谱图的 PSM 可能集中在同一折，导致训练集
 * 和验证集之间的信息泄露。用哈希排序后，同一谱图的 PSM 的哈希值相同，会被分到同一折，保证谱图级别的分组。
 * 
 * Orders ScoreHolders by a hash computed from the specFileNr and scan.
 * Uses the hash combination function h1 ^ (h2 << 1) suggested in
 * https://en.cppreference.com/w/cpp/utility/hash
 */
struct OrderScanHash {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    size_t hash_x = fast_uint_hash(x.pPSM->specFileNr) ^
                    (fast_uint_hash(x.pPSM->scan) << 1);
    size_t hash_y = fast_uint_hash(y.pPSM->specFileNr) ^
                    (fast_uint_hash(y.pPSM->scan) << 1);
    return hash_x < hash_y;
  }
};

/**
 * 按谱图分组，组内按分数排序，用于 TDC 去重⭐
 */
struct OrderScanMassCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    if (x.pPSM->specFileNr != y.pPSM->specFileNr)
      return x.pPSM->specFileNr < y.pPSM->specFileNr;
    if (x.pPSM->scan != y.pPSM->scan)
      return x.pPSM->scan < y.pPSM->scan;
    if (x.pPSM->expMass != y.pPSM->expMass)
      return x.pPSM->expMass < y.pPSM->expMass;
    return x.score > y.score;
  }
};

/**
 * scan, mass, label(降序), score(降序)⭐
 * label 降序：PSEUDO_TARGET(2) > TARGET(1) > UNDEFINED(0) > DECOY(-1)
 */
struct OrderScanMassLabelCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    if (x.pPSM->specFileNr != y.pPSM->specFileNr)
      return x.pPSM->specFileNr < y.pPSM->specFileNr;
    if (x.pPSM->scan != y.pPSM->scan)
      return x.pPSM->scan < y.pPSM->scan;
    if (x.pPSM->expMass != y.pPSM->expMass)
      return x.pPSM->expMass < y.pPSM->expMass;
    if (x.label != y.label)
      return x.label > y.label;
    return x.score > y.score;
  }
};

/**
 * scan, label⭐
 */
struct OrderScanLabel {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    if (x.pPSM->specFileNr != y.pPSM->specFileNr)
      return x.pPSM->specFileNr < y.pPSM->specFileNr;
    if (x.pPSM->scan != y.pPSM->scan)
      return x.pPSM->scan < y.pPSM->scan;
    return x.label > y.label;
  }
};

/**
 * 去重：scan+expMass⭐
 */
struct UniqueScanMassCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    return (x.pPSM->specFileNr == y.pPSM->specFileNr) &&
           (x.pPSM->scan == y.pPSM->scan) &&
           (x.pPSM->expMass == y.pPSM->expMass);
  }
};

/**
 * 去重：scan+label+expMass
 */
struct UniqueScanMassLabelCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    return (x.pPSM->specFileNr == y.pPSM->specFileNr) &&
           (x.pPSM->scan == y.pPSM->scan) && (x.label == y.label) &&
           (x.pPSM->expMass == y.pPSM->expMass);
  }
};

/**
 * 去重：scan+label
 */
struct UniqueScanLabel {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    return (x.pPSM->specFileNr == y.pPSM->specFileNr) &&
           (x.pPSM->scan == y.pPSM->scan) && (x.label == y.label);
  }
};

#endif /*SCORE_HOLDER_H_*/