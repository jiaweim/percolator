/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"
#include "SqtSanityCheck.h"

// Class doing the sanity check for non SQT file condition
// In most sence a place holder as very little logic is build-in in this case
SanityCheck::SanityCheck() :
  initPositives_(0), pTestset_(NULL), pTrainset_(NULL), concatenatedSearch_(true),
  test_fdr_(0.01), initial_train_fdr_(0.01) {
}

SanityCheck::~SanityCheck() {
}

/**
 * 是否忽略方向异常
 */
bool SanityCheck::overRule = false;
/**
 * 初始权重文件名
 */
string SanityCheck::initWeightFN = "";
/**
 * 默认特征编号和方向，正数表示特征越大越好，负数表示越小越好，0 表示不指定方向
 * 例如 3 表示第三个特征越大越好，-3 表示使用第三个特征，数值也也小越好
 * 实际访问时要 -1 
 */
int SanityCheck::initDefaultDir = 0;
/**
 * 默认特征名称
 */
string SanityCheck::initDefaultDirName = "";
vector<double> SanityCheck::default_weights = vector<double>();

/**
 * Returns an instance of the appropriate sanity check based on information
 * contained in otherCall.
 */
SanityCheck* SanityCheck::initialize(string otherCall){
  if(initWeightFN != "" || initDefaultDirName.size() > 0) {
    return new SanityCheck();
  } else if (otherCall.find(SqtSanityCheck::fingerPrint)!= string::npos){
    return new SqtSanityCheck();
  } else {
    return new SanityCheck();
  }
}

/**
 * 特征名称转换为特征编号
 */
void SanityCheck::checkAndSetDefaultDir() {
  if (!initDefaultDir && initDefaultDirName.size() > 0) {
    int sign = 1;
    if (initDefaultDirName[0] == '-') {
      initDefaultDirName.erase(0,1);
      sign = -1;
    }
    initDefaultDir = sign * DataSet::getFeatureNames().getFeatureNumber(initDefaultDirName);
    if (initDefaultDir == 0) {
      ostringstream temp;
      temp << "ERROR: Initial direction feature name \"" << initDefaultDirName << "\" not found" << std::endl;
      throw MyException(temp.str());
    }
  }
}

/**
 * ⭐2026-11-11
 * 交叉验证方向初始化的总调用方法。
 * @param testset 测试集，每个 fold 一个
 * @param trainset 训练集，每个 fold 一个
 * @param pNorm 归一化器
 * @param w 输出：每个 fold 的权重向量
 * @param test_fdr 最终评估用的 FDR
 * @param initial_train_fdr 方向搜索用的 FDR
 */
int SanityCheck::getInitDirection(vector<Scores>& testset,
                                  vector<Scores>& trainset,
                                  const Normalizer* pNorm,
                                  vector<vector<double> > & w,
                                  double test_fdr,
                                  double initial_train_fdr) {
  pTestset_ = &testset;
  pTrainset_ = &trainset;
  test_fdr_= test_fdr;
  initial_train_fdr_ = initial_train_fdr;
  if (initWeightFN.size() > 0) {
    // ww 暂存文件中的原始权重
    vector<double> ww(FeatureNames::getNumFeatures() + 1);
    ifstream weightStream(initWeightFN.data(), ios::in);
    if (weightStream.is_open()) {
      // 读取原始权重
      readWeights(weightStream, ww);
      weightStream.close();
      assert(pNorm); // 从文件读取权重时才需要
      pNorm->normalizeweight(ww, w[0]); // 将原始空间的 ww 映射到归一化空间，写入 w[0]
      for (size_t set = 1; set < w.size(); ++set) { // 将权重复制到所有 folds
        copy(w[0].begin(), w[0].end(), w[set].begin());
      }
    } else { // 读取权重文件出错
      std::cerr << "WARNING: Could not find weights input file " << initWeightFN
                << ". Using default weights instead." << std::endl;
      getDefaultDirection(w);
    }
  } else { // 完全没有提供初始权重
    getDefaultDirection(w);
  }

  // 3. 在测试集上评估初始方向
  initPositives_ = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    // 对每个 fold,用 w[set] 给测试集打分，累加 q<test_fdr 的 target 数
    initPositives_ += (*pTestset_)[set].calcScoresAndQvals(w[set], test_fdr);
  }
  return initPositives_;
}

/**
 * 为每个交叉验证 fold 初始化 SVM 的搜索方向向量 w。方向向量决定 SVM 从哪个初始权重开始迭代训练。⭐2026-11-11
 */
void SanityCheck::getDefaultDirection(vector<vector<double> >& w) {
    
  //If I have not been given a initial direction
  if (!initDefaultDir) {
    // 3. 从数据学习
    //if the default_weights from pin.xml are not present
    if (default_weights.size() == 0) {
      // Set init direction to be the most discriminative direction
      for (size_t set = 0; set < w.size(); ++set) {
        if (VERB > 1) std::cerr << "Split " << set + 1 << ":\t";
        calcInitDirection(w[set], set);
      }
    } else {
      // 2. 次优先级，从 pin 文件读取的权重
      // I want to assign the default vector that is present in the input file
      for (size_t set = 0; set < w.size(); ++set) {
        for (size_t ix = 0; ix < w[set].size(); ix++) {
          w[set][ix] = 0;
          if (ix < default_weights.size()){
              w[set][ix] = default_weights[ix];
          }          
        }
      }
    }
  } else {
    // 1. 优先级最高，手动设置的初始权重方向
    for (size_t set = 0; set < w.size(); ++set) {
      for (size_t ix = 0; ix < w[set].size(); ix++) {
        w[set][ix] = 0;
      }
      w[set][static_cast<std::size_t>(abs(initDefaultDir) - 1)] = (initDefaultDir < 0 ? -1 : 1);
    }
  }
}

/**
 * ⭐2026-11-11
 */
void SanityCheck::calcInitDirection(vector<double>& wSet, size_t set) {
  (*pTrainset_)[set].getInitDirection(initial_train_fdr_, wSet);
}

bool SanityCheck::validateDirection(vector<vector<double> >& w) {
  if (!pTestset_) {
    cerr << "SanityCheck wrongly configured" << endl;
    return false;
  }
  int overFDR = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    overFDR += (*pTestset_)[set].calcScoresAndQvals(w[set], test_fdr_);
  }
  if (overFDR <= 0) {
    cerr << "No targets found with q<" << test_fdr_ << endl;
    resetDirection(w);
    return false;
  }
  if (initPositives_ > overFDR) {
    cerr << "Less identifications (" << overFDR << " vs " << initPositives_
        << ") after percolator processing than before processing" << endl;
    resetDirection(w);
    return false;
  }
  if (initDefaultDir) {
    for (size_t set = 0; set < w.size(); ++set) {
      if (w[set][static_cast<std::size_t>(abs(initDefaultDir) - 1)] * initDefaultDir <= 0) {
        if (VERB > 1) {
          cerr
            << "Warning, wrong sign of the weight for main scoring direction"
            << endl;
        }
        resetDirection(w);
        return false;
      }
    }
  }
  return true;
}

/**
 * 从权重文件读取 SVM 的权重向量 w：包括特征权重和 bias (截距)⭐2026-11-11
 */
void SanityCheck::readWeights(istream& weightStream, vector<double>& w) {
  std::string line;
  // Find first line with weights, but skip this and get second line containing raw features
  while (std::getline(weightStream, line)) {
    if (line[0] == '-' || (line[0] >= '0' && line[0] <= '9')) break;
  }
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream >> w[ix];
  }
  cerr << "Read weights from file" << endl;
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    cerr << w[ix] << "\t";
  }
  cerr << endl;
}

void SanityCheck::resetDirection(vector<vector<double> >& w) {
  if (!overRule) {
    cerr << "Resetting score vector, using default vector. Use --override flag to prevent this." << endl;
    getDefaultDirection(w);
  }
}

