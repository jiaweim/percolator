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
#ifndef FEATURENAMES_H_
#define FEATURENAMES_H_

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <cassert>
#include <cctype>
#include <iterator>

using namespace std;

/**
 * FeatureNames 类用于管理特征名称及其相关信息。⭐
 */
class FeatureNames {
  public:
    FeatureNames();
    virtual ~FeatureNames();

    std::string getFeatureNames();
    /**
     * 返回指定索引的特征名称
     * @param index 特征索引
     */
    inline std::string getFeatureName(unsigned int index) { return featureNames.at(index); }
    /**
     * 返回特征数量
     */
    static inline size_t getNumFeatures() {
      return numFeatures;
    }
    /**
     * 设置全局特征数量，只在第一次调用时生效，之后的调用会被忽略
     * @param nf 特征数量
     */
    static inline void setNumFeatures(size_t nf) {
      if (!numFeatures) {
        numFeatures = nf;
      }
    }
    /**
     * 重置全局特征数量为 0，通常用于读取新数据集或调用 SetHandler::reset() 时
     */
    static inline void resetNumFeatures() {
      numFeatures = 0;
    }

    /**
     * 初始化特征元数据并设置全局特征数量。
     */
    void initFeatures();
    
    /**
     * 将特征名称插入到 featureNames 向量中，如果已经存在则不插入。
     * @param featureName 特征名称
     */
    void insertFeature(const string& featureName) {
      if(std::find(featureNames.begin(), featureNames.end(), featureName)==featureNames.end())
	      featureNames.push_back(featureName);
    }

    /**
     * 根据特征名称返回其索引(1-based)，如果不存在则返回 0。
     * @param featureName 特征名称，不区分大小写
     * @return 特征索引
     */
    int getFeatureNumber(const string& featureName);

    /**
     * 返回数据集中检测到的最小电荷状态
     */
    int getMinCharge() {
      return minCharge;
    }
    /**
     * 返回数据集中检测到的最大电荷状态
     */
    int getMaxCharge() {
      return maxCharge;
    }

    /**
     * 为扩展功能保留的函数。属于遗留接口，可忽略。
     */
    int getChargeFeatNum() {
      assert(1==0);
      return chargeFeatNum;
    }
    int getEnzFeatNum() {
      assert(1==0);
      return enzFeatNum;
    }
    int getNumSPFeatNum() {
      assert(1==0);
      return numSPFeatNum;
    }
    int getPtmFeatNum() {
      assert(1==0);
      return ptmFeatNum;
    }
    int getIntraSetFeatNum() {
      assert(1==0);
      return intraSetFeatNum;
    }
    int getQuadraticFeatNum() {
      assert(1==0);
      return quadraticFeatNum;
    }

  protected:
  
    /**
     * 特征名称列表
     */
    vector<string> featureNames;
    /**
     * 全局特征数量，所有数据集共享
     */
    static size_t numFeatures;
    /**
     * 数据集中检测到的最小和最大电荷状态
     */
    int minCharge, maxCharge;
    /**
     * 为扩展功能保留的特征索引。属于遗留接口，可忽略。
     */
    int chargeFeatNum, enzFeatNum, numSPFeatNum, ptmFeatNum,
        intraSetFeatNum, quadraticFeatNum;
};

#endif /*FEATURENAMES_H_*/
