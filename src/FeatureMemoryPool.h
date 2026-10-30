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
#ifndef FEATURE_MEMORY_POOL_H_
#define FEATURE_MEMORY_POOL_H_

/* Adapted from https://www.thinkmind.org/download.php?articleid=computation_tools_2012_1_10_80006 */

#include <vector>
#include <iostream>

/**
 * 二维连续内存池，为特征分配内存。
 * 在 Java 中没必要复现。⭐
 * 
 * 为 PSM 的特征数组提供连续内存池。它避免为每个 PSM 逐个执行 `new double[]`，从而减少内存碎片化和分配开销。
 */
class FeatureMemoryPool {
 private:
 /**
  * 每个内存块最多包含的 double 数量
  */
   static const unsigned int kBlockSize = 65536; // in number of doubles, e.g. 0.5MB if sizeof(double) = 8
   /**
    * 每个内存块可以容纳多少特征行
    * 每个特征行包含多少个 double
    * 已经创建并纳入池管理的特征行数量
    */
   unsigned int numRowsPerBlock_, numFeatures_, initializedRows_;
   /**
    * 每个内存块的开始地址，即 double 数组的地址
    */
   std::vector<double*> memStarts_;
   /**
    * 当前空闲特征行的地址
    */
   std::vector<double*> freeRows_;
   /**
    * 内存池是否已经初始化
    */
   bool isInitialized_;
 public:
 /**
  * 默认构造函数，未分配内存
  */
  FeatureMemoryPool() : numRowsPerBlock_(0), numFeatures_(0), 
                        initializedRows_(0), isInitialized_(false) {}

  /**
   * 析构函数，释放内存
   */
  ~FeatureMemoryPool() { destroyPool(); }

  /**
   * 创建指定特征数量的内存池
   * @param numFeatures 特征数量
   */
  void createPool(size_t numFeatures);
  void createNewBlock();
  void destroyPool();
  
  inline bool isInitialized() const { return isInitialized_; }

  double* addressFromIdx(unsigned int i) const;

  double* allocate();
  void deallocate(double* p);
};

#endif /* FEATURE_MEMORY_POOL_H_ */
