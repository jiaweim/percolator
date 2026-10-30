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

#include "FeatureMemoryPool.h"

/**
 * 创建特定特征数量的内存池。
 * @param numFeatures 特征数量
 */
void FeatureMemoryPool::createPool(size_t numFeatures) {
  numFeatures_ = static_cast<unsigned int>(numFeatures);
  numRowsPerBlock_ = kBlockSize / numFeatures_; // 计算每个内存块能容纳多少特征行
  isInitialized_ = true;
}

/**
 * 创建一个新的连续内存块。
 */
void FeatureMemoryPool::createNewBlock() {
  double* memStart = new double[numFeatures_ * numRowsPerBlock_](); // 创建新的连续内存块
  memStarts_.push_back(memStart);
}

/**
 * 销毁内存池并释放全部内存块。
 */
void FeatureMemoryPool::destroyPool() {
  for (size_t i = 0; i < memStarts_.size(); ++i) {
    if (memStarts_.at(i) != NULL) {
      delete[] memStarts_.at(i);
      memStarts_.at(i) = NULL;
    }
  }
  isInitialized_ = false;
}

/**
 * 根据全局特征行编号 i，计算该特征行的内存地址。
 */
double* FeatureMemoryPool::addressFromIdx(unsigned int i) const {
  return memStarts_.at(i / numRowsPerBlock_) + (i % numRowsPerBlock_) * numFeatures_;
}

/**
 * 分配一个可使用的特征行，并返回首地址。
 */
double* FeatureMemoryPool::allocate() {
  if (freeRows_.size() == 0) { // 没有空闲特征行，新建特征行。
    if (initializedRows_ >= numRowsPerBlock_ * memStarts_.size()) {
      createNewBlock();
    }
    freeRows_.push_back(addressFromIdx(initializedRows_));
    initializedRows_++;
  }
  double* ret = freeRows_.back();
  freeRows_.pop_back();
  return ret;
}

/**
 * 回收一个特征行地址，使以后可以重新分配。
 */
void FeatureMemoryPool::deallocate(double* p) {
  freeRows_.push_back(p);
}
