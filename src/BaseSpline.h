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

#ifndef BASESPLINE_H_
#define BASESPLINE_H_

#include <assert.h>
#include "Transform.h" // 数据 log, logit 变换
#include "PackedVector.h" // 稀疏向量
#include "PackedMatrix.h" // 稀疏矩阵

/**
 * 用于平滑非参数拟合与概率/PEP估计的基础样条回归类
 */
class BaseSpline {
  public:
    BaseSpline(){};
    virtual ~BaseSpline(){};
    /**
     * 样条求值
     * @param xx 自变量值
     * @return 样条函数值
     */
    double splineEval(double xx);
    /**
     * 判断alpha搜索是否收敛。
     */
    static double convergeEpsilon;
    /**
     * 判断IRLS迭代是否收敛。
     */
    static double stepEpsilon;
    /**
     * evaluateSlope()中样条斜率得分的权重。
     */
    static double weightSlope;
    /**
     * alpha搜索时的缩放系数。
     */
    static double scaleAlpha;
    
    void roughnessPenaltyIRLS();
    void roughnessPenaltyIRLS_Old();
    void iterativeReweightedLeastSquares(double alpha);
    void predict(const vector<double>& x, vector<double>& predict);
    void setData(const vector<double>& x);
    double predict(double xx) {
      return splineEval(xx);
    }
    static void solveInPlace(PackedMatrix& mat, PackedVector& res);
  protected:
    virtual void calcPZW() {}
    virtual void initg() {
      int n = static_cast<int>(x.size());
      g = PackedVector(n);
      gnew = PackedVector(n);
      w = PackedVector(n);
      z = PackedVector(n,0.5);
      gamma = PackedVector(n-2);
    }
    virtual void limitg() {}
    virtual void limitgamma() {}
    void initiateQR();
    double crossValidation(double alpha);
    double evaluateSlope(double alpha);
    pair<double, double> alphaLinearSearch(double min_p, double max_p,
                                           double p1, double p2,
                                           double cv1, double cv2);
    double alphaLinearSearchBA(double min_p, double max_p,
                               double p1, double p2,
                               double cv1, double cv2);
    void testPerformance();
    Transform transf;

    PackedMatrix Q, Qt, R;
    PackedVector gnew, w, z, dx;
    PackedVector g, gamma;
    vector<double> x;
};

#endif /*BASESPLINE_H_*/
