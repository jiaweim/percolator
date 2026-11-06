# 总流程

## 总流程概览

```text
读取 PIN
  -> 创建 target/decoy DataSet
  -> 可选 reservoir sampling
  -> 全局拟合并归一化特征
  -> 按 spectrum/scan 分成训练/测试 fold
  -> 确定初始方向
  -> grid search Cpos/Cneg
  -> 多轮迭代训练
  -> 合并 fold 分数
  -> 校正分数范围
  -> 计算 q-value、PEP/输出
```

## 1. 读取 PIN 文件

入口：

- `Caller::readInputFile`
- `SetHandler::readTab()`
- `SetHandler::readAndScoreTab()`

调用链：

```cpp
SetHandler::readTab(dataStream, pCheck_)
```

相关实现位于：

- `SetHandler::readTab`
- `SetHandler::readAndScoreTab`

执行内容：

1. 解析 PIN header。
2. 识别特征列、`Label`、`ScanNr`、Peptide、Proteins等列。
3. 为每个 PSM 建立 `PSMDescription`。
4. 根据 `Label` 创建：
   - target PSM
   - decoy PSM
5. 创建 `SanityCheck`。
6. 将 target 和 decoy 放入两个 `DataSet`。

核心方法：

```cpp
DataSet::readPsm()
SetHandler::readPSMs()
DataSet::registerPsm()
```

其中 `readPSMs()` 负责创建目标和 decoy 数据集，参见 `SetHandler::readPSMs`。

---

## 2. 可选 reservoir sampling

如果用户指定了 `-N <maxPSMs>`，Percolator 不一定立即使用全部数据训练。

相关入口位于 `Caller.cpp`。

方法：

```cpp
SetHandler::readAndScoreTab()
SetHandler::addQueueToSets()
SetHandler::samplePSMs()
```

行为：

```text
若 PSM 数量 <= N：使用全部数据
若 PSM 数量 > N：只保留最多 N 个 PSM用于训练
```

如果指定了 `-N`，最终会重新读取全部 PIN：

```cpp
setHandler.readAndScoreTab(fileStream, rawWeights, allScores, pCheck_)
```

因此 `-N` 通常表示：

- 用一个子集训练
- 用得到的平均模型对原始完整数据集重新打分

相关代码位于 `Caller.cpp`。

---

## 3. 收集所有特征并拟合归一化器

入口：

```cpp
setHandler.normalizeFeatures(pNorm_);
```

相关实现位于 `SetHandler::normalizeFeatures`。

流程：

```cpp
for (DataSet* subset : subsets_)
    subset->fillFeatures(featuresV);

pNorm = Normalizer::getNormalizer();
pNorm->setSet(featuresV, ..., numFeatures, 0);
pNorm->normalizeSet(featuresV, rtFeaturesV);
```

相关方法：

- `SetHandler::normalizeFeatures()`
- `DataSet::fillFeatures()`
- `Normalizer::getNormalizer()`
- `Normalizer::setSet()`
- `Normalizer::normalizeSet()`
- `Normalizer::normalize()`

归一化参数由所有 target 和 decoy 合并计算，不按标签分别拟合：

| Normalizer | `sub` | `div` |
|---|---|---|
| UNI | 每个特征最小值 | 最大值减最小值 |
| STDV | 每个特征均值 | 每个特征总体标准差 |
| NONORM | 不变 | 1 |

实现位于：

- `UniNormalizer::setSet`
- `StdvNormalizer::setSet`

归一化公式：

```cpp
normalized = (original - sub[index]) / div[index]
```

所有特征原地归一化：

```cpp
Normalizer::normalize(features, features, offset, numFeatures)
```

---

## 4. 生成 Scores 和 spectrum-based 交叉验证 fold

进入交叉验证前：

```cpp
CrossValidation crossValidation(...);

crossValidation.preIterationSetup(
    allScores, pCheck_, pNorm_, setHandler.getFeaturePool());
```

入口位于 `Caller.cpp`。

主要方法：

```cpp
CrossValidation::preIterationSetup()
Scores::createXvalSetsBySpectrum()
Scores::recalculateSizes()
```

`preIterationSetup()` 在 `CrossValidation.cpp` 中执行：

1. 为多个 fold 和 nested fold 建立 `AlgIn`。
2. 按 `ScanNr` 将所有 PSM 分配到训练/测试集合。
3. 初始化每个 fold 的权重向量。
4. 确定初始方向。
5. 初始化正则化参数候选值。

默认一般是 3-fold。划分方法必须按 spectrum 的 `scan` 分组，不能让同一张谱图同时进入训练集和测试集。

相关实现位于 `Scores::createXvalSetsBySpectrum`。

---

## 5. 确定初始方向和初始权重

`SanityCheck::getInitDirection()` 尝试确定一个初始 score vector。

优先级：

```text
已有 weight 文件
  -> DEFAULT_DIRECTION
  -> 自动选择最具区分性的单特征
  -> 输入 PIN 的 default weights
```

相关实现位于 `SanityCheck.cpp`。

使用的方法：

- `SanityCheck::getInitDirection()`
- `SanityCheck::getDefaultDirection()`
- `SanityCheck::calcInitDirection()`
- `Scores::getInitDirection()`
- `SanityCheck::readWeights()`

最具区分性的特征选择逻辑位于 `Scores::getInitDirection`。

它分别尝试每个特征：

```text
feature
-feature
```

并选择在指定初始 FDR 下产生最多 target 的特征方向。

如果模型权重不合理，`SanityCheck::validateDirection()` 会执行方向验证：

```cpp
if (!validateDirection(weights)) {
    resetDirection(weights);
}
```

否则可能出现方向导致 target q-value 下降异常。

---

## 6. 初始化 Cpos/Cneg 网格搜索

`CrossValidation::initializeGridSearch()` 根据 target/decoy 比例生成候选参数。

涉及参数：

- `cpos`
- `cneg`，通常间接表示为 `cfrac`
- nested CV
- 目标/decoy 类别权重

实现位于 `CrossValidation::initializeGridSearch`。

候选值一般不会固定为单一 `C=1`，而会在多个候选中选择每个 fold 的最优组合。

---

## 7. 每一轮迭代训练

入口：

```cpp
crossValidation.train(pNorm_);
```

实现位于 `CrossValidation::train`。

每轮执行：

```text
使用当前权重计算所有 train/test PSM 分数
  -> 更新 target/decoy q-value
  -> 选择高置信 target
  -> 保留 decoy
  -> 在 nested fold 中训练 SVM
  -> 选择最优 Cpos/Cneg
  -> 更新三个主 fold 的权重
```

核心方法：

- `CrossValidation::doStep()`
- `Scores::calcScoresAndQvals()`
- `Scores::generateNegativeTrainingSet()`
- `Scores::generatePositiveTrainingSet()`
- `Scores::createXvalSetsBySpectrum()`
- `CrossValidation::trainCpCnPair()`
- `CrossValidation::mergeCpCnPairs()`

正负训练集生成：

```cpp
Scores::generateNegativeTrainingSet()
```

该方法将全部 decoy 标记为：

```cpp
label = -1
```

实现位于 `Scores.cpp`。

```cpp
Scores::generatePositiveTrainingSet()
```

该方法根据当前 q-value 和 `selectionFDR` 选择高置信 target：

```text
q-value < selectionFDR
label = +1
```

实现位于 `Scores.cpp`。

默认情况下：

```text
目标训练集 = 当前模型认为可靠的 target
负类训练集 = decoy
```

这使 Percolator 成为逐轮更新的半监督算法。

---

## 8. Nested CV 选择正则化参数

每个主 CV fold 内还会执行 nested CV。

步骤：

```cpp
trainScores_[set].createXvalSetsBySpectrum(nestedTrainScores,
                                            nestedTestScores,
                                            nestedXvalBins_...)
```

然后：

1. 对每个 nested fold 使用候选 `cpos/cfrac`。
2. 用 `Scores::generatePositiveTrainingSet()` 生成正样本。
3. 用 `Scores::generateNegativeTrainingSet()` 生成负样本。
4. 调用 `CrossValidation::trainCpCnPair()`。
5. 在 nested test set 上计算 q-value。
6. 选择产生最多 target 的 `(cpos,cfrac)`。
7. 将选中的参数应用到对应主 fold。

核心实现位于：

- `CrossValidation::trainCpCnPair`
- `CrossValidation::mergeCpCnPairs`

---

## 9. 收敛判断并停止迭代

`CrossValidation::train()` 默认最多执行一定轮数，每轮结束后检查：

```text
本轮发现的 target 增量
上两轮发现的 target 增量
```

满足快速收敛条件时停止。

还会检查：

- 相邻两轮改进是否已经很小
- 发现的 target 是否停止增长
- 当前 fold 的方向是否合理
- 是否发生性能退化

相关逻辑位于 `CrossValidation::train`。

---

## 10. 合并所有 fold 的测试分数

迭代结束后调用：

```cpp
crossValidation.postIterationProcessing(allScores, pCheck_);
```

主要方法：

```cpp
Scores::merge()
Scores::normalizeScores()
Normalizer::endScoreNormalizeWeights()
Scores::postMergeStep()
Scores::calcQvals()
```

`postIterationProcessing()` 位于 `CrossValidation.cpp`。

执行：

1. 对每个 fold 的测试分数合并。
2. 计算每个 fold 的 q-value。
3. 每个 fold 分数进行统一仿射重标定。
4. 合并所有 fold。
5. 重新排序。
6. 重新统计 target/decoy 数量。
7. 棥查分离度和 pi0。
8. 验证最终模型方向。

分数重标定方法是：

```cpp
Scores::normalizeScores()
```

目标：

```text
q = trainFDR 时分数为 0
decoy 中位数分数为 -1
```

相关实现位于 `Scores.cpp`。

同时会修改权重：

```cpp
newWeight = (oldWeight / scale)
newIntercept = (oldIntercept - subScore) / scale
```

由：

```cpp
Normalizer::endScoreNormalizeWeights()
```

实现。

---

## 11. 可选 RESET

如果启用 RESET，`createXvalSetsBySpectrum()` 会额外：

1. 随机选择一部分 decoy 进入测试时身份。
2. 将这些 decoy 标记为 `PSEUDO_TARGET`。
3. 调整 target/decoy 比例。
4. 更新测试和训练集合的 null-target win probability。

代码位于 `Scores::createXvalSetsBySpectrum`。

相关概念：

```text
normal decoy -> 负样本
pseudo target -> RESET 伪目标
```

它主要改变训练集分布，不属于基础 3-fold Cross Validation 的必须步骤。

---

## 12. 计算最终结果

在合并所有 fold 后，Percolator 继续执行：

```cpp
allScores.postMergeStep();
allScores.calcQvals(selectionFdr_);
allScores.normalizeScores(selectionFdr_, rawWeights);
```

核心入口位于 `Caller.cpp`。

方法包括：

- `Scores::postMergeStep()`
- `Scores::calcQvals()`
- `Scores::normalizeScores()`
- `Scores::merge()`
- `calcAndOutputResult()`

`calcAndOutputResult()` 负责：

- 生成 PSM-level q-value
- 生成 peptide-level q-value
- 生成 protein-level q-value
- TDA/TDC FDR 统计
- PEP、posterior probability
- 写 pin、psms、peptides、proteins 输出

入口位于 `Caller::calcAndOutputResult`。

## Java 当前应至少补齐的关键流程

要接近 Percolator，而不只是普通二分类：

1. 按 `ScanNr` 分成 3-fold。
2. 每个 fold 使用独立归一化参数和模型。
3. 从初始特征方向或权重文件初始化。
4. 每轮使用 q-value 选择正样本。
5. 所有 decoy作为负样本。
6. 在 nested fold 中网格搜索 `Cpos/Cfrac`。
7. 收敛后合并每个 fold 的测试分数。
8. 按 decoy 中位数和目标 FDR 阈值重标定分数。
9. 最后计算 Percolator-compatible q-value。
10. 如果启用 `-N`，训练子集后重新读取并打分完整数据集。

其中最关键的不是 target/decoy 合并，而是：

```text
未参与某一 fold 训练的 PSM 才对该 fold 输出测试分数
```

以及：

```text
迭代选择高置信 target，而不是一次性用全部数据训练
```