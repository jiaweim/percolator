#include "SetHandler.h"

// ⭐2026-11-12

/**
 * 构造函数
 * maxPSMs_ == 0：读取全部 PSM。
 * maxPSMs_ > 0：最多保留指定数量的 PSM。
 */
SetHandler::SetHandler(unsigned int maxPSMs) : maxPSMs_(maxPSMs) {}

/**
 * 删除所有 Dataset
 */
SetHandler::~SetHandler() {
  reset();
}

/**
 * ⭐2026-11-09
 * 遍历 subsets，删除所有 Dataset
 */
void SetHandler::reset() {
  for (unsigned int ix = 0; ix < subsets_.size(); ix++) {
    if (subsets_[ix] != NULL) {
      delete subsets_[ix];
    }
    subsets_[ix] = NULL;
  }
  subsets_.clear();
  // 清空数据集列表，充值全局特征名称列表。
  DataSet::resetFeatureNames();
}

/**
 * ⭐2026-11-09
 * 查找标签对应的数据集。
 * 
 * Gets the vector index of the DataSet matching the label
 * @param label DataSet label
 * @return index of matching DataSet
 */
unsigned int SetHandler::getSubsetIndexFromLabel(LabelType label) {
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    if (subsets_[ix]->getLabel() == label) return ix;
  }
  ostringstream temp;
  temp << "Error: No DataSet found with label " << label << std::endl;
  throw MyException(temp.str());
}

/**
 * ⭐2026-11-09，2026-11-12
 * 添加一个数据集。
 * Insert DataSet object into this SetHandler
 * @param ds pointer to DataSet to be inserted
 */
void SetHandler::push_back_dataset( DataSet * ds ) {
  subsets_.push_back(ds);
}

/**
 * ⭐2026-11-09
 * 找到对应数据集后，把其中的 PSM 转换为 ScoreHolder 对象，加入 scores。
 * @param scores ScoreHolder 对象列表
 */
void SetHandler::populateScoresWithPSMs(vector<ScoreHolder> &scores, LabelType label) {
  subsets_[getSubsetIndexFromLabel(label)]->fillScores(scores);
}

/**
 * ⭐2026-11-06，2026-11-12
 * 特征归一化。
 */
void SetHandler::normalizeFeatures(Normalizer*& pNorm) {
  // 分别保存普通特征和保留时间特征指针
  std::vector<double*> featuresV, rtFeaturesV;
  for (unsigned int ix = 0; ix < subsets_.size(); ++ix) {
    subsets_[ix]->fillFeatures(featuresV);
  }
  pNorm = Normalizer::getNormalizer();

  pNorm->setSet(featuresV,
      rtFeaturesV,
      FeatureNames::getNumFeatures(),
      0);
  pNorm->normalizeSet(featuresV, rtFeaturesV);
}

/**
 * ⭐2026-11-12
 */
LabelType const SetHandler::getLabel(int setPos) {
  assert(setPos >= 0 && setPos < (signed int)subsets_.size());
  return subsets_[static_cast<std::size_t>(setPos)]->getLabel();
}

/**
 * ⭐2026-11-05，2026-11-12
 * 读取 pin 文件入口
 */
int SetHandler::readTab(istream& dataStream, SanityCheck*& pCheck) {
  std::vector<double> noWeights; // 空权重列表
  Scores noScores(true); // 空 Scores
  return readAndScoreTab(dataStream, noWeights, noScores, pCheck);
}

/**
 * ⭐2026-11-05，2026-11-12
 * 从 headerLine 中解析可选字段，返回可选字段数量。
 * @param headerLine 表头行
 */
int SetHandler::getOptionalFields(const std::string& headerLine, 
    std::vector<OptionalField>& optionalFields) {
  TabReader reader(headerLine);
  reader.skip(2u); // discard id, label
  bool hasScannr = false, hasRt = false, hasDm = false;
  while (!reader.error()) {
    std::string optionalHeader = reader.readString();
    // transform to lower case for case insensitive matching
    std::transform(optionalHeader.begin(), optionalHeader.end(), 
                   optionalHeader.begin(), ::tolower);
    if (optionalHeader == "scannr") {
      optionalFields.push_back(SCANNR);
      hasScannr = true;
    } else if (optionalHeader == "expmass") {
      optionalFields.push_back(EXPMASS);
    } else if (optionalHeader == "calcmass") {
      optionalFields.push_back(CALCMASS);
    } else if ((optionalHeader == "rt" || optionalHeader == "retentiontime")) {
      optionalFields.push_back(RETTIME);
      hasRt = true;
    } else if ((optionalHeader == "filename" || optionalHeader == "spectrafile")) {
      optionalFields.push_back(FILENAME);
    } else {
      break;
    }
  }
  if (!hasScannr && VERB > 0) {
    cerr << "\nWARNING: Tab delimited input does not contain ScanNr column," <<
            "\n         scan numbers will be assigned automatically.\n" << endl;
  }
  return static_cast<int>(optionalFields.size());
}

/**
 * ⭐2026-11-05，2026-11-12
 * 判断 PIN 第二行是否为 DefaultDirection 行。
 * @param defaultDirectionLine PIN 第二行
 */
bool SetHandler::isDefaultDirectionLine(const std::string& defaultDirectionLine) {
  TabReader reader(defaultDirectionLine);
  
  std::string psmid = reader.readString();

  std::transform(psmid.begin(), psmid.end(), psmid.begin(), ::tolower);
  return (psmid == "defaultdirection");
}

/**
 * ⭐2026-11-05，2026-11-12
 * 通过第一行 PSM 计算特征数量。
 * @param line 第一行 PSM
 * @param optionalFieldCount 可选字段数量
 */
int SetHandler::getNumFeatures(const std::string& line, int optionalFieldCount) {
  TabReader reader(line);

  // 2: SpecId + Label
  reader.skip(static_cast<std::size_t>(2 + optionalFieldCount)); // remove id, label and optional fields
  
  reader.readDouble();
  int numFeatures = 0;
  while (!reader.error()){
    ++numFeatures;
    reader.readDouble();
  }

  return numFeatures;
}

/**
 * ⭐2026-11-05，2026-11-12
 * 从表头行中解析特征名称。
 * @param headerLine 表头行
 * @param numFeatures 特征数量
 * @param optionalFieldCount 可选字段数量
 * @param featureNames 特征名称容器
 */
void SetHandler::getFeatureNames(const std::string& headerLine, 
    int numFeatures, int optionalFieldCount, FeatureNames& featureNames) {
  TabReader reader(headerLine);
  // removes enumerator, label and if present optional fields
  reader.skip(static_cast<std::size_t>(2 + optionalFieldCount));
  int numFeatLeft = numFeatures;
  while (!reader.error() && numFeatLeft > 0) {
    std::string tmp = reader.readString();
    if (!tmp.empty()) { // Ensure that the string is not empty or consider other validity checks
      featureNames.insertFeature(tmp);
      numFeatLeft--;
    }
  }
  // Check if we have read the correct number of features
  if (numFeatLeft != 0) {
      cerr << "Not enough features read" << endl;
      exit(-1);
  } 
  featureNames.initFeatures();
  assert(numFeatures == DataSet::getNumFeatures());
}

/**
 * ⭐2026-11-05，2026-11-12
 * 从 DefaultDirection 获取每个特征的初始权重
 * @param defaultDirectionLine DefaultDirection 行
 * @param optionalFieldCount 可选字段数量
 * @param init_values 初始值容器
 * @return 是否存在默认值
 */
bool SetHandler::getInitValues(const std::string& defaultDirectionLine, 
    int optionalFieldCount, std::vector<double>& init_values) {
  TabReader reader(defaultDirectionLine);
  // removes enumerator, label and if present optional fields
  reader.skip(static_cast<std::size_t>(2 + optionalFieldCount));
  
  bool hasDefaultValues = false;
  unsigned int ix = 0;
  double a = reader.readDouble();
  while (!reader.error()) {
    if (a != 0.0) hasDefaultValues = true;
    if (VERB > 2) {
      std::cerr << "Initial direction for " << 
                   DataSet::getFeatureNames().getFeatureName(ix) << " is " << 
                   a << std::endl;
    }
    init_values.push_back(a);
    a = reader.readDouble();
    ix++;
  }
  return hasDefaultValues;
}

/**
 * ⭐2026-11-12
 * 输出数据集
 */
void SetHandler::writeTab(const string& dataFN, SanityCheck * pCheck) {
  ofstream dataStream(dataFN.data(), ios::out);
  dataStream << "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\t";
  if (PSMDescription::hasSpectrumFileName()) {
    dataStream << "filename\t";
  }
  dataStream << DataSet::getFeatureNames().getFeatureNames()
      << "\tPeptide\tProteins" << std::endl;
  vector<double> initial_values = pCheck->getDefaultWeights();
  if (initial_values.size() > 0) {
    dataStream << "DefaultDirection\t-\t-\t-\t-";
    for (size_t i = 0; i < initial_values.size(); ++i) {
      dataStream << '\t' << initial_values[i];
    }
    dataStream << std::endl;
  }
  for (std::vector<DataSet*>::iterator it = subsets_.begin();
         it != subsets_.end(); ++it) {
    (*it)->writeTabData(dataStream);
  }
}

/**
 * ⭐2026-11-12
 * 解析和注册 PSM
 * @param dataStream 输入流
 * @param psmLine 第一套 PSM 行
 * @param hasInitialValueRow 是否存在 DefaultDirection 行
 * @param concatenatedSearch 是否为连接搜索
 * @param optionalFields 可选字段列表
 */
void SetHandler::readPSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, bool& concatenatedSearch,
    std::vector<OptionalField>& optionalFields) {

  // 分为 target 和 decoy 的目的：
  // - 统计 target 和 decoy 数量
  // - 计算 target/decoy 比例
  // - 生成正负训练样本
  // - 按标签筛选目标
  // - 单独输出 target/decoy 分数
  // - 拆分 CV fold
  // - 根据 decoy 估计 FDR
  DataSet* targetSet = new DataSet();
  assert(targetSet);
  targetSet->setLabel(LabelType::TARGET);
  DataSet* decoySet = new DataSet();
  assert(decoySet);
  decoySet->setLabel(LabelType::DECOY);
  
  unsigned int lineNr = (hasInitialValueRow ? 3u : 2u);
  if (psmLine.size() == 0) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, could not find any PSMs." << std::endl;
    if (NO_TERMINATE) {
      cerr << temp.str() << "No-terminate flag set: ignoring error and continuing "
          << "without PSMs." << std::endl;
      
    } else {
      throw MyException(temp.str());
    }
  } else if (maxPSMs_ > 0u) { // reservoir sampling to create subset of size maxPSMs_
    std::priority_queue<PSMDescriptionPriority> subsetPSMs;
    // ScanId -> (priority, isDecoy)
    std::map<ScanId, std::pair<size_t, bool> > scanIdLookUp;
    unsigned int upperLimit = UINT_MAX;
    do {
      if (lineNr % 1000000 == 0 && VERB > 1) {
        std::cerr << "Processing line " << lineNr << std::endl;
      }
      psmLine = rtrim(psmLine);
      
      int label = 0;
      ScanId scanId = getScanId(psmLine, label, optionalFields, lineNr);
      bool isDecoy = (label == -1);
      size_t randIdx;
      if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
        if (concatenatedSearch && isDecoy != scanIdLookUp[scanId].second) {
          concatenatedSearch = false;
        }
        randIdx = scanIdLookUp[scanId].first;
      } else {
        randIdx = PseudoRandom::lcg_rand();
        scanIdLookUp[scanId].first = randIdx;
        scanIdLookUp[scanId].second = isDecoy;
      }
      
      if (subsetPSMs.size() < maxPSMs_ || randIdx < upperLimit) {
        PSMDescriptionPriority psmPriority;
        bool readProteins = false;
        psmPriority.label = DataSet::readPsm(psmLine, lineNr, optionalFields, 
                                 readProteins, psmPriority.psm, featurePool_, decoyPrefix_);
        psmPriority.priority = randIdx;
        subsetPSMs.push(psmPriority);
        if (subsetPSMs.size() > maxPSMs_) {
          PSMDescriptionPriority del = subsetPSMs.top();
          upperLimit = static_cast<unsigned int>(del.priority);
          featurePool_.deallocate(del.psm->features);
          PSMDescription::deletePtr(del.psm);
          subsetPSMs.pop();
        }
      }
      ++lineNr;
    } while (getline(dataStream, psmLine));
    
    addQueueToSets(subsetPSMs, targetSet, decoySet);
  } else { // simply read all PSMs
    std::map<ScanId, bool> scanIdLookUp; // ScanId -> isDecoy
    do {
      if (lineNr % 1000000 == 0 && VERB > 1) {
        std::cerr << "Reading line " << lineNr << std::endl;
      }
      psmLine = rtrim(psmLine);
      int label = 0;
      ScanId scanId = getScanId(psmLine, label, optionalFields, lineNr);
      bool isDecoy = (label == -1);
      if (scanIdLookUp.find(scanId) != scanIdLookUp.end()) {
        if (concatenatedSearch && isDecoy != scanIdLookUp[scanId]) {
          concatenatedSearch = false;
        }
      } else {
        scanIdLookUp[scanId] = isDecoy;
      }
      if (label == 1) {
        targetSet->readPsm(psmLine, lineNr, optionalFields, featurePool_, decoyPrefix_);
      } else if (label == -1) {
        decoySet->readPsm(psmLine, lineNr, optionalFields, featurePool_, decoyPrefix_);
      } else {
        std::cerr << "Warning: the PSM on line " << lineNr
            << " has a label not in {1,-1} and will be ignored." << std::endl;
      }
      ++lineNr;
    } while (getline(dataStream, psmLine));
  }
  
  if (VERB > 1) {
    std::cerr << "Found " << lineNr - (hasInitialValueRow ? 3u : 2u) << " PSMs" << std::endl;
  }
  
  push_back_dataset(targetSet);
  push_back_dataset(decoySet);
}

/**
 * ⭐2026-11-12
 * 将采样序列中的 PSM 按标签分发到 target 和 decoy Dataset.
 */
void SetHandler::addQueueToSets(
    std::priority_queue<PSMDescriptionPriority>& subsetPSMs,
    DataSet* targetSet, DataSet* decoySet) {
  while (!subsetPSMs.empty()) {
    PSMDescriptionPriority psmPriority = subsetPSMs.top();
    if (psmPriority.label == LabelType::TARGET) {
      targetSet->registerPsm(psmPriority.psm);
    } else if (psmPriority.label == LabelType::DECOY) {
      decoySet->registerPsm(psmPriority.psm);
    } else {
      std::cerr << "Warning: the PSM " << psmPriority.psm->getId()
          << " has a label not in {1,-1} and will be ignored." << std::endl;
      featurePool_.deallocate(psmPriority.psm->features);
      PSMDescription::deletePtr(psmPriority.psm);
    }
    subsetPSMs.pop();
  }
}

/**
 * ⭐2026-11-05, 2026-11-12
 * 从 PSM 行中解析 ScanId，返回 ScanId 对象。
 * @param psmLine PSM 行
 * @param label 输出参数，PSM 标签
 * @param optionalFields 可选字段列表
 * @param lineNr 当前行号
 * @return ScanId 对象
 */
ScanId SetHandler::getScanId(const std::string& psmLine, int& label,
    std::vector<OptionalField>& optionalFields, unsigned int lineNr) {
  ScanId scanId;
  TabReader reader(psmLine);
  
  reader.skip(); // 跳过 psmId
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr 
        << ". Could not read PSMid." << std::endl;
    throw MyException(temp.str());
  }
  
  label = reader.readInt();
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM on line " << lineNr 
        << ". Could not read label." << std::endl;
    throw MyException(temp.str());
  }
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        scanId.first = reader.readInt();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number on line " 
              << lineNr << ". Check if scan number is an integer." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        scanId.second = reader.readDouble();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading experimental mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case CALCMASS: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading calculated mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case RETTIME: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading retention time on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case DELTAMASS: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading delta mass on line " 
              << lineNr << ". Check if experimental mass is a floating point number." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } case FILENAME: {
        reader.skip();
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading spectra file name on line " 
              << lineNr << "." << std::endl;
          throw MyException(temp.str());
        }
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) scanId.first = static_cast<int>(lineNr);
  
  return scanId;
}

/**
 * ⭐2026-11-12
 * 读取 tab 文件实现
 * @param dataStream 输入数据流
 * @param rawWeights 已有模型权重，非空时，读取 PSM 时直接计算分数，空时，读取 PSM 并进行校验。暂时不需要
 * @param allScores 直接读取并打分时保存打分结果
 * @param pCheck 校验对象，读取 PSM 并进行校验时返回
 * @return 0 读取失败，1 读取成功
 */
int SetHandler::readAndScoreTab(istream& dataStream, 
    std::vector<double>& rawWeights, Scores& allScores, SanityCheck*& pCheck) {
  if (!dataStream) {
    std::cerr << "ERROR: Cannot open data stream." << std::endl;
    return 0;
  }
  std::string psmLine, headerLine, defaultDirectionLine;
  
  // 1. 读取表头
  getline(dataStream, headerLine); // line with feature names
  headerLine = rtrim(headerLine);
  if (headerLine.substr(0,5) == "<?xml") {
    std::cerr << "ERROR: Cannot read Tab delimited input from data stream.\n" << 
       "Input file seems to be in XML format, which is deprecated since percolator 3.9." << 
       std::endl;
    return 0;
  }
  
  // 2. 解析可选字段
  // Checking for optional headers "ScanNr", "ExpMass", "CalcMass", "Rt"/"RetentionTime" and "dM"/"DeltaMass"
  std::vector<OptionalField> optionalFields;
  int optionalFieldCount = getOptionalFields(headerLine, optionalFields);

  // 3. 判断第二行是否为 DefaultDirection 行
  // parse second line for default direction
  getline(dataStream, defaultDirectionLine);
  defaultDirectionLine = rtrim(defaultDirectionLine);
  bool hasInitialValueRow = isDefaultDirectionLine(defaultDirectionLine);

  // 4. 从第一条 PSM
  if (hasInitialValueRow) {
    getline(dataStream, psmLine);
  } else {
    psmLine = defaultDirectionLine;
  }
  psmLine = rtrim(psmLine);
  // 5. 计算特征数量
  int numFeatures = getNumFeatures(psmLine, optionalFieldCount);
  FeatureNames& featureNames = DataSet::getFeatureNames();
  // 6. 解析特征名称 fill in the feature names from the header line
  getFeatureNames(headerLine, numFeatures, optionalFieldCount, featureNames);

  // 7. 检查特征是否足够
  if (numFeatures < 1) {
    ostringstream oss;
    oss << "ERROR: Reading tab file, too few features present." << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error and using a pseudo-feature of zeroes." << std::endl;
      featurePool_.createPool(1u);
    } else {
      throw MyException(oss.str());
    }
  } else {
    featurePool_.createPool(DataSet::getNumFeatures());
  }  
  
  // 8. 读取初始权重 fill in the default weights if present
  std::vector<double> init_values;
  bool hasDefaultValues = false;
  if (hasInitialValueRow) {
    hasDefaultValues = getInitValues(defaultDirectionLine, optionalFieldCount, 
                                     init_values);
  }
  if (hasDefaultValues && init_values.size() > static_cast<std::size_t>(numFeatures)) {
    ostringstream oss;
    oss << "ERROR: Reading tab file, too many default values present." << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error and trimming default value vector." << std::endl;
      init_values.resize(static_cast<std::size_t>(numFeatures));
    } else {
      throw MyException(oss.str());
    }
  }

  // 9.1 已有权重，直接读取并打分，跳过
  // read in the data
  if (rawWeights.size() > 0) {
    readAndScorePSMs(dataStream, psmLine, hasInitialValueRow, optionalFields, rawWeights, allScores);
  } else {
    // 9.2 没有权重，读取 PSM 并进行校验，普通的训练集读取
    // detect if the input came from separate target and decoy searches or 
    // from a concatenated search by looking for scan+expMass combinations
    // that have both at least one target and decoy PSM
    bool concatenatedSearch = true;
    
    readPSMs(dataStream, psmLine, hasInitialValueRow, concatenatedSearch, optionalFields);
    
    pCheck = new SanityCheck();
    pCheck->checkAndSetDefaultDir();
    if (hasDefaultValues) pCheck->addDefaultWeights(init_values); 
    pCheck->setConcatenatedSearch(concatenatedSearch);
  }
  return 1;
}

/**
 * ⭐2026-11-05, 2026-11-12
 * 已有权重，读取 PSM 并打分
 */
void SetHandler::readAndScorePSMs(istream& dataStream, std::string& psmLine, 
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields, 
    std::vector<double>& rawWeights, Scores& allScores) {
  unsigned int lineNr = (hasInitialValueRow ? 3u : 2u);
  bool readProteins = true;
  do {
    if (lineNr % 1000000 == 0 && VERB > 1) {
      std::cerr << "Processing line " << lineNr << std::endl;
    }
    psmLine = rtrim(psmLine);
    ScoreHolder sh;
    sh.label = DataSet::readPsm(psmLine, lineNr, optionalFields, readProteins, sh.pPSM, featurePool_, decoyPrefix_);
    allScores.scoreAndAddPSM(sh, rawWeights, featurePool_);
    ++lineNr;
  } while (getline(dataStream, psmLine));
  
  if (VERB > 1) {
    std::cerr << "Found " << lineNr - (hasInitialValueRow ? 3u : 2u) << " PSMs" << std::endl;
  }
}

/**
 * ⭐2026-11-12
 * 删除空格
 */
std::string& SetHandler::rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
  return s;
}
