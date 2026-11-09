#include "DataSet.h"
#include <cmath>

/**
 * static，PSM 特征名称表，所有 DataSet 共享同一个特征名称表
 */
FeatureNames DataSet::featureNames_;
/**
 * decoy 前缀不匹配的警告是否已经触发过，只输出一次 decoy 前缀警告。
 */
bool DataSet::decoyWarningTripped_ = false;

/**
 * 对象创建时不预分配 PSM，只初始化空 vector 和标签默认状态。
 * 构造前必须调用 setLabel()：
 */
DataSet::DataSet() {}

/**
 * 逐个清理 PSM，特征数组由 FeatureMemoryPool 管理，不需要手动释放
 */
DataSet::~DataSet() {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription::deletePtr(*it);
  }
}

/**
 * 将当前数据写入 Tab 文件
 */
bool DataSet::writeTabData(ofstream& out) {
  unsigned int nf = static_cast<unsigned int>(FeatureNames::getNumFeatures()); // 特征数量
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) { // 遍历 PSM
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    // 输出 PSM-id, label，scan, expMass, calcMass, filename.
    out << psm->getId() << '\t' << label_ << '\t' << psm->scan << '\t' 
        << psm->expMass << '\t' << psm->calcMass;
    // 输出文件名
    if (psm->hasSpectrumFileName()) {
      out << '\t' << psm->getSpectrumFileName();
    }
    // 输出特征列
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    // 输出肽段和蛋白质
    out << '\t' << psm->getFullPeptide() << '\t';
    psm->printProteins(out);
    out << endl;
  }
  return true;
}

/**
 * 将全部特征输出到 cerr，用于调试。格式为 "特征索引:特征值"，每个 PSM 一行。
 */
void DataSet::print_features() {
  for (std::size_t i = 0; i < getSize(); i++) { // 遍历 PSM
    for (std::size_t j = 0; j < FeatureNames::getNumFeatures(); j++) { // 遍历特征
      cerr << j + 1 << ":" << psms_[i]->features[j] << " "; // 1-based 特征
    }
    cerr << endl;
  }
}

/**
 * 打印前 10 个 PSM 的特征，用于调试。格式为 "特征索引:特征值"，每个 PSM 一行。
 */
void DataSet::print_10features() {
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
  for (std::size_t i = 0; i < 10; i++) {
    for (std::size_t j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << psms_[i]->features[j] << "\t";
    }
    cerr << endl;
  }
}

/**
 * ⭐2026-11-09
 * 将数据集的 PSM 转换为训练或排序所需的 ScoreHolder 对象列表。
 * Dataset 负责提供 PSM，ScoreHolder 负责保存打分和统计结果。
 */
void DataSet::fillScores(std::vector<ScoreHolder>& scores) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    scores.push_back(ScoreHolder(.0, label_, psm)); // .0 为初始打分
  }
}

/**
 * ⭐2026-11-09
 * 导出特征数组指针列表，用于训练或排序。每个 PSM 的特征数组由 FeatureMemoryPool 管理。
 */
void DataSet::fillFeatures(std::vector<double*>& features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features);
  }
}

/**
 * 从 PIN 读取一行，得到 PSMDescription 对象，并注册到 DataSet
 * 
 * Read in psm details from a string out of tab delimited file
 * @param dataStream filestream of tab delimited file, only passed to close on exception
 * TODO: remove dataStream parameter and return int with type of error to handle upstream.
 * @param line tab delimited string containing the psm details
 */
void DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, FeatureMemoryPool& featurePool,
    std::string decoyPrefix) { 
  PSMDescription* myPsm = NULL;
  bool readProteins = true;
  readPsm(line, lineNr, optionalFields, readProteins, myPsm, featurePool, decoyPrefix);
  registerPsm(myPsm);
}

/**
 * 最重要的方法，用于读取 PIN 文件中的一行，解析出 PSMDescription 对象，并返回标签类型
 * 
 * @param line 待解析 tab 行
 * @param lineNr 行号，用于错误信息定位
 * @param optionalFields 可选字段列表
 * @param readProteins 是否读取蛋白质信息
 * @param myPsm 返回的 PSMDescription 对象指针
 * @param featurePool 特征内存池，用于分配特征数组
 * @param decoyPrefix decoy 前缀，用于判断蛋白质是否为 decoy
 * @return 标签类型（TARGET/DECOY）
 */
LabelType DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool, std::string decoyPrefix) {
  TabReader reader(line); // 负责按 TAB 分割字段
  std::string tmp;
  
  myPsm = new PSMDescription();
  myPsm->setId(reader.readString()); // 1. 读取 PSM-id
  // 1 表示 TARGET，其他表示 DECOY
  LabelType label = reader.readInt() == 1 ? LabelType::TARGET : LabelType::DECOY; // 2. 读取 label
  
  // 3. 读取可选字段，可选字段顺序必须与文件列顺序一致
  bool hasScannr = false; // 只检查 scan number 是否存在，其他可选字段不存在不报错
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        myPsm->scan = static_cast<unsigned int>(reader.readInt());
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number of PSM " 
              << myPsm->getId() << ". Check if scan number is an integer." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        myPsm->expMass = reader.readDouble();
        break;
      } case CALCMASS: {
        myPsm->calcMass = reader.readDouble();
        break;
      } case RETTIME: {
        myPsm->setRetentionTime(reader.readDouble());
        break;
      } case DELTAMASS: { // 废弃字段，遇到时抛出异常
          ostringstream temp;
          temp << "ERROR: The parser tries to allocate a deprecated delta mass field for the PSM " 
              << myPsm->getId() << std::endl;
          throw MyException(temp.str());
        break;
      } case FILENAME: { // 谱图文件名
        myPsm->setSpectrumFileName(reader.readString());
        break;
      } default: { // 未知字段，抛出异常
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  // 没有 scan number，使用 PIN 行号代替
  if (!hasScannr) myPsm->scan = lineNr;
  
  // 3. 读取特征向量
  unsigned int numFeatures = static_cast<unsigned int>(FeatureNames::getNumFeatures()); // 获取全局特征数量
  double* featureRow = featurePool.allocate(); // 从内存池申请一个特征数组
  myPsm->features = featureRow;
  for (unsigned int j = 0; j < numFeatures; j++) {
    featureRow[j] = reader.readDouble(); // 读取特征
    if (!isfinite(featureRow[j])) { // NaN, Inf, -Inf 会破坏 SVM 训练，必须拒绝
      ostringstream oss;
      oss << "ERROR: Reached strange feature with val=" << featureRow[j]
           << " col=" << j << " for PSM with id " << myPsm->getId() << endl;
      if (NO_TERMINATE) { // NO_TERMINATE 模式下，打印警告并将该特征值设为 0
        std::cerr << oss.str();
        std::cerr << "No-terminate flag set: setting value to 0 and ignoring the error." << std::endl;
        featureRow[j] = 0.0;
      } else { // 否则抛出异常，终止程序
        throw MyException(oss.str() + "Terminating.\n");
      }
    }
  }

  // 读完特征后检查，错误检查放在循环外，如果 readDouble 因为字段不足或格式错误而失败，err 被设置，但上面的循环仍然会继续，直到读取完所有特征列
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading in feature vector of PSM " 
      << myPsm->getId() << ". Check if there are enough features on this line and "
      << "if they are all floating point numbers or integers." << std::endl;
    throw MyException(temp.str());
  }
  
  // 4. 解析肽段
  std::string peptide_seq = reader.readString();
  myPsm->setPeptide(peptide_seq);
  if (reader.error()) { // 肽段字段不能缺失
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM " << myPsm->getId() 
      << ". Check if a peptide and at least one protein are specified." << std::endl;
    throw MyException(temp.str());
  } else if (ProteinProbEstimator::getCalcProteinLevelProb()) {
    // MT: we only need the peptide sequences to be well formatted if protein inference is applied
    // 如果启用蛋白质推断，进行额外检查
    if (peptide_seq.size() < 5) { // 长度至少为 5，才能包含两个氨基酸的肽段和两个氨基酸的前后缀：X.序列.Y
      ostringstream temp;
      temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq 
        << " with PSM id " << myPsm->getId() << " is too short." << std::endl;
      throw MyException(temp.str());
    } else if (peptide_seq.at(1) != '.' && peptide_seq.at(peptide_seq.size()-1) != '.') { 
      // 检查肽段是否包含前后缀
      // 典型格式：K.PEPTIDER.P
      ostringstream temp;
      temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq 
        << " with PSM id " << myPsm->getId() << " does not contain one or two of its"
        << " flanking amino acids." << std::endl;
      throw MyException(temp.str());
    }
  }
  
  // 5. 蛋白质解析
  if (readProteins) {
    std::vector<std::string> proteins;

    if (PSMDescription::getProteinNameSeparator() == "\t") {
      // 每个蛋白质一个 tab 字段
      while (!reader.error()) {
        std::string tmp = reader.readString();
        if (tmp.size() > 0) proteins.push_back(tmp);
      }
    } else {
      // 自定义分隔符模式，比如 ';' 或 '|'，将蛋白质字段拆分为多个蛋白质
      std::string names = reader.readString();
      size_t pos = 0;
      std::string token;
      while ((pos = names.find(PSMDescription::getProteinNameSeparator())) != std::string::npos) {
        token = names.substr(0, pos);
        if (token.size() > 0) proteins.push_back(token);
        names.erase(0, pos + PSMDescription::getProteinNameSeparator().length());
      }
      if (names.size() > 0) proteins.push_back(names);

      // check if end of line was reached after reading the protein string
      // if this is not the case, it means the user probably incorrectly provided proteins as a tab separated input
      if (!reader.error()) {
        ostringstream temp;
        temp << "ERROR: Reading tab file, error reading PSM " << myPsm->getId() 
          << ". Found tab character after reading proteins even though " << PSMDescription::getProteinNameSeparator() 
          << " was specified as --protein-name-separator." << std::endl;
        if (NO_TERMINATE) {
          std::cerr << temp.str();
          std::cerr << "No-terminate flag set: ignoring error and continuing "
                    << "without PSMs." << std::endl;
        } else {
          throw MyException(temp.str());
        }
      }
    }
    proteins.swap(myPsm->proteinIds); // 保存，交换 vector，可以避免额外复制
  }

  // 检查 decoy 前缀是否匹配
  if (label == LabelType::DECOY) {
    // 检查每个蛋白质是否以 decoy 前缀开头，如果没有匹配，输出警告，但只输出一次
    for (auto const& proteinId: myPsm->proteinIds) { 
      bool startsWithDecoyPrefix = (proteinId.rfind(decoyPrefix, 0) == 0);
      if (!startsWithDecoyPrefix && VERB > 1 && !decoyWarningTripped_) {
        std::cerr << "Warning: protein decoy prefix " << decoyPrefix 
                  << " doesn't match the decoy protein identifier " 
                  << proteinId << "." << std::endl;
        decoyWarningTripped_ = true;
      }
    }
  }
  return label;
}

/**
 * 将 PSM 注册到数据集中
 * @param myPsm PSMDescription 对象指针
 */
void DataSet::registerPsm(PSMDescription* myPsm) {
  // 检查是否已指定标签
  switch (label_) {
    case LabelType::TARGET: { break; };
    case LabelType::DECOY: { break; };
    default:  { throw MyException("ERROR : Reading PSM, class DataSet has not been initiated\
    to neither target nor decoy label\n");}
  }
  // 标签合法，加入数据集
  psms_.push_back(myPsm);
}
