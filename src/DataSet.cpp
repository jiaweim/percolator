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

#include "DataSet.h"
#include <cmath>

FeatureNames DataSet::featureNames_;
bool DataSet::decoyWarningTripped_ = false;

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

/*const double* DataSet::getFeatures(const int pos) const {
  return &feature[pos];
}*/

bool DataSet::writeTabData(ofstream& out) {
  unsigned int nf = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    out << psm->getId() << '\t' << label_ << '\t' << psm->scan << '\t' 
        << psm->expMass << '\t' << psm->calcMass;
    if (psm->hasSpectrumFileName()) {
      out << '\t' << psm->getSpectrumFileName();
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    out << '\t' << psm->getFullPeptide() << '\t';
    psm->printProteins(out);
    out << endl;
  }
  return true;
}

void DataSet::print_features() {
  for (std::size_t i = 0; i < getSize(); i++) {
    for (std::size_t j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << j + 1 << ":" << psms_[i]->features[j] << " ";
    }
    cerr << endl;
  }
}

void DataSet::print_10features() {
  cerr << DataSet::getFeatureNames().getFeatureNames() << endl;
  for (std::size_t i = 0; i < 10; i++) {
    for (std::size_t j = 0; j < FeatureNames::getNumFeatures(); j++) {
      cerr << psms_[i]->features[j] << "\t";
    }
    cerr << endl;
  }
}

// TODO: find a way to make these three functions generic
void DataSet::fillScores(std::vector<ScoreHolder>& scores) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    scores.push_back(ScoreHolder(.0, label_, psm));
  }
}

void DataSet::fillFeatures(std::vector<double*>& features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features);
  }
}

/**
 * 从 PIN 读取一行，得到 PSMDescription 对象，并注册到 DataSet 中
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
  TabReader reader(line);
  std::string tmp;
  
  myPsm = new PSMDescription();
  myPsm->setId(reader.readString()); // 1. id
  LabelType label = reader.readInt() == 1 ? LabelType::TARGET : LabelType::DECOY; // 2. label
  
  // 读取可选字段，可选字段顺序必须与文件列顺序一致
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
      } case DELTAMASS: {
          ostringstream temp;
          temp << "ERROR: The parser tries to allocate a deprecated delta mass field for the PSM " 
              << myPsm->getId() << std::endl;
          throw MyException(temp.str());
        break;
      } case FILENAME: {
        myPsm->setSpectrumFileName(reader.readString());
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) myPsm->scan = lineNr;
  
  // 读取特征向量
  unsigned int numFeatures = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  double* featureRow = featurePool.allocate(); // 从内存池分配特征数组
  myPsm->features = featureRow;
  for (unsigned int j = 0; j < numFeatures; j++) {
    featureRow[j] = reader.readDouble();
    if (!isfinite(featureRow[j])) { // NaN, Inf, -Inf 会破坏 SVM 训练，必须拒绝
      ostringstream oss;
      oss << "ERROR: Reached strange feature with val=" << featureRow[j]
           << " col=" << j << " for PSM with id " << myPsm->getId() << endl;
      if (NO_TERMINATE) { // NO_TERMINATE 模式下，打印警告并将该特征值设为 0
        std::cerr << oss.str();
        std::cerr << "No-terminate flag set: setting value to 0 and ignoring the error." << std::endl;
        featureRow[j] = 0.0;
      } else {
        throw MyException(oss.str() + "Terminating.\n");
      }
    }
  }
  // 错误检查放在循环外，如果 readDouble 因为字段不足或格式错误而失败，err 被设置，但上面的循环仍然会继续，直到读取完所有特征列
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading in feature vector of PSM " 
      << myPsm->getId() << ". Check if there are enough features on this line and "
      << "if they are all floating point numbers or integers." << std::endl;
    throw MyException(temp.str());
  }
  
  // 解析肽段
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
    } else if (peptide_seq.at(1) != '.' && peptide_seq.at(peptide_seq.size()-1) != '.') { // 检查肽段是否包含前后缀
      ostringstream temp;
      temp << "ERROR: Reading tab file, the peptide sequence " << peptide_seq 
        << " with PSM id " << myPsm->getId() << " does not contain one or two of its"
        << " flanking amino acids." << std::endl;
      throw MyException(temp.str());
    }
  }
  
  // 蛋白质解析
  if (readProteins) {
    std::vector<std::string> proteins;
    if (PSMDescription::getProteinNameSeparator() == "\t") {
      // 每个蛋白质一个 tab 字段
      while (!reader.error()) {
        std::string tmp = reader.readString();
        if (tmp.size() > 0) proteins.push_back(tmp);
      }
    } else {
      // 所有蛋白质在一个字段中，用指定的分隔符分隔
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
    proteins.swap(myPsm->proteinIds); // shrink to fit
  }

  // 检查 decoy 前缀是否匹配
  if (label == LabelType::DECOY) {
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
  switch (label_) {
    case LabelType::TARGET: { break; };
    case LabelType::DECOY: { break; };
    default:  { throw MyException("ERROR : Reading PSM, class DataSet has not been initiated\
    to neither target nor decoy label\n");}
  }
  psms_.push_back(myPsm);
}
