#include "ScoreHolder.h"
#include "Scores.h"

#include <algorithm>
#include <boost/algorithm/string/trim.hpp>
#include <iomanip>
#include <regex>

/**
 * 通用排序：⭐
 * 1. score，降序
 * 2. scan，降序
 * 3. expMass，降序
 * 4. label，降序
 */
bool operator>(const ScoreHolder& one, const ScoreHolder& other) {
  return (one.score > other.score) ||
         (one.score == other.score && one.pPSM->scan > other.pPSM->scan) ||
         (one.score == other.score && one.pPSM->scan == other.pPSM->scan &&
          one.pPSM->expMass > other.pPSM->expMass) ||
         (one.score == other.score && one.pPSM->scan == other.pPSM->scan &&
          one.pPSM->expMass == other.pPSM->expMass && one.label > other.label);
}

/**
 * 通用排序：⭐
 * 1. score，升序
 * 2. scan，升序
 * 3. expMass，升序
 * 4. label，升序
 */
bool operator<(const ScoreHolder& one, const ScoreHolder& other) {
  return (one.score < other.score) ||
         (one.score == other.score && one.pPSM->scan < other.pPSM->scan) ||
         (one.score == other.score && one.pPSM->scan == other.pPSM->scan &&
          one.pPSM->expMass < other.pPSM->expMass) ||
         (one.score == other.score && one.pPSM->scan == other.pPSM->scan &&
          one.pPSM->expMass == other.pPSM->expMass && one.label < other.label);
}

/**
 * 过滤掉 ASCII 码 <32 的不可打印字符（控制字符），防止它们破坏 XML 输出格式。
 */
string getRidOfUnprintablesAndUnicode(string inpString) {
  string outputs = "";
  for (unsigned int jj = 0; jj < inpString.size(); jj++) {
    signed char ch = inpString[jj];
    // NOTE signed char ranges -128 to 127
    if (((int)ch) >= 32) {
      outputs += ch;
    }
  }
  return outputs;
}

/**
 * 
 */
void ScoreHolder::printPSM(ostream& os, bool printDecoys, bool printExpMass) {
  if (!isDecoy() || printDecoys) {
    os << "    <psm p:psm_id=\"" << pPSM->getId() << "\"";
    if (printDecoys) {
      if (isDecoy())
        os << " p:decoy=\"true\"";
      else
        os << " p:decoy=\"false\"";
    }
    os << ">" << endl;
    os << "      <svm_score>" << fixed << score << "</svm_score>" << endl; // svm_score，精确到小数点后 6 位
    os << "      <q_value>" << scientific << q << "</q_value>" << endl; // q_value，科学计数法表示
    os << "      <pep>" << scientific << pep << "</pep>" << endl; // pep，科学计数法表示

    if (printExpMass) {
      os << "      <exp_mass>" << fixed << std::setprecision(4) << pPSM->expMass // exp_mass，精确到小数点后 4 位
         << "</exp_mass>" << endl;
    }

    os << "      <calc_mass>" << fixed << std::setprecision(3) << pPSM->calcMass // calc_mass，精确到小数点后 3 位
       << "</calc_mass>" << endl;
    /* Remove this tag for now */
    /* if (isfinite(pPSM->getRetentionTime())) {
          os << "      <retentionTime>" << fixed << std::setprecision (3) <<
       pPSM->getRetentionTime() << "</retentionTime>" << endl;
        } */

    if (pPSM->getPeptideSequence().size() > 0) {
      string n = pPSM->getFlankN();
      string c = pPSM->getFlankC();
      string centpep = pPSM->getPeptideSequence();
      os << "      <peptide_seq n=\"" << n << "\" c=\"" << c << "\" seq=\""
         << centpep << "\"/>" << endl;
    }

    std::vector<std::string>::const_iterator pidIt = pPSM->proteinIds.begin();
    for (; pidIt != pPSM->proteinIds.end(); ++pidIt) {
      os << "      <protein_id>" << getRidOfUnprintablesAndUnicode(*pidIt)
         << "</protein_id>" << endl;
    }

    os << "      <p_value>" << scientific << p << "</p_value>" << endl;
    os << "    </psm>" << endl;
  }
}

/**
 * 与 printPSM 类似，但输出的是 peptide 信息.
 * 
 */
void ScoreHolder::printPeptide(ostream& os,
                               bool printDecoys,
                               bool printExpMass,
                               Scores& fullset) {
  if (!isDecoy() || printDecoys) {
    os << "    <peptide p:peptide_id=\"" << pPSM->getPeptideSequence() << "\"";
    if (printDecoys) {
      if (isDecoy())
        os << " p:decoy=\"true\"";
      else
        os << " p:decoy=\"false\"";
    }
    os << ">" << endl;

    os << "      <svm_score>" << fixed << score << "</svm_score>" << endl;
    os << "      <q_value>" << scientific << q << "</q_value>" << endl;
    os << "      <pep>" << scientific << pep << "</pep>" << endl;

    if (printExpMass) {
      os << "      <exp_mass>" << fixed << std::setprecision(4) << pPSM->expMass
         << "</exp_mass>" << endl;
    }
    os << "      <calc_mass>" << fixed << std::setprecision(3) << pPSM->calcMass
       << "</calc_mass>" << endl;

    std::vector<std::string>::const_iterator pidIt = pPSM->proteinIds.begin();
    for (; pidIt != pPSM->proteinIds.end(); ++pidIt) {
      os << "      <protein_id>" << getRidOfUnprintablesAndUnicode(*pidIt)
         << "</protein_id>" << endl;
    }

    os << "      <p_value>" << scientific << p << "</p_value>" << endl;
    os << "      <psm_ids>" << endl;

    // output all psms that contain the peptide
    std::vector<PSMDescription*>::const_iterator psmIt =
        fullset.getPsms(pPSM).begin();
    for (; psmIt != fullset.getPsms(pPSM).end(); ++psmIt) {
      os << "        <psm_id>" << (*psmIt)->getId() << "</psm_id>" << endl;
    }
    os << "      </psm_ids>" << endl;
    os << "    </peptide>" << endl;
  }
}

/**
 * 从 PSM ID 中提取电荷信息
 */
std::string ScoreHolder::getCharge(std::string id) {
  /* Aposymb_Proteome_DIA_RAW_A01_Q1.00148.00148.1_2 -> 1_2*/
  std::string last_element(id.substr(id.rfind("_") - 1));
  /* 1_2 -> 1 */
  return last_element.substr(0, last_element.find("_"));
}

/**
 * 兼容 TPP 的 pepXML 输出格式，输出 spectrum_query 标签
 */
void ScoreHolder::printPepXML(ostream& os,
                              map<char, float>& aaWeight,
                              int index) {
  /* std::cerr << pepXMLBaseName << std::endl; */
  std::string id = pPSM->getId();
  /* Get scan ids */
  unsigned int scan = pPSM->scan;
  unsigned int native_id = scan - 1;
  /* Get charge */
  std::string assumed_charge = getCharge(id);
  /* Get RT */
  double RT = pPSM->getRetentionTime();
  /*  uncalibrated_precursor_neutral_mass ? */
  double expMass = pPSM->expMass;
  /*  precursor_neutral_mass ? */
  double calcMass = pPSM->calcMass;

  os << "        <spectrum_query spectrum=\"" << id
     << "\" precursor_neutral_mass=\"" << expMass << "\" assumed_charge=\""
     << assumed_charge << "\" end_scan=\"" << scan << "\" index=\"" << index
     << "\" retention_time_sec=\"" << fixed << std::setprecision(3)
     << pPSM->getRetentionTime() << "\" start_scan=\"" << scan << "\">" << endl;
  std::string centpep = pPSM->getPeptideSequence();
  std::string trimmed_pep = boost::algorithm::trim_left_copy_if(
      centpep, boost::algorithm::is_any_of("n"));
  regex r("\\[(.*?)\\]");
  std::string peptide_sequence = regex_replace(trimmed_pep, r, "");
  os << "            <search_result>" << endl;

  /* Print protein information */
  size_t n_protein = 0;
  /* Placeholders */
  int hit_rank = 1;
  int massdiff = 1;
  /* num_tot_proteins */
  int num_tot_proteins = pPSM->proteinIds.size();

  std::vector<std::string>::const_iterator pidIt = pPSM->proteinIds.begin();
  for (; pidIt != pPSM->proteinIds.end(); ++pidIt) {
    if (n_protein == 0) {
      /*  set calc_neutral_pep_mass  as calcMass as placeholder for now */
      os << "                <search_hit calc_neutral_pep_mass=\"" << calcMass
         << "\" num_tot_proteins=\"" << num_tot_proteins << "\" hit_rank=\""
         << hit_rank << "\" massdiff=\"" << massdiff << "\" peptide=\""
         << peptide_sequence << "\" protein=\""
         << getRidOfUnprintablesAndUnicode(*pidIt) << "\">" << endl;
    } else {
      os << "                    <alternative_protein protein=\""
         << getRidOfUnprintablesAndUnicode(*pidIt) << "\"/>" << endl;
    }
    n_protein++;
  }
  /* Print modification */
  string subject(trimmed_pep);
  smatch match;
  long int mod_pos = 0;
  float mod_weight;
  size_t n_mod = 0;

  while (regex_search(subject, match, r)) {
    mod_pos += match.position(0);
    // suffix to find the rest of the string.
    if (mod_pos == 0) {
      mod_weight = std::stof(match.str(1)) + 1.0074;
      if (n_mod == 0) {
        os << "                    <modification_info mod_nterm_mass=\""
           << round(mod_weight * 1000) / 1000 << "\">" << endl;
      }
    } else {
      mod_weight =
          aaWeight[peptide_sequence.at(mod_pos - 1)] + std::stof(match.str(1));
      if (n_mod == 0) {
        os << "                    <modification_info>" << endl;
      }
      os << "                        <mod_aminoacid_mass mass=\""
         << round(mod_weight * 1000) / 1000 << "\" position=\""
         << std::to_string(mod_pos) << "\" />" << endl;
    }
    subject = match.suffix().str();
    n_mod++;
  }
  if (n_mod != 0) {
    os << "                    </modification_info>" << endl;
  }

  /* Print Percolator information */
  os << "                    <analysis_result analysis=\"peptideprophet\">"
     << endl;
  os << "                        <peptideprophet_result probability=\""
     << scientific << 1.0 - pep << "\" />" << endl;
  os << "                    </analysis_result>" << endl;
  os << "                </search_hit>" << endl;
  os << "            </search_result>" << endl;
  os << "        </spectrum_query>" << endl;
}