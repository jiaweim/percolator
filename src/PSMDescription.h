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
#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Enzyme.h"

/*
 * PSMDescription⭐2026-11-02
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 */
class PSMDescription {
   public:
    PSMDescription();
    PSMDescription(const std::string& peptide);

    virtual ~PSMDescription();
    static void deletePtr(PSMDescription* psm);
    virtual void deleteRetentionFeatures() {}

    void clear() { proteinIds.clear(); }
    double* getFeatures() { return features; }

    // TODO: move these static functions somewhere else
    static std::string removePTMs(const std::string& peptideSeq);
    static std::string removeFlanks(const std::string& peptideSeq) {
        return peptideSeq.substr(2, peptideSeq.size() - 4);
    }
    static bool ptrLess(PSMDescription* one, PSMDescription* other) {
        return *one < *other;
    }
    static bool ptrEqual(PSMDescription* one, PSMDescription* other) {
        return *one == *other;
    }

    /**
     * @return 肽段序列，不包含翼氨酸酸
     */
    std::string getPeptideSequence() { return peptide.substr(2, peptide.size() - 4); }
    /**
     * @return 肽段序列，包含翼氨酸酸
     */
    std::string& getFullPeptideSequence() { return peptide; }
    /**
     * @return 肽段的 N 端翼氨酸
     */
    std::string getFlankN() { return peptide.substr(0, 1); }
    /**
     * @return 肽段的 C 端翼氨酸
     */
    std::string getFlankC() { return peptide.substr(peptide.size() - 1, peptide.size()); }

    friend std::ostream& operator<<(std::ostream& out, PSMDescription& psm);
    void printProteins(std::ostream& out);

    /**
     * 比较运算符，用于排序 PSMDescription 对象，首先按肽段序列排序，如果肽段序列相同，则按保留时间排序
     * @param other 另一个 PSMDescription 对象
     */
    bool operator<(const PSMDescription& other) const {
        return (peptide < other.peptide) ||
               (peptide == other.peptide && getRetentionTime() < other.getRetentionTime());
    }

    bool operator==(const PSMDescription& other) const {
        return (peptide == other.peptide);
    }

    inline void setId(const std::string& id) { id_ = id; }
    inline std::string& getId() { return id_; }

    std::string& getFullPeptide() { return peptide; }
    void setPeptide(const std::string& pep_seq) { peptide = pep_seq; }
    PSMDescription* getAParent() { return this; }
    void checkFragmentPeptides(
        std::vector<PSMDescription*>::reverse_iterator other,
        std::vector<PSMDescription*>::reverse_iterator theEnd) {
        (void)other;
        (void)theEnd;
    }
    static inline void setProteinNameSeparator(const std::string sep) {
        proteinNameSeparator_ = sep;
    }
    static inline const std::string& getProteinNameSeparator() { return proteinNameSeparator_; }

    /**
     * 在池中线性查找谱图文件名，如果存在则返回索引，否则添加到池中并返回新索引
     * @param fileName 谱图文件名
     */
    void setSpectrumFileName(std::string fileName) {
        size_t index(0);
        auto specFilePos = std::find(spectraFileNames_.begin(), spectraFileNames_.end(), fileName);
        if (specFilePos != spectraFileNames_.end()) {
            index = specFilePos - spectraFileNames_.begin();
        } else {
            index = spectraFileNames_.end() - spectraFileNames_.begin();
            spectraFileNames_.push_back(fileName);
        }
        specFileNr = index;
    }
    /**
     * 
     */
    inline const std::string getSpectrumFileName() {
        std::string fn("");
        if (hasSpectrumFileName())
            fn = spectraFileNames_.at(specFileNr);
        return fn;
    }
    /**
     * 检查文件池是否为空
     */
    inline bool static hasSpectrumFileName() { return !spectraFileNames_.empty(); }

    void setRetentionFeatures(double* retentionFeatures) {(void) retentionFeatures; }
    double* getRetentionFeatures() { return NULL; }

    void setParentFragment(PSMDescription*) {}
    PSMDescription* getParentFragment() { return NULL; }

    /**
     * 设置保留时间
     * @param retentionTime 保留时间
     */
    inline void setRetentionTime(const double retentionTime) {
        retentionTime_ = retentionTime;
    }
    /**
     * @return 保留时间
     */
    inline double getRetentionTime() const { return retentionTime_; }

    /**
     * 特征向量数组的指针，默认 NULL
     */
    double* features;  // owned by a FeatureMemoryPool instance, no need to delete
    /**
     * 实验质量（0.）、计算质量（0.）、保留时间
     */
    double expMass, calcMass, retentionTime_;
    /**
     * scan number（0u），若无 scan，则用行号代替
     */
    unsigned int scan;
    /**
     * 谱图文件在 spectraFileNames_ 中的索引（0u）
     */
    unsigned int specFileNr;
    /**
     * 该 PSM 匹配的蛋白质 ID 列表
     */
    std::vector<std::string> proteinIds;

   protected:
   /**
    * PSM identifier
    */
    std::string id_;
    /**
     * 肽段序列，包含翼氨酸酸
     */
    std::string peptide;
    /**
     * 蛋白质名称分隔符，用于将多个蛋白质 ID 连接成一个字符串，默认为 \t
     */
    static std::string proteinNameSeparator_;
    /**
     * 所有谱图文件名的列表，所有 PSMDescription 共享
     */
    static vector<std::string> spectraFileNames_;
};

inline std::ostream& operator<<(std::ostream& out, PSMDescription& psm) {
    out << "Peptide: " << psm.peptide << endl;
    out << "Spectrum scan number: " << psm.scan << endl;
    out << endl;
    return out;
}

#endif /*PSMDESCRIPTION_H_*/
