#ifndef TABFILEVALIDATOR_H_
#define TABFILEVALIDATOR_H_

#include <string>
#include <vector>
#include <random>
#include <algorithm>

#include "Globals.h"
#include "TabReader.h"

/**
 * TabFileValidator 类用于验证 PIN 文件的格式和内容⭐
 */
class TabFileValidator {
  public:
    static bool isTabFile(std::string fileName);
    static bool isTabFiles(std::vector<std::string> fileNames);
    static std::string getDecoyPrefix(std::vector<std::string> fileNames);
    static std::string detectDecoyPrefix(std::string fileName);
    static bool validateTabFiles(std::vector<std::string> fileNames, std::string &decoyPrefix);
    static void getProteinAndLabelColumnIndices(std::string fileName, int &proteinIndex, int &labelIndex);
    static std::string findDecoyPrefix(std::string fileName, int proteinIndex, int labelIndex);
    static std::string getLongestCommonPrefix(std::vector<std::string> strings);
};

#endif /*TABFILEVALIDATOR_H_*/