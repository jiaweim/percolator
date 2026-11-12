#include "ValidateTabFile.h"

/**
 * ⭐2026-11-12
 * 合并多个 PIN 文件：以列数最多的 header 为标准 schema，逐个文件读取，对确实的列补 0，合并写入一个临时文件。
 * @param fileNames 多个 PIN 文件
 * @return 合并后临时文件的路径
 */
std::string ValidateTabFile::concatenateMultiplePINs(std::vector<std::basic_string<char>> fileNames) {
  // Validate multiple tab-files, and then concatenate them.
  // 最终返回的临时文件路径
  std::string inputFN;
  std::ifstream pinFileStream;

  // 标准 schema，记录 columnIndex->columnName
  std::map<int,std::string> columnMap = {};
  // 临时 map，每个文件解析 header 时用它暂存，然后与 columnMap 比大小
  std::map<int,std::string> tmpMap;

  /* Keep track on col index */
  int column_index;

  // 1. 遍历每个文件，找出最完整的 header
  for(const std::string &text : fileNames) {
      
    // 打开文件，读第一行
    pinFileStream.open(text.c_str(), std::ios::in);
    std::string headerRow;
    getline(pinFileStream, headerRow);

    // 记录 columnIndex->columnName
    TabReader reader(headerRow);
    column_index = 0;
    while (!reader.error()) {
      std::string optionalHeader = reader.readString();

      tmpMap[column_index] = optionalHeader;
      column_index++;
    }

    // 如果当前文件的 header 列数更多，就用 tmpMap 替换 columnMap
    if (columnMap.size() < tmpMap.size()) {
      columnMap = tmpMap;
    }
    tmpMap.clear();
    pinFileStream.close();
  }

  // 2. 创建临时文件
  string tcf = "";
  string tcd = "";

  TmpDir tmpDir;
  tmpDir.createTempFile(tcf, tcd); // 创建临时文件，tcf 为文件路径，tcd 是目录路径
  inputFN = tcf; // 记录路径，最后返回
  ofstream outFile;
  outFile.open(inputFN); // 打开输出流

  // 第一个文件负责写 header，后续文件跳过 header 和 defaultDirection
  bool firstFile = true;
  // 内层 while 的控制变量
  /* If go to next column when going though missing charge states */
  bool nextColumn;

  // 核心假设，当前文件的 header 是 columnMap 的子序列，即所有文件内容差不多，仅 charge 有的多有的少
  /* Start appending pin files to one pin file. */
  for(const string &text : fileNames) {
    if (VERB > 0) {
      std::cerr << "Reading file: " << text.c_str() << std::endl;
    }
    // 当前文件确实的 columnIndex，每个文件独立计算，因为不同文件的 header 可能不同
    std::vector<int> missingCols;
        
    // 打开文件，读 header
    pinFileStream.open(text.c_str(), std::ios::in);
    std::string headerRow;
    getline(pinFileStream, headerRow);
    TabReader reader(headerRow);

    // header 对齐与缺失列识别
    /* Keep track on column index */
    column_index = 0;
    /* Go through all columns to print header and search for missing charge state */
    while (!reader.error()) {
      // 遍历当前文件 header 的每个字段
      std::string optionalHeader = reader.readString();
      
      // 逐个扫描，与 columnMap 比对
      /* Check if charge state is missing */
      // 这里默认缺失的是 charge
      if (columnMap[column_index] != optionalHeader) { // 当前文件的第 column_index 字段与 columnMap 的 column_index 不匹配
        /* Charge state is missing at  'column_index' */
        missingCols.push_back(column_index);
        nextColumn=true;
        if (firstFile) {
          outFile << columnMap[column_index] << "\t";
        }
        /* Check if there are any adjacent missing charge states */
        while(nextColumn) {
          column_index++;
          if (columnMap[column_index] == optionalHeader) {
            /* No missing charge state, continue */
            nextColumn=false;
          } else {
            /* Missing charge state */
            missingCols.push_back(column_index);
            if (firstFile) {
              outFile << columnMap[column_index] << "\t";
            }             
          }
        }
      } 
      if (firstFile) {
        outFile << columnMap[column_index] << "\t";
      }
      column_index++;
    }

    /* Go through rest of the rows in pin-file */
    std::string nextRow;
    getline(pinFileStream, nextRow);

    if (firstFile) {
      outFile <<  "\n";
    } else {
      std::string defaultDirectionString = "defaultdirection";
      if (nextRow.size() >= defaultDirectionString.size()) {
        std::string psmid = nextRow.substr(0, defaultDirectionString.size());
        std::transform(psmid.begin(), psmid.end(), psmid.begin(), ::tolower);
        if (psmid == defaultDirectionString) {
          getline(pinFileStream, nextRow); // skip row with default direction
        }
      }
    }

    do {
      TabReader readerRow(nextRow);

      /* Keep track on column for missing charge states */
      int col = 0;
      while (!readerRow.error()) {           
        std::string value = readerRow.readString();
        /* Check for missing charge state */
        if (std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
          nextColumn=true;
          outFile << 0 << "\t";
          col++;
          /* Check for adjecent missing charge states */
          while(nextColumn) {
            if(std::find(missingCols.begin(), missingCols.end(), col) != missingCols.end()) {
              /* Missing adjecent charge state  */
              outFile << 0 << "\t";
              col++;
            } else {
              /* Charge states  */
              nextColumn = false;
            }
          }  
        }
        outFile << value << "\t";
        col++;
      }
      outFile << "\n";
    } while (getline(pinFileStream, nextRow));
    pinFileStream.close();
    firstFile = false;      
  }
  outFile.close();
  return inputFN;
}
