#pragma once

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "base/core.h"
#include "base/flow.h"

namespace lms {
namespace flow {

using lms::core::ModelMeta;
using lms::flow::Flow;

template <typename T>
class FlowBuilder {
 public:
  FlowBuilder& from_directory(const std::string& dir) {
    dir_ = dir;
    return *this;
  }

  void BuildAll(const std::string& flow_file, const ModelMeta<T>& meta) {
    BuildFlowCsv(flow_file, meta.time_interval_s_);
  }

  const Flow<T>& flow() const { return flow_; }

 private:
  std::filesystem::path dir_;
  Flow<T> flow_;

  void BuildFlowCsv(const std::string& file_name, std::size_t time_interval_s) {
    std::ifstream file((dir_ / file_name).string());
    if (!file.is_open()) {
      throw std::runtime_error("FlowBuilder: Flow data file not found: " + file_name);
    }

    std::string line;
    if (!std::getline(file, line)) {
      return;
    }

    std::stringstream header_ss(line);
    std::string cell;
    std::vector<std::string> headers;
    while (std::getline(header_ss, cell, ',')) {
      headers.push_back(cell);
    }

    int q_col_idx = -1;
    for (std::size_t i = 0; i < headers.size(); ++i) {
      if (headers[i] == "Q" || headers[i] == "q") {
        q_col_idx = static_cast<int>(i);
        break;
      }
    }

    if (q_col_idx == -1) {
      throw std::runtime_error("FlowBuilder: Flow column 'Q' or 'q' not found in file: " + file_name);
    }

    flow_.duration = 0;
    flow_.time_interval_s = time_interval_s;
    flow_.data.clear();

    while (std::getline(file, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line);
      std::string value;
      int col = 0;
      T val = 0;
      bool found = false;
      while (std::getline(ss, value, ',')) {
        if (col == q_col_idx) {
          val = static_cast<T>(std::stod(value));
          found = true;
          break;
        }
        ++col;
      }
      if (found) {
        flow_.data.push_back(val);
        ++flow_.duration;
      }
    }
  }
};

}  // namespace flow
}  // namespace lms